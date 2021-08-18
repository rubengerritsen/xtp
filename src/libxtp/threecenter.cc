/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/threecenter.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/make_libint_work.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/symmetric_matrix.h"

// include libint last otherwise it overrides eigen
#include <libint2.hpp>

namespace votca {
namespace xtp {

void TCMatrix_dft::Fill(const AOBasis& auxbasis, const AOBasis& dftbasis) {
  {
    AOCoulomb auxAOcoulomb;
    auxAOcoulomb.Fill(auxbasis);
    inv_sqrt_ = auxAOcoulomb.Pseudo_InvSqrt(1e-8);
    removedfunctions_ = auxAOcoulomb.Removedfunctions();
  }
  matrix_ = std::vector<Symmetric_Matrix>(auxbasis.AOBasisSize());

#pragma omp parallel for schedule(dynamic, 4)
  for (Index i = 0; i < auxbasis.AOBasisSize(); i++) {
    matrix_[i] = Symmetric_Matrix(dftbasis.AOBasisSize());
  }

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> dftshells = dftbasis.GenerateLibintBasis();
  std::vector<libint2::Shell> auxshells = auxbasis.GenerateLibintBasis();
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(
      libint2::Operator::coulomb,
      std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
      static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 0);
  engines[0].set(libint2::BraKet::xs_xx);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::vector<Index> shell2bf = dftbasis.getMapToBasisFunctions();
  std::vector<Index> auxshell2bf = auxbasis.getMapToBasisFunctions();

#pragma omp parallel for schedule(dynamic)
  for (Index is = dftbasis.getNumofShells() - 1; is >= 0; is--) {

    libint2::Engine& engine = engines[OPENMP::getThreadId()];
    const libint2::Engine::target_ptr_vec& buf = engine.results();
    const libint2::Shell& dftshell = dftshells[is];
    Index start = shell2bf[is];
    std::vector<Eigen::MatrixXd> block(dftshell.size());
    for (Index i = 0; i < Index(dftshell.size()); i++) {
      Index size = start + i + 1;
      block[i] = Eigen::MatrixXd::Zero(auxbasis.AOBasisSize(), size);
    }

    for (Index aux = 0; aux < auxbasis.getNumofShells(); aux++) {
      const libint2::Shell& auxshell = auxshells[aux];
      Index aux_start = auxshell2bf[aux];

      for (Index dis = 0; dis <= is; dis++) {

        const libint2::Shell& shell_col = dftshells[dis];
        Index col_start = shell2bf[dis];
        engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 0>(
            auxshell, libint2::Shell::unit(), dftshell, shell_col);

        if (buf[0] == nullptr) {
          continue;
        }
        Eigen::TensorMap<Eigen::Tensor<const double, 3, Eigen::RowMajor> const>
            result(buf[0], auxshell.size(), dftshell.size(), shell_col.size());

        for (size_t left = 0; left < dftshell.size(); left++) {
          for (size_t auxf = 0; auxf < auxshell.size(); auxf++) {
            for (size_t col = 0; col < shell_col.size(); col++) {
              // symmetry
              if ((col_start + col) > (start + left)) {
                break;
              }
              block[left](aux_start + auxf, col_start + col) =
                  result(auxf, left, col);
            }
          }
        }
      }
    }

    for (Index i = 0; i < Index(block.size()); ++i) {
      Eigen::MatrixXd temp = inv_sqrt_ * block[i];
      for (Index mu = 0; mu < temp.rows(); ++mu) {
        for (Index j = 0; j < temp.cols(); ++j) {
          matrix_[mu](i + start, j) = temp(mu, j);
        }
      }
    }
  }

  return;
}

void TCMatrix_gwbse::Initialize(Index basissize, Index mmin, Index mmax,
                                Index nmin, Index nmax) {

  // here as storage indices starting from zero
  nmin_ = nmin;
  nmax_ = nmax;
  ntotal_ = nmax - nmin + 1;
  mmin_ = mmin;
  mmax_ = mmax;
  mtotal_ = mmax - mmin + 1;
  auxbasissize_ = basissize;

  // vector has mtotal elements
  // largest object should be allocated in multithread fashion
  matrix_ = std::vector<Eigen::MatrixXd>(mtotal_);
#pragma omp parallel for schedule(dynamic, 4)
  for (Index i = 0; i < mtotal_; i++) {
    matrix_[i] = Eigen::MatrixXd::Zero(ntotal_, auxbasissize_);
  }
}

/*
 * Modify 3-center matrix elements consistent with use of symmetrized
 * Coulomb interaction using either CUDA or Openmp.
 */
void TCMatrix_gwbse::MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix) {
  OpenMP_CUDA gemm;
  gemm.setOperators(matrix_, matrix);
#pragma omp parallel
  {
    Index threadid = OPENMP::getThreadId();
#pragma omp for schedule(dynamic)
    for (Index i = 0; i < msize(); i++) {
      gemm.MultiplyRight(matrix_[i], threadid);
    }
  }
}
/*
 * Fill the 3-center object by looping over shells of GW basis set and
 * calling FillBlock, which calculates all 3-center overlap integrals
 * associated to a particular shell, convoluted with the DFT orbital
 * coefficients
 */
void TCMatrix_gwbse::Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
                          const Eigen::MatrixXd& dft_orbitals) {
  // needed for Rebuild())
  auxbasis_ = &auxbasis;
  dftbasis_ = &dftbasis;
  dft_orbitals_ = &dft_orbitals;

  Fill3cMO(auxbasis, dftbasis, dft_orbitals);

  AOOverlap auxoverlap;
  auxoverlap.Fill(auxbasis);
  AOCoulomb auxcoulomb;
  auxcoulomb.Fill(auxbasis);
  Eigen::MatrixXd inv_sqrt = auxcoulomb.Pseudo_InvSqrt_GWBSE(auxoverlap, 5e-7);
  removedfunctions_ = auxcoulomb.Removedfunctions();
  MultiplyRightWithAuxMatrix(inv_sqrt);

  return;
}

/*
 * Determines the 3-center integrals for a given shell in the aux basis
 * by calculating the 3-center repulsion integral of the functions in the
 * aux shell with ALL functions in the DFT basis set
 */
std::vector<Eigen::MatrixXd> TCMatrix_gwbse::ComputeAO3cBlock(
    const libint2::Shell& auxshell, const AOBasis& dftbasis,
    libint2::Engine& engine) const {
  std::vector<Eigen::MatrixXd> ao3c = std::vector<Eigen::MatrixXd>(
      auxshell.size(),
      Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize()));

  std::vector<libint2::Shell> dftshells = dftbasis.GenerateLibintBasis();
  std::vector<Index> shell2bf = dftbasis.getMapToBasisFunctions();

  const libint2::Engine::target_ptr_vec& buf = engine.results();
  // alpha-loop over the "left" DFT basis function
  for (Index row = 0; row < Index(dftshells.size()); row++) {

    const libint2::Shell& shell_row = dftshells[row];
    const Index row_start = shell2bf[row];
    // ThreecMatrix is symmetric, restrict explicit calculation to triangular
    // matrix
    for (Index col = 0; col <= row; col++) {
      const libint2::Shell& shell_col = dftshells[col];
      const Index col_start = shell2bf[col];

      engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 0>(
          auxshell, libint2::Shell::unit(), shell_col, shell_row);

      if (buf[0] == nullptr) {
        continue;
      }
      Eigen::TensorMap<Eigen::Tensor<const double, 3, Eigen::RowMajor> const>
          result(buf[0], auxshell.size(), shell_col.size(), shell_row.size());

      for (size_t aux_c = 0; aux_c < auxshell.size(); aux_c++) {
        for (size_t col_c = 0; col_c < shell_col.size(); col_c++) {
          for (size_t row_c = 0; row_c < shell_row.size(); row_c++) {
            // ao3c is col major result is row-major so this is ideal
            ao3c[aux_c](row_start + row_c, col_start + col_c) =
                result(aux_c, col_c, row_c);
          }  // ROW copy
        }    // COL copy
      }      // AUX copy

    }  // gamma-loop
  }    // alpha-loop

  for (Eigen::MatrixXd& mat : ao3c) {
    mat.triangularView<Eigen::Upper>() =
        mat.triangularView<Eigen::Lower>().transpose();
  }
  return ao3c;
}

void TCMatrix_gwbse::Fill3cMO(const AOBasis& auxbasis, const AOBasis& dftbasis,
                              const Eigen::MatrixXd& dft_orbitals) {

  const Eigen::MatrixXd dftm = dft_orbitals.middleCols(mmin_, mtotal_);
  const Eigen::MatrixXd dftn =
      dft_orbitals.middleCols(nmin_, ntotal_).transpose();

  OpenMP_CUDA transform;
  transform.setOperators(dftn, dftm);
  Index nthreads = OPENMP::getMaxThreads();

  std::vector<libint2::Shell> auxshells = auxbasis.GenerateLibintBasis();
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(
      libint2::Operator::coulomb,
      std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
      static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 0);
  engines[0].set(libint2::BraKet::xs_xx);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::vector<Index> auxshell2bf = auxbasis.getMapToBasisFunctions();

#pragma omp parallel
  {
    Index threadid = OPENMP::getThreadId();
#pragma omp for schedule(dynamic)
    for (Index aux = 0; aux < Index(auxshells.size()); aux++) {
      const libint2::Shell& auxshell = auxshells[aux];

      std::vector<Eigen::MatrixXd> ao3c =
          ComputeAO3cBlock(auxshell, dftbasis, engines[threadid]);

      // this is basically a transpose of AO3c and at the same time the ao->mo
      // transformation
      // we do not want to put it into  matrix_ straight away is because,
      //  matrix_ is shared between all threads and we want a nice clean access
      // pattern to it
      std::vector<Eigen::MatrixXd> block = std::vector<Eigen::MatrixXd>(
          mtotal_, Eigen::MatrixXd::Zero(ntotal_, ao3c.size()));

      Index dim = static_cast<Index>(ao3c.size());
      for (Index k = 0; k < dim; ++k) {
        transform.MultiplyLeftRight(ao3c[k], threadid);
        for (Index i = 0; i < ao3c[k].cols(); ++i) {
          block[i].col(k) = ao3c[k].col(i);
        }
      }

      // put into correct position
      for (Index m_level = 0; m_level < mtotal_; m_level++) {
        matrix_[m_level].middleCols(auxshell2bf[aux], auxshell.size()) =
            block[m_level];
      }  // m-th DFT orbital
    }    // shells of GW basis set
  }
}

}  // namespace xtp
}  // namespace votca
