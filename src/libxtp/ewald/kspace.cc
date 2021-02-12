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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <vector>

// Local includes
#include "kspace.h"

namespace votca {
namespace xtp {

KSpace::KSpace(const EwaldOptions& options, const UnitCell& unitcell,
               std::vector<EwdSegment>& ewaldSegments)
    : _unit_cell(unitcell), _ewaldSegments(ewaldSegments) {
  a1 = options.alpha;
  a2 = a1 * a1;
  fourPiVolume =
      4.0 * boost::math::constants::pi<double>() / _unit_cell.getVolume();
  cutoff = options.k_cutoff;
  cutoff2 = cutoff * cutoff;

  // compute max k-space vectors
  const Eigen::Matrix3d& inverseCellMatrix = _unit_cell.getInverseMatrix();
  for (Index i = 0; i < 3; ++i) {
    max_K[i] =
        static_cast<Index>(std::ceil(cutoff / inverseCellMatrix.col(i).norm()));
  }
  std::cout << "\nMax_K: [" << max_K[0] << ", " << max_K[1] << ", " << max_K[2]
            << "]\n"
            << std::endl;
  std::cout << std::endl;

  for (const EwdSegment& seg : _ewaldSegments) {
    for (const EwdSite& sit : seg) {
      mu_tot += sit.getStaticDipole();
    }
  }
}

void KSpace::computeStaticField() {
  computeKVectors();
  
#pragma omp parallel for
  for (Index i = 0; i < Index(_ewaldSegments.size()); ++i) {
    EwdSegment& seg = _ewaldSegments[i];
    for (EwdSite& site : seg) {
      for (const KVector& kvec : _kvector_list) {

        site.addToStaticField(
          (
            fourPiVolume *  ii *kvec.getVector() *std::exp(ii*kvec.getVector().dot(site.getPos())) * kvec.getAk() * std::conj(kvec.getSk()) ).real());
      }
    }
  }
}

void KSpace::testFunction(){;
}

/**************************************************
 * PRIVATE FUNCTIONS                              *
 **************************************************/

std::complex<double> KSpace::computeSk(const Eigen::Vector3d& kvector) const {
  std::complex<double> sk(0.0, 0.0);
  for (const EwdSegment& seg : _ewaldSegments) {
    for (const EwdSite& site : seg) {
      std::complex<double> generalizedCharge =
          (site.getCharge() + ii * kvector.dot(site.getTotalDipole()) -
           kvector.dot(site.getQuadrupole() * kvector));
      std::complex<double> expKR = std::exp(ii * kvector.dot(site.getPos()));
      sk += generalizedCharge * expKR;
    }
  }
  return sk;
}

double KSpace::computeAk(const Eigen::Vector3d& kvector) const {
  /* Compute the A_k factor, i.e. k^(-2) exp(-k^2/(4\alpha^2)) */
  double k_squared = kvector.squaredNorm();
  return std::exp(-k_squared / (4 * a2)) / k_squared;
}

void KSpace::computeKVectors() {
  for (Index ix = -max_K[0]; ix <= max_K[0]; ++ix) {
    for (Index iy = -max_K[1]; iy <= max_K[1]; ++iy) {
      for (Index iz = -max_K[2]; iz <= max_K[2]; ++iz) {
        if (ix == 0 && iy == 0 && iz == 0) {
          continue;
        }
        Eigen::Vector3d kvector = _unit_cell.getKVector(ix, iy, iz);
        _kvector_list.push_back(KVector(kvector));
      }
    }
  }
  std::sort(_kvector_list.begin(), _kvector_list.end());
  std::cout << "Part 2" << std::endl;
#pragma omp parallel for
  for (Index i = 0; i < static_cast<Index>(_kvector_list.size()); ++i) {
    KVector& kvec = _kvector_list[i];
    kvec.setSk(computeSk(kvec.getVector()));
    kvec.setAk(computeAk(kvec.getVector()));
  }
}

}  // namespace xtp
}  // namespace votca