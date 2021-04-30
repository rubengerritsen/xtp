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
  a3 = a1 * a2;
  a4 = a2 * a2;
  a5 = a4 * a1;
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
}

void KSpace::computeStaticField() {
  computeKVectors();

#pragma omp parallel for
  for (Index i = 0; i < Index(_ewaldSegments.size()); ++i) {
    EwdSegment& seg = _ewaldSegments[i];
    for (EwdSite& site : seg) {
      for (const KVector& kvec : _kvector_list) {

        site.addToStaticField(
            (fourPiVolume * ii * kvec.getVector() *
             std::exp(ii * kvec.getVector().dot(site.getPos())) * kvec.getAk() *
             std::conj(kvec.getSk()))
                .real());
      }
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
            (fourPiVolume * ii * kvec.getVector() *
             std::exp(ii * kvec.getVector().dot(site.getPos())) * kvec.getAk() *
             std::conj(kvec.getSk()))
                .real());
      }
    }
  }
}

void KSpace::computeShapeField(Shape shape) {
  Eigen::Vector3d dip_tot = Eigen::Vector3d::Zero();
  Eigen::Vector3d shapeField = Eigen::Vector3d::Zero();

  // Compute total dipole
  for (const EwdSegment& seg : _ewaldSegments) {
    for (const EwdSite& sit : seg) {
      dip_tot += sit.getCharge() * sit.getPos();
      dip_tot += sit.getStaticDipole();
    }
  }

  switch (shape) {
    case Shape::xyslab:
      shapeField[2] = fourPiVolume * dip_tot[2];
      break;
    case Shape::cube:
    case Shape::sphere:
      shapeField = (fourPiVolume / 3.0) * dip_tot;
      break;
    default:
      throw std::runtime_error("Shape not implemented.");
  }

  // Apply the field to the sites
  for (EwdSegment& seg : _ewaldSegments) {
    for (EwdSite& sit : seg) {
      sit.addToStaticField(shapeField);
    }
  }
}

void KSpace::computeInducedShapeField(Shape shape) {
  Eigen::Vector3d dip_tot = Eigen::Vector3d::Zero();
  Eigen::Vector3d shapeField = Eigen::Vector3d::Zero();

  // Compute total dipole
  for (const EwdSegment& seg : _ewaldSegments) {
    for (const EwdSite& sit : seg) {
      dip_tot += sit.getInducedDipole();
    }
  }

  switch (shape) {
    case Shape::xyslab:
      shapeField[2] = fourPiVolume * dip_tot[2];
      break;
    case Shape::cube:
    case Shape::sphere:
      shapeField = (fourPiVolume / 3.0) * dip_tot;
      break;
    default:
      throw std::runtime_error("Shape not implemented.");
  }

  // Apply the field to the sites
  for (EwdSegment& seg : _ewaldSegments) {
    for (EwdSite& sit : seg) {
      sit.addToInducedField(shapeField);
    }
  }
}

void KSpace::computeIntraMolecularCorrection() {
  for (EwdSegment& seg : _ewaldSegments) {
    for (EwdSite& sit1 : seg) {
      for (EwdSite& sit2 : seg) {
        sit1.addToStaticField(-staticFieldAtBy(sit1, sit2));
      }
    }
  }
}

void KSpace::computeIntraMolecularInducedCorrection() {
  for (EwdSegment& seg : _ewaldSegments) {
    for (EwdSite& sit1 : seg) {
      sit1.addToInducedField(-inducedFieldAtBy(sit1, sit1));
    }

  }
}

void KSpace::testFunction() { ; }

/**************************************************
 * PRIVATE FUNCTIONS                              *
 **************************************************/

void KSpace::computeScreenedInteraction() {
  double rSqrtPiExp = rSqrtPi * std::exp(-a2 * R2);
  // Note KSpace screening is with erf
  rR1s = std::erf(a1 * R1) * rR1;
  rR3s = rR2 * (rR1s - 2.0 * a1 * rSqrtPiExp);
  rR5s = rR2 * (3.0 * rR3s - (4.0 * a3) * rSqrtPiExp);
  rR7s = rR2 * (5.0 * rR5s - (8.0 * a5) * rSqrtPiExp);
}

void KSpace::computeDistanceVariables(Eigen::Vector3d distVec) {
  dr = distVec;
  R1 = dr.norm();
  R2 = R1 * R1;
  rR1 = 1.0 / R1;
  rR2 = rR1 * rR1;
}

Eigen::Vector3d KSpace::staticFieldAtBy(EwdSite& site, const EwdSite& nbSite) {
  computeDistanceVariables(_unit_cell.minImage(site, nbSite));
  computeScreenedInteraction();

  Eigen::Vector3d field = Eigen::Vector3d::Zero();
  Index rank = nbSite.getRank();
  if (R1 < 1e-2) {
    field += 4.0 / 3.0 * a3 * rSqrtPi * nbSite.getStaticDipole();
  } else {
    // charge
    field += -nbSite.getCharge() * dr * rR3s;
    if (rank > 0) {  // dipole
      field += nbSite.getStaticDipole() * rR3s;
      field += -rR5s * dr * dr.dot(nbSite.getStaticDipole());
      if (rank > 1) {  // quadrupole
        // Using that the trace of a quadrupole contributes nothing, we can skip
        // that part
        field += rR5s * 2 * nbSite.getQuadrupole() * dr;
        field += -rR7s * dr * dr.dot(nbSite.getQuadrupole() * dr);
      }
    }
  }
  return field;
}

Eigen::Vector3d KSpace::inducedFieldAtBy(EwdSite& site, const EwdSite& nbSite) {
  computeDistanceVariables(_unit_cell.minImage(site, nbSite));
  computeScreenedInteraction();

  Eigen::Vector3d field = Eigen::Vector3d::Zero();
  Index rank = nbSite.getRank();
  if (R1 < 1e-2) {
    field += 4.0 / 3.0 * a3 * rSqrtPi * nbSite.getInducedDipole();
  } else {
    field += nbSite.getInducedDipole() * rR3s;
    field += -rR5s * dr * dr.dot(nbSite.getInducedDipole());
  }
  return field;
}

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
#pragma omp parallel for
  for (Index i = 0; i < static_cast<Index>(_kvector_list.size()); ++i) {
    KVector& kvec = _kvector_list[i];
    kvec.setSk(computeSk(kvec.getVector()));
    kvec.setAk(computeAk(kvec.getVector()));
  }
}

}  // namespace xtp
}  // namespace votca