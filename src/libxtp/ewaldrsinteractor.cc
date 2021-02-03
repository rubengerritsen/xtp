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

#include "votca/xtp/ewaldrsinteractor.h"
#include <vector>

namespace votca {
namespace xtp {

void EwaldRSInteractor::updateVariables(Eigen::Vector3d distVec) {
  dr = distVec;
  R1 = dr.norm();
  R2 = R1 * R1;
  rR1 = 1.0 / R1;
  rR2 = rR1 * rR1;
}

void EwaldRSInteractor::computeScreenedInteraction() {
  rSqrtPiExp = rSqrtPi * std::exp(-a2 * R2);

  rR1s = std::erfc(a1 * R1) * rR1;
  rR3s = rR2 * (rR1s + 2.0 * a1 * rSqrtPiExp);
  rR5s = rR2 * (3.0 * rR3s + (4.0 * a3) * rSqrtPiExp);
  rR7s = rR2 * (5.0 * rR5s + (8.0 * a5) * rSqrtPiExp);
}

void EwaldRSInteractor::staticField(EwaldSite& site, const EwaldSite& nbSite,
                                    const Eigen::Vector3d shift) {
  updateVariables(_unit_cell.minImage(site, nbSite) + shift);
  computeScreenedInteraction();

  Eigen::Vector3d field = Eigen::Vector3d::Zero();
  Index rank = nbSite.getRank();
  // charge
  field += -nbSite.getCharge() * dr * rR3s;
  if (rank > 0) {  // dipole
    field += nbSite.getStaticDipole() * rR3s;
    field += -rR5s * dr * dr.dot(nbSite.getStaticDipole());
    if (rank > 1) {  // quadrupole
      // Using that the quadrupole is traceless we can skip that part
      field += rR5s * 2 * nbSite.getQuadrupole() * dr;
      Eigen::Matrix3d dyadic = dr * dr.transpose();
      field +=
          -rR7s * dr * (dyadic.array() * nbSite.getQuadrupole().array()).sum();
    }
  }

  site.addToStaticField(field);
}

}  // namespace xtp
}  // namespace votca