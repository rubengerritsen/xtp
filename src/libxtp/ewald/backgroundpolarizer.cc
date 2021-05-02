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

#include "backgroundpolarizer.h"
#include "kspace.h"
#include "rspace.h"
#include <iostream>
#include <vector>

namespace votca {
namespace xtp {

Index BackgroundPolarizer::computeSystemSize(
    std::vector<EwdSegment>& ewaldSegments) const {
  Index systemSize = 0;
  for (const auto& seg : ewaldSegments) {
    systemSize += 3 * seg.size();
  }
  return systemSize;
}

void BackgroundPolarizer::Polarize(std::vector<EwdSegment>& ewaldSegments) {

  Index systemSize = computeSystemSize(ewaldSegments);

  // Setup the real and reciprocal space
  RSpace rspace(_options, _unit_cell, ewaldSegments);
  KSpace kspace(_options, _unit_cell, ewaldSegments);

  // Generate Permanent fields (F^a_alpha)
  rspace.computeStaticField();
  kspace.computeStaticField();
  kspace.computeShapeField(Shape::cube);
  kspace.computeIntraMolecularCorrection();

  /*******************************************************
   * We turn the problem into a linear problem (A + B)x = b
   * b = the permanent field
   * x = are the induced dipoles
   * A = the "interaction" tensor that turns x into corresponding induced field
   *     at that location
   * B = The inverse polarization matrix (a super matrix containing on it's
   *     diagonal the 3x3 inverse polarization matrices)
   * *****************************************************/

  // Get static field from the sites and convert it into a big 1D vector
  // The  same for the initial guess
  Eigen::VectorXd staticField = Eigen::VectorXd::Zero(systemSize);
  Eigen::VectorXd initialGuess = Eigen::VectorXd::Zero(systemSize);
  Index index = 0;
  for (auto& seg : ewaldSegments) {
    for (auto& site : seg) {
      site.induceDirect();  // compute induced dipole based on static field
      Eigen::Vector3d E = site.getStaticField();
      Eigen::Vector3d induced_dipole = site.getInducedDipole();
      staticField.segment<3>(index) = E;
      initialGuess.segment<3>(index) = induced_dipole;
      index += 3;
    }
  }

  // Set up the dipole interaction matrix
  Eigen::MatrixXd inducedDipoleInteraction =
      rspace.getInducedDipoleInteraction();
  inducedDipoleInteraction += kspace.getInducedDipoleInteraction();

  // Calculate induction field
  rspace.computeIntraMolecularField();
  rspace.computeInducedField();

  kspace.computeInducedField();
  kspace.computeInducedShapeField();
  kspace.computeIntraMolecularInducedCorrection();

  // // Update induced dipoles with SOR
  // double wSOR = 0.35;
  // for (const EwdSegment& seg : ewaldSegements) {
  //   for (const EwdSite& site : seg) {
  //     site.updateInducedDipoles(wSOR);
  //   }
  // }

  // Convergence check
}
}  // namespace xtp
}  // namespace votca
