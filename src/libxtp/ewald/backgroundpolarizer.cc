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
#include <fstream>
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

  XTP_LOG(Log::error, _log) << _unit_cell << std::endl;
  XTP_LOG(Log::error, _log) << "Setup real and reciprocal space" << std::endl;
  RSpace rspace(_options, _unit_cell, ewaldSegments, _log);
  KSpace kspace(_options, _unit_cell, ewaldSegments, _log);

  XTP_LOG(Log::error, _log)
      << "Compute real space permanent fields" << std::endl;
  rspace.computeStaticField();
  XTP_LOG(Log::error, _log)
      << "Compute reciprocal space permanent fields" << std::endl;
  // kspace.computeStaticField();
  // kspace.computeShapeField();
  // kspace.computeIntraMolecularCorrection();

  std::ofstream infile2;
  infile2.open("staticFieldXTP.txt"); 
  infile2 << "id x y z q Ex Ey Ez" << std::endl;
  for (const auto& seg : ewaldSegments) {
    for (const auto& site : seg) {
      infile2 << site << std::endl;
    }
  }

  exit(0);

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
  std::cout << "Computed initial guess" << std::endl;

  // Set up the dipole interaction matrix
  Eigen::MatrixXd inducedDipoleInteraction =
      rspace.getInducedDipoleInteraction();
  std::cout << "Setup real space part of dipole interaction matrix"
            << std::endl;
  inducedDipoleInteraction += kspace.getInducedDipoleInteraction();
  std::cout << "Setup reciprocal space part of dipole interaction matrix"
            << std::endl;

  // Add  the inverse polarization
  Index diagIndex = 0;
  for (auto& seg : ewaldSegments) {
    for (auto& site : seg) {
      inducedDipoleInteraction.block<3, 3>(diagIndex, diagIndex) +=
          site.getPolarizationMatrix().inverse();
      diagIndex += 3;
    }
  }
  std::cout << "Setup the inverse polarization matrix. \nStarting "
               "preconditioned conjugate gradient solver"
            << std::endl;

  // Solving the linear system
  Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg;
  cg.setMaxIterations(100);
  cg.setTolerance(1e-3);
  cg.compute(inducedDipoleInteraction);
  Eigen::VectorXd x = cg.solveWithGuess(staticField, initialGuess);

  std::cout << TimeStamp() << " CG: #iterations: " << cg.iterations()
            << ", estimated error: " << cg.error() << std::endl;

  if (cg.info() == Eigen::ComputationInfo::NoConvergence) {
    std::cout << "PCG iterations did not converge" << std::endl;
  }

  if (cg.info() == Eigen::ComputationInfo::NumericalIssue) {
    std::cout << "PCG had a numerical issue" << std::endl;
  }

  x = 0.05291 * x;  // convert to CTP units

  std::ofstream outfile;
  outfile.open("inducedDipoleXTP.txt");
  outfile << "idx idy idz" << std::endl;
  for (Index i = 0; i < systemSize; i += 3) {
    outfile << x[i + 0] << " " << x[i + 1] << " " << x[i + 2] << "\n";
  }
  outfile << std::endl;
}
}  // namespace xtp
}  // namespace votca
