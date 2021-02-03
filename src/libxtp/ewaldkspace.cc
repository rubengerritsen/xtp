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

#include "votca/xtp/ewaldkspace.h"
#include <vector>

namespace votca {
namespace xtp {

std::complex<double> EwaldKSpace::computeSk(
    const Eigen::Vector3d& kvector) const {
  std::complex<double> sk(0.0, 0.0);
  for (const EwaldSegment& seg : ewaldSegments) {
    for (const EwaldSite& site : seg) {
      std::complex<double> generalizedCharge =
          (site.getCharge() + ii * kvector.dot(site.getTotalDipole()) -
           kvector.dot(site.getQuadrupole() * kvector));
      std::complex<double>  expKR = std::exp(ii * kvector.dot(site.getPos()));
      sk += generalizedCharge * expKR;
    }
  }
}

double EwaldKSpace::computeAk(const Eigen::Vector3d& kvector) const {
  /* Compute the A_k factor, i.e. k^(-2) exp(-k^2/(4\alpha^2)) */
  double k_squared = kvector.squaredNorm();
  return std::exp(-k_squared / (4 * a2)) / k_squared;
}

}  // namespace xtp
}  // namespace votca