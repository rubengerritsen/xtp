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

#pragma once
#ifndef VOTCA_XTP_EWALDKSINTERACTOR_H
#define VOTCA_XTP_EWALDKSINTERACTOR_H

#include <boost/math/constants/constants.hpp>
#include <complex>
#include <vector>

// Local VOTCA includes

#include "votca/xtp/ewaldsegment.h"
#include "votca/xtp/ewaldunitcell.h"

namespace votca {
namespace xtp {

class EwaldKSpace {
 public:
  EwaldKSpace(double alpha, const EwaldUnitCell& unitcell,
              std::vector<EwaldSegment>& ewaldSegments)
      : alpha(alpha), _unit_cell(unitcell), ewaldSegments(ewaldSegments) {
    a1 = alpha;
    a2 = alpha * alpha;
  };
  ~EwaldKSpace() = default;

  std::complex<double> computeSk(const Eigen::Vector3d& kvector) const;
  double computeAk(const Eigen::Vector3d& kvector) const;

 private:
  EwaldUnitCell _unit_cell;
  std::vector<EwaldSegment> ewaldSegments;
  double alpha, a1, a2;  // alpha (splitting param) and its powers
  static constexpr std::complex<double> ii = std::complex<double>(0.0, 1.0);  // imaginary i
};
}  // namespace xtp
}  // namespace votca

#endif