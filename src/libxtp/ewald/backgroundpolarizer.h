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
#ifndef VOTCA_XTP_BACKGROUNDPOLARIZER_H
#define VOTCA_XTP_BACKGROUNDPOLARIZER_H
#include <vector>
// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/logger.h"

// Private VOTCA includes
#include "ewd_segment.h"
#include "unitcell.h"

namespace votca {
namespace xtp {

struct EwaldOptions {
  double k_cutoff = 0.4;
  double r_cutoff = 134.1;
  double alpha = 6 / 18.897259;
  double sharpness = 0.390;
};

class BackgroundPolarizer {
 public:
  BackgroundPolarizer(Logger& log, UnitCell& unitcell, EwaldOptions options)
      : _log(log), _unit_cell(unitcell), _options(options){};

  ~BackgroundPolarizer() = default;

  void Polarize(std::vector<EwdSegment>& ewaldSegments);

 private:
  Index computeSystemSize(std::vector<EwdSegment>& ewaldSegments) const;
  Logger& _log;
  UnitCell _unit_cell;
  EwaldOptions _options;
};
}  // namespace xtp
}  // namespace votca

#endif