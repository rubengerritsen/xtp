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
#include "votca/xtp/bgneighbourlist.h"
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/ewaldsegment.h"
#include "votca/xtp/ewaldunitcell.h"
#include "votca/xtp/logger.h"

namespace votca {
namespace xtp {

struct EwaldOptions {
  double realcutoff = 134.1;
  double ewaldsplitting = 1.0 / 18.9;
};

class BackgroundPolarizer {
 public:
  BackgroundPolarizer(Logger& log, EwaldUnitCell& unitcell,
                      EwaldOptions options)
      : _log(log), _unit_cell(unitcell), _options(options){};

  ~BackgroundPolarizer() = default;

  void Polarize(std::vector<EwaldSegment>& ewaldSegments);

 private:
  Logger& _log;

  BgNeighbourList _nbList;
  EwaldUnitCell _unit_cell;
  EwaldOptions _options;

  Index _rs_n1_max, _rs_n2_max, _rs_n3_max;

  void computeStaticFieldAt(Index segId,
                            std::vector<EwaldSegment>& ewaldSegments);
  void computeStaticFieldsRS(std::vector<EwaldSegment>& ewaldSegments);
  void computeNeighbourList(std::vector<EwaldSegment>& ewaldSegments);
};
}  // namespace xtp
}  // namespace votca

#endif