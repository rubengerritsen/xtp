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

void BackgroundPolarizer::Polarize(std::vector<EwdSegment>& ewaldSegments) {

  // Setup the real and reciprocal space
  RSpace rspace(_options, _unit_cell, ewaldSegments);
  KSpace kspace(_options, _unit_cell, ewaldSegments);

  std::cout << "Starting RSpace part" << std::endl;
  //rspace.computeStaticField();
  std::cout << "Starting KSpace part" << std::endl;
  //kspace.computeStaticField();
  std::cout << "Starting Shape Calculation" << std::endl;
  kspace.computeShapeField(Shape::cube);

  std::cout << " ID " << ewaldSegments[1100].getId() << std::endl;
  for (auto& site : ewaldSegments[1100]) {
    std::cout << site << std::endl;
  }
}
}  // namespace xtp
}  // namespace votca
