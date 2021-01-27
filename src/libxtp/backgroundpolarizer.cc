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

#include "votca/xtp/backgroundpolarizer.h"
#include <vector>

namespace votca {
namespace xtp {

void BackgroundPolarizer::computeStaticFieldAt(Index segId, std::vector<PolarSegment>& polarSegments){

  // compute realspace sum

  // compute reciprocal space sum

  // self interaction correction

  // aperiodic correction

  // shape contribution
}

void BackgroundPolarizer::computeStaticFields(std::vector<PolarSegment>& polarSegments){
  // Create Neighbours
  


  // Perform actual field computations
  for (Index segment = 0; segment < polarSegments.size(); ++segment){
    computeStaticFieldAt(segment, polarSegments);
  }
}


void BackgroundPolarizer::Polarize(std::vector<PolarSegment>& polarSegments) {

  computeStaticFields(polarSegments);

  //computeInducedFields(polarSegments);

}
}  // namespace xtp
}  // namespace votca
