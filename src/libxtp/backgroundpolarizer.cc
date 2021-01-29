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

void BackgroundPolarizer::computeStaticFieldAt(
    Index segId, std::vector<PolarSegment>& polarSegments) {

  // compute realspace sum

  // compute reciprocal space sum

  // self interaction correction

  // aperiodic correction

  // shape contribution
}

void BackgroundPolarizer::computeNeighbourList(
    std::vector<PolarSegment>& polarSegments) {
  std::array<Index, 3> maxCopies =
      _unit_cell.getNrOfRealSpaceCopiesForCutOff(_options.realcutoff);
  _nbList.setSize(polarSegments.size());

#pragma omp parallel for
  for (Index segId = 0; segId < polarSegments.size(); ++segId) {
    PolarSegment& currentSeg = polarSegments[segId];
    for (PolarSegment seg : polarSegments) {
      Eigen::Vector3d dr = _unit_cell.minImage(currentSeg, seg);
      // triple for-loop is over all unitcell copies
      for (Index n1 = -maxCopies[0]; n1 < maxCopies[0]; ++n1) {
        for (Index n2 = -maxCopies[1]; n2 < maxCopies[1]; ++n2) {
          for (Index n3 = -maxCopies[2]; n3 < maxCopies[2]; ++n3) {
            if (n1 == 0 && n2 == 0 && n3 == 0 &&
                currentSeg.getId() == seg.getId()) {
              continue;
            }
            // LVector is the vector pointing to the n1,n2,n3th box
            Eigen::Vector3d lvector = _unit_cell.getLVector(n1, n2, n3);
            Eigen::Vector3d dr_l = dr + lvector;
            double dist = dr_l.norm();
            if (dist < _options.realcutoff) {
              _nbList.addNeighbourTo(segId, Neighbour(seg.getId(), dr_l, dist));
            }
          }
        }
      }
    }
    _nbList.sortOnDistance(segId);
  }
}

void BackgroundPolarizer::computeStaticFields(
    std::vector<PolarSegment>& polarSegments) {
  std::cout << "Hey Hallo" << std::endl;
}

void BackgroundPolarizer::Polarize(std::vector<PolarSegment>& polarSegments) {

  computeNeighbourList(polarSegments);

  computeStaticFields(polarSegments);

  // computeInducedFields(polarSegments);
}
}  // namespace xtp
}  // namespace votca
