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
#include "votca/xtp/ewaldkspace.h"
#include "votca/xtp/ewaldrsinteractor.h"
#include <vector>

namespace votca {
namespace xtp {

void BackgroundPolarizer::computeStaticFieldAt(
    Index segId, std::vector<EwaldSegment>& ewaldSegments) {

  // compute realspace sum

  // compute reciprocal space sum

  // self interaction correction

  // aperiodic correction

  // shape contribution
}

void BackgroundPolarizer::computeNeighbourList(
    std::vector<EwaldSegment>& ewaldSegments) {
  std::array<Index, 3> maxCopies =
      _unit_cell.getNrOfRealSpaceCopiesForCutOff(_options.realcutoff);
  _nbList.setSize(ewaldSegments.size());

#pragma omp parallel for
  for (Index segId = 0; segId < ewaldSegments.size(); ++segId) {
    EwaldSegment& currentSeg = ewaldSegments[segId];
    for (EwaldSegment seg : ewaldSegments) {
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
              _nbList.addNeighbourTo(
                  segId, Neighbour(seg.getId(), dr_l, lvector, dist));
            }
          }
        }
      }
    }
    _nbList.sortOnDistance(segId);
  }
}

void BackgroundPolarizer::computeStaticFieldsRS(
    std::vector<EwaldSegment>& ewaldSegments) {

#pragma omp parallel for
  for (Index segId = 1100; segId < 1101; ++segId) {
    EwaldRSInteractor rs_interactor(_options.ewaldsplitting, _unit_cell);
    EwaldSegment& currentSeg = ewaldSegments[segId];
    for (const Neighbour& neighbour : _nbList.getNeighboursOf(segId)) {
      EwaldSegment& nbSeg = ewaldSegments[neighbour.getId()];
      for (EwaldSite& site : currentSeg) {
        for (EwaldSite& nbSite : nbSeg) {
          rs_interactor.staticField(site, nbSite, neighbour.getShift());
        }
      }
    }
  }
}

void BackgroundPolarizer::computeStaticFieldsKS(
    std::vector<EwaldSegment>& ewaldSegments) {
  EwaldKSpace ks_interactor(_options.ewaldsplitting, _unit_cell,
                                  ewaldSegments);
  ;
}

void BackgroundPolarizer::Polarize(std::vector<EwaldSegment>& ewaldSegments) {

  computeNeighbourList(ewaldSegments);

  computeStaticFieldsRS(ewaldSegments);

  computeStaticFieldsKS(ewaldSegments);

  std::cout << " ID " << ewaldSegments[1100].getId() << std::endl;
  for (auto& site : ewaldSegments[1100]) {
    std::cout << site << std::endl;
  }

  // computeInducedFields(ewaldSegments);
}
}  // namespace xtp
}  // namespace votca
