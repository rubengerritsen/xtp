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

#include "rspace.h"
#include <vector>

namespace votca {
namespace xtp {

void RSpace::computeStaticField() {
#pragma omp parallel for
  for (Index segId = 0; segId < _ewaldSegments.size(); ++segId) {
    EwdSegment& currentSeg = _ewaldSegments[segId];
    for (const Neighbour& neighbour : _nbList.getNeighboursOf(segId)) {
      EwdSegment& nbSeg = _ewaldSegments[neighbour.getId()];
      for (EwdSite& site : currentSeg) {
        for (EwdSite& nbSite : nbSeg) {
          site.addToStaticField(
              staticFieldAtBy(site, nbSite, neighbour.getShift()));
        }
      }
    }
  }
}

Eigen::MatrixXd RSpace::getInducedDipoleInteraction() {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero();
#pragma omp parallel for
  for (Index segId = 0; segId < _ewaldSegments.size(); ++segId) {
    // The first part can be done in the same way as the static field ...
    EwdSegment& currentSeg = _ewaldSegments[segId];
    for (const Neighbour& neighbour : _nbList.getNeighboursOf(segId)) {
      EwdSegment& nbSeg = _ewaldSegments[neighbour.getId()];
      Index startRow = segmentOffSet[segId];
      for (EwdSite& site : currentSeg) {
        Index startCol = segmentOffSet[neighbour.getId()];
        for (EwdSite& nbSite : nbSeg) {
          result.block<3, 3>(startRow, startCol) +=
              inducedDipoleInteractionAtBy(site, nbSite, neighbour.getShift());
          startCol += 3;
        }
        startRow += 3;
      }
    }
    // ... but we also need the dipole interaction within a segment
    Index startRow = segmentOffSet[segId];
    Index startCol = startRow;
    for (Index site_ind1 = 0; site_ind1 < currentSeg.size(); ++site_ind1) {
      for (Index site_ind2 = site_ind1 + 1; site_ind2 < currentSeg.size();
           ++site_ind2) {
        result.block<3, 3>(startRow + 3 * site_ind1,
                           startCol + 3 * site_ind2) +=
            inducedDipoleInteractionAtBy(currentSeg[site_ind1],
                                         currentSeg[site_ind2]);
        result.block<3, 3>(startRow + 3 * site_ind2,
                           startCol + 3 * site_ind1) +=
            inducedDipoleInteractionAtBy(currentSeg[site_ind2],
                                         currentSeg[site_ind1]);
      }
    }
  }
  return result;
}

void RSpace::computeInducedField() {
#pragma omp parallel for
  for (Index segId = 0; segId < _ewaldSegments.size(); ++segId) {
    EwdSegment& currentSeg = _ewaldSegments[segId];
    for (const Neighbour& neighbour : _nbList.getNeighboursOf(segId)) {
      EwdSegment& nbSeg = _ewaldSegments[neighbour.getId()];
      for (EwdSite& site : currentSeg) {
        for (EwdSite& nbSite : nbSeg) {
          site.addToInducedField(
              inducedFieldAtBy(site, nbSite, neighbour.getShift()));
        }
      }
    }
  }
}

void RSpace::computeIntraMolecularField() {
#pragma omp parallel for
  for (Index segId = 0; segId < _ewaldSegments.size(); ++segId) {
    EwdSegment& currentSeg = _ewaldSegments[segId];
    for (Index site_ind1 = 0; site_ind1 < currentSeg.size(); ++site_ind1) {
      for (Index site_ind2 = site_ind1 + 1; site_ind2 < currentSeg.size();
           ++site_ind2) {
        currentSeg[site_ind1].addToInducedField(
            inducedFieldAtBy(currentSeg[site_ind1], currentSeg[site_ind2]));
        currentSeg[site_ind2].addToInducedField(
            inducedFieldAtBy(currentSeg[site_ind2], currentSeg[site_ind1]));
      }
    }
  }
}

/**************************************************
 * PRIVATE FUNCTIONS                              *
 **************************************************/

void RSpace::computeDistanceVariables(Eigen::Vector3d distVec) {
  dr = distVec;
  R1 = dr.norm();
  R2 = R1 * R1;
  rR1 = 1.0 / R1;
  rR2 = rR1 * rR1;
}

void RSpace::computeScreenedInteraction() {
  // Note RSpace screening is with erfc
  rSqrtPiExp = rSqrtPi * std::exp(-a2 * R2);

  rR1s = std::erfc(a1 * R1) * rR1;
  rR3s = rR2 * (rR1s + 2.0 * a1 * rSqrtPiExp);
  rR5s = rR2 * (3.0 * rR3s + 4.0 * a3 * rSqrtPiExp);
  rR7s = rR2 * (5.0 * rR5s + 8.0 * a5 * rSqrtPiExp);
}

void RSpace::computeTholeVariables(const Eigen::Matrix3d& pol1,
                                   const Eigen::Matrix3d& pol2) {
  thole_u3 =
      (R1 * R2) / std::sqrt((1.0 / 3.0) * (pol1.array() * pol2.array()).sum());

  if (thole_u3 < 40) {
    double thole_exp = std::exp(-thole * thole_u3);
    double thole_u6 = thole_u3 * thole_u3;
    l3 = 1 - thole_exp;
    l5 = 1 - (1 + thole * thole_u3) * thole_exp;
    l7 = 1 - (1 + thole * thole_u3 + (3. / 5.) * thole2 * thole_u6);
    l9 = 1 - (1 + thole * thole_u3 + (18. / 35.) * thole2 * thole_u6 +
              (9. / 35.) * thole3 * thole_u6 * thole_u3) *
                 thole_exp;
  } else {
    l3 = l5 = l7 = l9 = 1.0;
  }
}

void RSpace::setupNeighbourList() {
  std::array<Index, 3> maxCopies =
      _unit_cell.getNrOfRealSpaceCopiesForCutOff(cutoff);
  _nbList.setSize(_ewaldSegments.size());

#pragma omp parallel for
  for (Index segId = 0; segId < static_cast<Index>(_ewaldSegments.size());
       ++segId) {
    EwdSegment& currentSeg = _ewaldSegments[segId];
    for (EwdSegment seg : _ewaldSegments) {
      Eigen::Vector3d minImage_dr = _unit_cell.minImage(currentSeg, seg);
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
            Eigen::Vector3d dr_l = minImage_dr + lvector;
            double dist = dr_l.norm();
            if (dist < cutoff) {
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

Eigen::Vector3d RSpace::staticFieldAtBy(EwdSite& site, const EwdSite& nbSite,
                                        const Eigen::Vector3d shift) {
  computeDistanceVariables(_unit_cell.minImage(site, nbSite) + shift);
  computeScreenedInteraction();

  Eigen::Vector3d field = Eigen::Vector3d::Zero();
  Index rank = nbSite.getRank();
  // charge
  field += -nbSite.getCharge() * dr * rR3s;
  if (rank > 0) {  // dipole
    field += nbSite.getStaticDipole() * rR3s;
    field += -rR5s * dr * dr.dot(nbSite.getStaticDipole());
    if (rank > 1) {  // quadrupole
      // Using that the trace of a quadrupole contributes nothing, we can skip
      // that part
      field += rR5s * 2 * nbSite.getQuadrupole() * dr;
      field += -rR7s * dr * dr.dot(nbSite.getQuadrupole() * dr);
    }
  }
  return field;
}

Eigen::Vector3d RSpace::inducedFieldAtBy(EwdSite& site, const EwdSite& nbSite,
                                         const Eigen::Vector3d shift) {
  computeDistanceVariables(_unit_cell.minImage(site, nbSite) + shift);
  computeScreenedInteraction();
  computeTholeVariables(site.getPolarizationMatrix(),
                        nbSite.getPolarizationMatrix());

  Eigen::Vector3d field = Eigen::Vector3d::Zero();
  field += nbSite.getInducedDipole() * l3 * rR3s;
  field += dr * nbSite.getInducedDipole().dot(dr) * l5 * rR5s;
  return field;
}


Eigen::Matrix3d RSpace::inducedDipoleInteractionAtBy(
    EwdSite& site, const EwdSite& nbSite, const Eigen::Vector3d shift) {
  computeDistanceVariables(_unit_cell.minImage(site, nbSite) + shift);
  computeScreenedInteraction();
  computeTholeVariables();

  Eigen::Matrix3d interaction = Eigen::Matrix3d::Zero();
  interaction.diagonal().array() += l3 * rR3s;
  interaction += dr * dr.transpose() * l5 * rR5s;
  return interaction;
}

}  // namespace xtp
}  // namespace votca