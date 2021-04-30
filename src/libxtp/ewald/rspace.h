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
#ifndef VOTCA_XTP_RSPACE_H
#define VOTCA_XTP_RSPACE_H

#include <boost/math/constants/constants.hpp>
#include <vector>

// Local VOTCA includes
#include "backgroundpolarizer.h"
#include "ewd_nblist.h"
#include "ewd_segment.h"
#include "unitcell.h"

namespace votca {
namespace xtp {

class RSpace {
 public:
  RSpace(const EwaldOptions& options, const UnitCell& unitcell,
         std::vector<EwdSegment>& ewaldSegments)
      : _unit_cell(unitcell), _ewaldSegments(ewaldSegments) {
    cutoff = options.r_cutoff;
    a1 = options.alpha;
    a2 = a1 * a1;
    a3 = a1 * a2;
    a4 = a2 * a2;
    a5 = a4 * a1;
    thole = options.sharpness;
    thole2 = thole * thole;
    thole3 = thole * thole2;

    setupNeighbourList();

    systemSize = 0;
    for (const auto& seg : ewaldSegments) {
      segmentOffSet.push_back(systemSize);
      systemSize += 3 * seg.size();
    }
  };

  ~RSpace() = default;

  void computeStaticField();

  void computeInducedField();

  void computeIntraMolecularField();

  Eigen::MatrixXd getInducedDipoleInteraction();

 private:
  void computeDistanceVariables(Eigen::Vector3d distVec);

  void computeScreenedInteraction();

  void computeTholeVariables();

  void setupNeighbourList();

  Eigen::Vector3d staticFieldAtBy(
      EwdSite& site, const EwdSite& nbSite,
      const Eigen::Vector3d shift = Eigen::Vector3d::Zero());

  Eigen::Vector3d inducedFieldAtBy(
      EwdSite& site, const EwdSite& nbSite,
      const Eigen::Vector3d shift = Eigen::Vector3d::Zero());

  /****************************/
  /* VARIABLES                */
  /****************************/
  std::vector<Index> segmentOffSet;
  Index systemSize;
  double cutoff;
  double a1, a2, a3, a4, a5;  // alpha (splitting param) and its powers
  double l3, l5, l7, l9;
  double thole, thole2, thole3, thole_u3;
  UnitCell _unit_cell;
  Eigen::Vector3d dr = Eigen::Vector3d::Zero();
  std::vector<EwdSegment>& _ewaldSegments;
  EwdNbList _nbList;
  double rR1, rR2;  // reciprocal (i.e. 1.0/ ...) distance and powers
  double R1, R2;    // distance and powers
  double rSqrtPiExp;
  double pi = boost::math::constants::pi<double>();
  double rSqrtPi = 1.0 / std::sqrt(pi);

  // rRns = reciprocal R, of order n, screened with erfc
  double rR1s, rR3s, rR5s, rR7s;
};
}  // namespace xtp
}  // namespace votca

#endif