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
#ifndef VOTCA_XTP_BGNEIGHBOURLIST_H
#define VOTCA_XTP_BGNEIGHBOURLIST_H

#include <vector>
#include <algorithm>
// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"

namespace votca {
namespace xtp {

class Neighbour {
 public:
  Neighbour(Index id, Eigen::Vector3d dr, double dist)
      : segId(id), dr(dr), dist(dist){};

  ~Neighbour() = default;

  Index getId() const { return segId; }

  bool operator<(const Neighbour& other){
    if(this->dist < other.getDist()){
      return true;
    }
    return false;
  } 

  bool operator>(const Neighbour& other){
    if(this->dist > other.getDist()){
      return true;
    }
    return false;
  } 

  bool operator==(const Neighbour& other){
    if(this->segId == other.getId() && this->dr.isApprox(other.getDr(), 1e-5)){
      return true;
    }
    return false;
  }

  const Eigen::Vector3d& getDr() const { return dr; }

  const double getDist() const { return dist; }

 private:
  Index segId;
  Eigen::Vector3d dr;
  double dist;
};

class BgNeighbourList {
 public:
  BgNeighbourList() = default;
  ~BgNeighbourList() = default;

  void setSize(Index size) { _nbList.resize(size); }

  Index getSize() const { return _nbList.size(); }

  void addNeighbourTo(Index segId, Neighbour nb) {
    _nbList[segId].push_back(nb);
  }

  void sortOnDistance(Index segId){
    std::sort(_nbList[segId].begin(), _nbList[segId].end());
  }

  const std::vector<Neighbour>& getNeighboursOf(Index segId) const {
    return _nbList[segId];
  }

 private:
  std::vector<std::vector<Neighbour>> _nbList;
};
}  // namespace xtp
}  // namespace votca
#endif