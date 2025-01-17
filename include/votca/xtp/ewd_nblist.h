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
#ifndef VOTCA_XTP_EWDNBLIST_H
#define VOTCA_XTP_EWDNBLIST_H

#include <algorithm>
#include <vector>
// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"

namespace votca {
namespace xtp {

class Neighbour {
 public:
  Neighbour(Index id, Eigen::Vector3d dr, Eigen::Vector3d shift, double dist)
      : segId(id), dr_(dr), shift_(shift), dist_(dist){};

  ~Neighbour() = default;

  Index getId() const { return segId; }

  bool operator<(const Neighbour& other) {
    if (this->dist_ < other.getDist()) {
      return true;
    }
    return false;
  }

  bool operator>(const Neighbour& other) {
    if (this->dist_ > other.getDist()) {
      return true;
    }
    return false;
  }

  bool operator==(const Neighbour& other) {
    if (this->segId == other.getId() &&
        this->dr_.isApprox(other.getDr(), 1e-5)) {
      return true;
    }
    return false;
  }

  const Eigen::Vector3d& getDr() const { return dr_; }

  const Eigen::Vector3d& getShift() const { return shift_; }

  double getDist() const { return dist_; }

 private:
  Index segId;
  Eigen::Vector3d dr_;
  Eigen::Vector3d shift_;
  double dist_;
};

class EwdNbList {
 public:
  EwdNbList() = default;
  ~EwdNbList() = default;

  void setSize(Index size) {
    size_ = size;
    nbList_.resize(size);
  }

  Index getSize() const { return nbList_.size(); }

  void addNeighbourTo(Index segId, Neighbour nb) {
    if (nb.getId() > size_) {
      throw std::runtime_error("Trying to add neighbour with index " +
                               std::to_string(nb.getId()) +
                               " but max size is " + std::to_string(size_));
    }
    nbList_[segId].push_back(nb);
  }

  void sortOnDistance(Index segId) {
    std::sort(nbList_[segId].begin(), nbList_[segId].end());
  }

  const std::vector<Neighbour>& getNeighboursOf(Index segId) const {
    return nbList_[segId];
  }

 private:
  std::vector<std::vector<Neighbour>> nbList_;
  Index size_;
};
}  // namespace xtp
}  // namespace votca
#endif