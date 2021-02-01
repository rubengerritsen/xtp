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
#ifndef VOTCA_XTP_EWALDSEGMENT_H
#define VOTCA_XTP_EWALDSEGMENT_H

#include <vector>

// Local VOTCA includes

#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/ewaldsite.h"

namespace votca {
namespace xtp {

class EwaldSegment {
 public:
  EwaldSegment(const PolarSegment& pol) {
    for (const PolarSite& psite : pol) {
      EwaldSite esite(psite);
      _sites.push_back(esite);
    }
    _id = pol.getId();
    _position = pol.getPos();
  }
  ~EwaldSegment() = default;

  const Eigen::Vector3d& getPos() const { return _position; }

  Index getId() const {return _id;}

  const EwaldSite& at(Index index) const { return _sites.at(index); }
  EwaldSite& at(Index index) { return _sites.at(index); }

  const EwaldSite& operator[](Index index) const { return _sites[index]; }
  EwaldSite& operator[](Index index) { return _sites[index]; }

  typename std::vector<EwaldSite>::iterator begin() { return _sites.begin(); }
  typename std::vector<EwaldSite>::iterator end() { return _sites.end(); }

  typename std::vector<EwaldSite>::const_iterator begin() const {
    return _sites.begin();
  }
  typename std::vector<EwaldSite>::const_iterator end() const {
    return _sites.end();
  }

 private:
  Index _id;
  std::vector<EwaldSite> _sites;
  Eigen::Vector3d _position;
};
}  // namespace xtp
}  // namespace votca
#endif