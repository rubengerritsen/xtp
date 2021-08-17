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

// Private VOTCA includes
#include "ewd_site.h"

namespace votca {
namespace xtp {

class EwdSegment {
 public:
  EwdSegment(const PolarSegment& pol) {
    for (const PolarSite& psite : pol) {
      EwdSite esite(psite);
      _sites.push_back(esite);
    }
    _id = pol.getId();
    _position = pol.getPos();
  }
  ~EwdSegment() = default;

  const Eigen::Vector3d& getPos() const { return _position; }

  Index getId() const {return _id;}

  const EwdSite& at(Index index) const { return _sites.at(index); }
  EwdSite& at(Index index) { return _sites.at(index); }

  const EwdSite& operator[](Index index) const { return _sites[index]; }
  EwdSite& operator[](Index index) { return _sites[index]; }

  typename std::vector<EwdSite>::iterator begin() { return _sites.begin(); }
  typename std::vector<EwdSite>::iterator end() { return _sites.end(); }

  typename std::vector<EwdSite>::const_iterator begin() const {
    return _sites.begin();
  }
  typename std::vector<EwdSite>::const_iterator end() const {
    return _sites.end();
  }

  Index size() const { return _sites.size();}

  void calcPos() {
    tools::Elements element;
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    double totalmass = 0.0;
    for (const auto& site : _sites) {
      double mass = element.getMass(site.getElement());
      totalmass += mass;
      pos += mass * site.getPos();
    }
    _position = pos / totalmass;
  }

 private:
  Index _id;
  std::vector<EwdSite> _sites;
  Eigen::Vector3d _position;
};
}  // namespace xtp
}  // namespace votca
#endif