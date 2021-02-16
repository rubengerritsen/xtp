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
#ifndef VOTCA_XTP_EWALDSITE_H
#define VOTCA_XTP_EWALDSITE_H

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/logger.h"

namespace votca {
namespace xtp {

class EwdSite {
 public:
  EwdSite(const PolarSite& pol);
  ~EwdSite() = default;

  friend std::ostream& operator<<(std::ostream& out, const EwdSite& site) {
    out << site.getId() << std::fixed << std::setprecision(2)
        << " pos: " << 0.05291 * site.getPos().transpose() << " q: " << site.getCharge()
        << std::scientific << std::setprecision(8) << " D: " <<  0.05291 *site.getStaticDipole().transpose() 
        << " E: " << 5.142e11 * site.getStaticField().transpose();
    return out;
  }

  Index getId() const { return _id; }

  const Eigen::Vector3d& getPos() const { return _position; }

  const Eigen::Vector3d& getStaticDipole() const { return _dipole_static; }

  const Eigen::Vector3d getTotalDipole() const {return _dipole_static + _dipole_induced;}

  double getCharge() const { return _charge; }

  const Eigen::Matrix3d& getQuadrupole() const { return _quadrupole; }

  Index getRank() const { return _rank; }

  const Eigen::Vector3d& getStaticField() const { return _field_static; }

  void addToStaticField(Eigen::Vector3d field) { _field_static += field; }

 private:
  Index _id;
  Index _rank;
  Eigen::Vector3d _position;
  double _charge;
  Eigen::Vector3d _dipole_static;
  Eigen::Vector3d _dipole_induced = Eigen::Vector3d::Zero();
  Eigen::Matrix3d _quadrupole;
  Eigen::Vector3d _field_static = Eigen::Vector3d::Zero();
};
}  // namespace xtp
}  // namespace votca
#endif