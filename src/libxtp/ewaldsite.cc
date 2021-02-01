/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/backgroundpolarizer.h"
#include "votca/xtp/ewaldunitcell.h"
#include "votca/xtp/segmentmapper.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "votca/xtp/ewaldsite.h"

namespace votca {
namespace xtp {

EwaldSite::EwaldSite(const PolarSite& pol) {
  _id = pol.getId();
  _rank = pol.getRank();
  _position = pol.getPos();
  _charge = pol.getCharge();
  _dipole_static = pol.getStaticDipole();
  _dipole_induced = Eigen::Vector3d::Zero();
  // the 1/3 saves factors in further calculations
  // NB the quadrupole in the Ewald site should by mulplied by 3 to obtain the
  // true quadrupole
  _quadrupole = (1.0 / 3.0) * pol.CalculateCartesianMultipole();
}

}  // namespace xtp
}  // namespace votca