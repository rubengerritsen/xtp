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
#ifndef VOTCA_XTP_KSPACE_H
#define VOTCA_XTP_KSPACE_H

#include <boost/math/constants/constants.hpp>
#include <complex>
#include <vector>

// Local VOTCA includes

#include "backgroundpolarizer.h"
#include "ewd_segment.h"
#include "unitcell.h"

namespace votca {
namespace xtp {

class KVector;

class KSpace {
 public:
  KSpace(const EwaldOptions& options, const UnitCell& unitcell,
         std::vector<EwdSegment>& ewaldSegments);
  ~KSpace() = default;

  void computeStaticField();

  void testFunction();

 private:
  std::complex<double> computeSk(const Eigen::Vector3d& kvector) const;
  double computeAk(const Eigen::Vector3d& kvector) const;
  void computeKVectors();
  double a1, a2;  // alpha (splitting param) and its powers
  UnitCell _unit_cell;
  std::vector<EwdSegment>& _ewaldSegments;
  std::vector<KVector> _kvector_list;
  double fourPiVolume;
  Eigen::Vector3d  mu_tot = Eigen::Vector3d::Zero();
  double cutoff, cutoff2;
  const std::complex<double> ii =
      std::complex<double>(0.0, 1.0);  // imaginary i
  std::array<Index, 3> max_K;
};

/**
 * \brief Class that contains everything related to a single k-vector.
 * Its:
 *  - k-vector
 *  - A(k) value
 *  - S(k) value (the structure factor)
 * Implements comparison operators for std::sort.
 */
class KVector {
 public:
  KVector(const Eigen::Vector3d& kvector) : _kvector(kvector){};
  ~KVector() = default;

  const Eigen::Vector3d& getVector() const { return _kvector; }
  double getAk() const { return _Ak; }
  std::complex<double> getSk() const { return _Sk; }

  void setAk(double Ak) { _Ak = Ak; }
  void setSk(std::complex<double> Sk) { _Sk = Sk; }

  friend std::ostream& operator<<(std::ostream& out, const KVector& kvector) {
    out << std::scientific << std::setprecision(5)
        << "vec: " << kvector.getVector().transpose()
        << " Ak: " << kvector.getAk() << " Sk: " << kvector.getSk();
    return out;
  }

  bool operator<(const KVector& other) {
    if (this->_kvector.norm() < other.getVector().norm()) {
      return true;
    }
    return false;
  }

  bool operator>(const KVector& other) {
    if (this->_kvector.norm() > other.getVector().norm()) {
      return true;
    }
    return false;
  }

  bool operator==(const KVector& other) {
    if (this->_kvector.norm() == other.getVector().norm()) {
      return true;
    }
    return false;
  }

 private:
  Eigen::Vector3d _kvector;
  double _Ak = 0;
  std::complex<double> _Sk;
};

}  // namespace xtp
}  // namespace votca

#endif