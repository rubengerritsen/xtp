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
#ifndef VOTCA_XTP_BACKGROUND_H
#define VOTCA_XTP_BACKGROUND_H
#include <vector>

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/region.h"
#include "votca/xtp/segid.h"

// Private VOTCA includes
#include "ewaldoptions.h"
#include "ewd_segment.h"
#include "kspace.h"
#include "rspace.h"
#include "unitcell.h"

namespace votca {
namespace xtp {

class Background {
 public:
  Background(Logger& log, const Eigen::Matrix3d& uc_matrix,
             const EwaldOptions options,
             std::vector<PolarSegment>& polar_background);

  Background(Logger& log) : log_(log), kspace(log), rspace(log) {}

  ~Background() = default;

  void Polarize();

  void ApplyBackgroundFields(
      std::vector<std::unique_ptr<votca::xtp::Region>>& regions,
      const std::vector<std::vector<SegId>>& region_seg_ids);

  void writeToStateFile(std::string& state_file);

  void readFromStateFile(const std::string& state_file);

 private:
  Index computeSystemSize(std::vector<EwdSegment>& ewaldSegments) const;
  Eigen::VectorXd solveLinearSystem(Eigen::MatrixXd A, Eigen::VectorXd b,
                                    Eigen::VectorXd guess);

  void bgFieldAtSegment(PolarSegment& seg, std::vector<SegId> pCloud_indices);
  Logger& log_;
  UnitCell unit_cell_;
  EwaldOptions options_;
  std::vector<EwdSegment> ewald_background_;
  RSpace rspace;
  KSpace kspace;
};
}  // namespace xtp
}  // namespace votca

#endif