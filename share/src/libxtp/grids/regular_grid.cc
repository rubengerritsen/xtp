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

// Local VOTCA includes
#include "votca/xtp/regular_grid.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/radial_euler_maclaurin_rule.h"
#include "votca/xtp/sphere_lebedev_rule.h"

namespace votca {
namespace xtp {

void Regular_Grid::GridSetup(const Eigen::Array3d& stepsizes,
                             const Eigen::Array3d& padding,
                             const QMMolecule& atoms, const AOBasis& basis) {

  std::pair<Eigen::Vector3d, Eigen::Vector3d> extension =
      atoms.CalcSpatialMinMax();
  Eigen::Array3d min = extension.first.array();
  Eigen::Array3d max = extension.second.array();
  Eigen::Array3d doublesteps = (max - min + 2 * padding) / stepsizes + 1.0;
  Eigen::Array<votca::Index, 3, 1> steps = (doublesteps.ceil()).cast<Index>();

  // needed to symmetrize grid around molecule
  Eigen::Array3d padding_sym =
      (doublesteps - steps.cast<double>()) * stepsizes * 0.5 + padding;
  GridSetup(steps, padding_sym, atoms, basis);
}

void Regular_Grid::GridSetup(const Eigen::Array<Index, 3, 1>& steps,
                             const Eigen::Array3d& padding,
                             const QMMolecule& atoms, const AOBasis& basis) {

  std::pair<Eigen::Vector3d, Eigen::Vector3d> extension =
      atoms.CalcSpatialMinMax();
  Eigen::Array3d min = extension.first.array();
  Eigen::Array3d max = extension.second.array();
  _startingpoint = min - padding;
  _stepsizes = (max - min + 2 * padding) / (steps - 1).cast<double>();
  _steps = steps;
  const Index gridboxsize = 500;

  GridBox gridbox;
  for (Index i = 0; i < steps.x(); i++) {
    double x = _startingpoint.x() + double(i) * _stepsizes.x();
    for (Index j = 0; j < steps.y(); j++) {
      double y = _startingpoint.y() + double(j) * _stepsizes.y();
      for (Index k = 0; k < steps.z(); k++) {
        double z = _startingpoint.z() + double(k) * _stepsizes.z();
        GridContainers::Cartesian_gridpoint point;
        point.grid_weight = 1.0;
        point.grid_pos = Eigen::Vector3d(x, y, z);
        gridbox.addGridPoint(point);
        if (gridbox.size() == gridboxsize) {
          _grid_boxes.push_back(gridbox);
          gridbox = GridBox();
        }
      }
    }
  }
  if (gridbox.size() > 0) {
    _grid_boxes.push_back(gridbox);
  }
#pragma omp parallel for
  for (Index i = 0; i < getBoxesSize(); i++) {
    _grid_boxes[i].FindSignificantShells(basis);
    _grid_boxes[i].PrepareForIntegration();
  }

  _totalgridsize = 0;
  for (auto& box : _grid_boxes) {
    _totalgridsize += box.size();
  }
}

}  // namespace xtp
}  // namespace votca
