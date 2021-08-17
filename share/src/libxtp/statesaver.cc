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

// Third party includes
#include <boost/interprocess/sync/file_lock.hpp>

// VOTCA includes
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/statesaver.h"
#include "votca/xtp/topology.h"

namespace votca {
namespace xtp {

std::vector<Index> StateSaver::getFrames() const {
  CheckpointFile cpf(_hdf5file, CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader();
  std::vector<Index> frames = std::vector<Index>{};
  try {
    r(frames, "frames");
  } catch (std::runtime_error&) {
    ;
  }
  return frames;
}
void StateSaver::WriteFrame(const Topology& top) {
  if (!tools::filesystem::FileExists(_hdf5file)) {
    std::cout << "Creating statefile " << _hdf5file << std::endl;
    CheckpointFile cpf(_hdf5file, CheckpointAccessLevel::CREATE);
  }
  boost::interprocess::file_lock flock(_hdf5file.c_str());
  flock.lock();
  if (!TopStepisinFrames(top.getStep())) {
    std::vector<Index> frames = this->getFrames();
    frames.push_back(top.getStep());
    CheckpointFile cpf(_hdf5file, CheckpointAccessLevel::MODIFY);
    CheckpointWriter w = cpf.getWriter();
    w(frames, "frames");
    std::cout << "Frame with id " << top.getStep() << " was not in statefile "
              << _hdf5file << " ,adding it now." << std::endl;
  }

  CheckpointFile cpf(_hdf5file, CheckpointAccessLevel::MODIFY);
  CheckpointWriter w = cpf.getWriter("/frame_" + std::to_string(top.getStep()));
  top.WriteToCpt(w);

  flock.unlock();

  std::cout << "Wrote MD topology (step = " << top.getStep()
            << ", time = " << top.getTime() << ") to " << _hdf5file
            << std::endl;
  std::cout << "... ";

  std::cout << ". " << std::endl;
  return;
}

Topology StateSaver::ReadFrame(Index frameid) const {
  if (!tools::filesystem::FileExists(_hdf5file)) {
    throw std::runtime_error("Statefile " + _hdf5file + " does not exist.");
  }
  std::cout << "Import MD Topology (i.e. frame " << frameid << ")"
            << " from " << _hdf5file << std::endl;
  std::cout << "...";
  boost::interprocess::file_lock flock(_hdf5file.c_str());
  flock.lock();
  if (!TopStepisinFrames(frameid)) {
    throw std::runtime_error("Frame with id " + std::to_string(frameid) +
                             " is not in statefile.");
  }

  CheckpointFile cpf(_hdf5file, CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader("/frame_" + std::to_string(frameid));
  Topology top;
  top.ReadFromCpt(r);
  flock.unlock();
  std::cout << ". " << std::endl;
  return top;
}

bool StateSaver::TopStepisinFrames(Index frameid) const {
  std::vector<Index> frames = this->getFrames();
  return std::find(frames.begin(), frames.end(), frameid) != frames.end();
}

}  // namespace xtp
}  // namespace votca
