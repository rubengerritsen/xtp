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
/// For an earlier history see ctp repo commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <libint2/initialize.h>

// Local VOTCA includes
#include "votca/xtp/parallelxjobcalc.h"

using boost::format;

namespace votca {
namespace xtp {

template <typename JobContainer>
bool ParallelXJobCalc<JobContainer>::Evaluate(const Topology &top) {
  libint2::initialize();
  // INITIALIZE PROGRESS OBSERVER
  std::string progFile = _jobfile;
  std::unique_ptr<JobOperator> master = std::unique_ptr<JobOperator>(
      new JobOperator(-1, top, *this, _openmp_threads));
  master->getLogger().setReportLevel(Log::current_level);
  master->getLogger().setMultithreading(true);
  master->getLogger().setPreface(Log::info, "\nMST INF");
  master->getLogger().setPreface(Log::error, "\nMST ERR");
  master->getLogger().setPreface(Log::warning, "\nMST WAR");
  master->getLogger().setPreface(Log::debug, "\nMST DBG");
  _progObs->InitFromProgFile(progFile, *(master.get()));

  // CREATE + EXECUTE THREADS (XJOB HANDLERS)
  std::vector<std::unique_ptr<JobOperator>> jobOps;

  for (Index id = 0; id < _nThreads; id++) {
    jobOps.push_back(std::unique_ptr<JobOperator>(
        new JobOperator(id, top, *this, _openmp_threads)));
  }

  for (Index id = 0; id < _nThreads; ++id) {
    CustomizeLogger(*jobOps[id]);
  }

  if (!_maverick) {
    std::cout << std::endl;  // REQUIRED FOR PROGRESS BAR IN OBSERVER
  }

  for (Index id = 0; id < _nThreads; id++) {
    jobOps[id]->Start();
  }

  for (Index id = 0; id < _nThreads; id++) {
    jobOps[id]->WaitDone();
  }

  if (!_maverick) {
    for (Index id = 0; id < _nThreads; id++) {
      std::cout << std::endl << (jobOps[id]->getLogger()) << std::flush;
    }
  }

  jobOps.clear();

  // SYNC REMAINING COMPLETE JOBS
  _progObs->SyncWithProgFile(*(master.get()));
  libint2::finalize();
  return true;
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::JobOperator::Run() {
  OPENMP::setMaxThreads(_openmp_threads);
  while (true) {
    Job *job = _master._progObs->RequestNextJob(*this);

    if (job == nullptr) {
      break;
    } else {
      Result res = this->_master.EvalJob(_top, *job, *this);
      this->_master._progObs->ReportJobDone(*job, res, *this);
    }
  }
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::ParseCommonOptions(
    const tools::Property &options) {
  std::cout << "\n... ... Initialized with " << _nThreads << " threads.\n";

  _maverick = (_nThreads == 1) ? true : false;

  std::cout << "\n... ... Using " << _openmp_threads << " openmp threads for "
            << _nThreads << "x" << _openmp_threads << "="
            << _nThreads * _openmp_threads << " total threads." << std::flush;
  _jobfile = options.get(".job_file").as<std::string>();
  _mapfile = options.get(".map_file").as<std::string>();
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::CustomizeLogger(QMThread &thread) {

  // CONFIGURE LOGGER
  Logger &log = thread.getLogger();
  log.setReportLevel(Log::current_level);
  log.setMultithreading(_maverick);

  log.setPreface(Log::info,
                 (format("\nT%1$02d INF ...") % thread.getId()).str());
  log.setPreface(Log::error,
                 (format("\nT%1$02d ERR ...") % thread.getId()).str());
  log.setPreface(Log::warning,
                 (format("\nT%1$02d WAR ...") % thread.getId()).str());
  log.setPreface(Log::debug,
                 (format("\nT%1$02d DBG ...") % thread.getId()).str());
}

template <typename JobContainer>
tools::Property ParallelXJobCalc<JobContainer>::UpdateGWBSEOptions(
    const tools::Property &options) {
  tools::Property gwbse_options = options.get(".gwbse_options");
  gwbse_options.get(".gwbse").add("basisset",
                                  options.get("basisset").as<std::string>());
  gwbse_options.get(".gwbse").add("auxbasisset",
                                  options.get("auxbasisset").as<std::string>());
  gwbse_options.get(".gwbse.vxc")
      .add("functional", options.get("functional").as<std::string>());

  return gwbse_options;
}

template <typename JobContainer>
tools::Property ParallelXJobCalc<JobContainer>::UpdateDFTOptions(
    const tools::Property &options) {
  tools::Property package_options = options.get(".dftpackage");
  package_options.get("package").add("basisset",
                                     options.get("basisset").as<std::string>());
  package_options.get("package").add(
      "auxbasisset", options.get("auxbasisset").as<std::string>());
  package_options.get("package").add(
      "functional", options.get("functional").as<std::string>());

  return package_options;
}
// REGISTER PARALLEL CALCULATORS
template class ParallelXJobCalc<std::vector<Job>>;

}  // namespace xtp
}  // namespace votca
