/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#ifndef VOTCA_XTP_JOBAPPLICATION
#define VOTCA_XTP_JOBAPPLICATION

#include <votca/xtp/xtpapplication.h>

#include <votca/xtp/progressobserver.h>
#include <votca/xtp/topology.h>

#include "statesaversqlite.h"
#include <votca/xtp/jobcalculator.h>

namespace votca {
namespace xtp {

class JobApplication : public XtpApplication {
 public:
  JobApplication();
  void Initialize();
  bool EvaluateOptions();
  void Run(void);

  void BeginEvaluate(
      int nThreads,
      ProgObserver<std::vector<Job *>, Job *, Job::JobResult> *obs);
  bool EvaluateFrame();
  void AddCalculator(JobCalculator *calculator);

 protected:
  bool _generate_input = false;
  bool _run = false;
  bool _import = false;
  Topology _top;
  std::vector<std::unique_ptr<JobCalculator> > _calculators;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_JOBAPPLICATION
