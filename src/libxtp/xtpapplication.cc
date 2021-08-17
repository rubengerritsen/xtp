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
#include <boost/format.hpp>

// VOTCA includes
#include <stdexcept>
#include <votca/tools/globals.h>
#include <votca/tools/propertyiomanipulator.h>

// Local VOTCA includes
#include "votca/tools/property.h"
#include "votca/tools/tokenizer.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/version.h"
#include "votca/xtp/xtpapplication.h"

namespace votca {
namespace xtp {
XtpApplication::XtpApplication() = default;
XtpApplication::~XtpApplication() = default;

/**
 * \brief Adds program options to the executable
 *
 * Every executable requires option file for calculators it is running
 * It is thus a part of the base XtpApplication class
 *
 */
void XtpApplication::Initialize(void) {

  namespace propt = boost::program_options;
  AddProgramOptions()(
      "options,o", propt::value<std::string>(),
      std::string("  " + CalculatorType() + " options").c_str());

  AddProgramOptions()("nthreads,t", propt::value<Index>()->default_value(1),
                      "  number of threads to create");

  AddProgramOptions()(
      "execute,e", propt::value<std::string>(),
      std::string("Name of " + CalculatorType() + " to run").c_str());
  AddProgramOptions()(
      "list,l",
      std::string("Lists all available " + CalculatorType() + "s").c_str());
  AddProgramOptions()(
      "description,d", propt::value<std::string>(),
      std::string("Short description of a " + CalculatorType() + "s").c_str());

  AddProgramOptions()(
      "cmdoptions,c", propt::value<std::vector<std::string>>()->multitoken(),
      "Modify options via command line by e.g. '-c xmltag.subtag=value'. "
      "Use whitespace to separate multiple options");

  AddProgramOptions()(
      "printoptions,p", propt::value<std::string>(),
      std::string("Prints xml options of a " + CalculatorType()).c_str());

#ifdef USE_CUDA
  AddProgramOptions()("gpus,g", propt::value<Index>()->default_value(-1),
                      "  Number of gpus to use");
#endif

  AddCommandLineOptions();
}

bool XtpApplication::EvaluateOptions() {

  if (OptionsMap().count("list")) {
    std::cout << "Available XTP" + CalculatorType() + "s: \n";
    for (const auto& name : CalculatorNames()) {
      PrintDescription(std::cout, name, "xtp/xml", Application::HelpShort);
    }
    StopExecution();
    return true;
  }
#ifdef USE_CUDA
  OpenMP_CUDA::SetNoGPUs(OptionsMap()["gpus"].as<Index>());
#endif

  if (OptionsMap().count("description")) {
    std::string calcname = OptionsMap()["description"].as<std::string>();
    if (CalcExists(calcname)) {
      PrintDescription(std::cout, calcname, "xtp/xml", Application::HelpLong);
    } else {
      std::cout << CalculatorType() << " " << calcname << " does not exist\n";
    }

    StopExecution();
    return true;
  }
  if (OptionsMap().count("printoptions")) {
    std::string calcname = OptionsMap()["printoptions"].as<std::string>();
    if (CalcExists(calcname)) {
      tools::Property opt;
      std::cout << "XML options for " << calcname << std::endl;
      opt.LoadFromXML(tools::GetVotcaShare() + "/xtp/xml/" + calcname + ".xml");
      std::cout << opt.get("options") << std::endl;
    } else {
      std::cout << CalculatorType() << " " << calcname << " does not exist\n";
    }

    StopExecution();
    return true;
  }

  CheckRequired("execute", "Nothing to do here: Abort.");
  std::string calcname = OptionsMap()["execute"].as<std::string>();

  if (CalcExists(calcname)) {
    CreateCalculator(calcname);
  } else {
    std::cout << CalculatorType() << " " << calcname << " does not exist\n";
    StopExecution();
  }

  EvaluateSpecificOptions();

  if (OptionsMap().count("options")) {
    std::string optionsFile = OptionsMap()["options"].as<std::string>();
    options_.LoadFromXML(optionsFile);
  } else {
    // Empty user options
    tools::Property& opts = options_.add("options", "");
    opts.add(calcname, "");
  }

  if (OptionsMap().count("cmdoptions")) {
    for (const std::string& opt :
         OptionsMap()["cmdoptions"].as<std::vector<std::string>>()) {
      std::vector<std::string> entries = tools::Tokenizer(opt, "=").ToVector();
      if (entries.size() != 2) {
        throw std::runtime_error(opt + " is not well formated!");
      } else {
        options_.getOradd("options." + calcname + "." + entries[0]).value() =
            entries[1];
      }
    }
  }

  return true;
}

void XtpApplication::Run() {

  std::string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }
  xtp::HelpTextHeader(name);

  execute();
}

void XtpApplication::ShowHelpText(std::ostream& out) {
  std::string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }
  xtp::HelpTextHeader(name);
  HelpText(out);
  out << "\n\n" << VisibleOptions() << std::endl;
}

}  // namespace xtp
}  // namespace votca
