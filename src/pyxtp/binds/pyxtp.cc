/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#include "pyxtp.hpp"
#include <iostream>
#include <pybind11/pybind11.h>
#include <vector>

using namespace votca;
namespace py = pybind11;

namespace pyxtp {

/**
 * @brief Construct a new pybind11 module object to invoke votca-xtp*
 *
 */
int call_calculator(const std::string& name, int nThreads) {
  pyxtp::PyXTP pyxtp;
  pyxtp.Initialize(name, nThreads);
  return 42;
}

void printer(const std::vector<std::string>& keys) {
  for (const auto& key : keys) {
    std::cout << "key: " << key << "\n";
  }
  std::cout << "size: " << keys.size() << "\n";
}

void PyXTP::Initialize(const std::string& name, int nThreads) {
  std::cout << "Votca Factory:\n";
  xtp::Calculatorfactory::RegisterAll();
  auto pt = xtp::Calculators().Create(name);
  std::vector<std::string> keys = xtp::Calculators().getKeys();
  printer(keys);
  std::cout << "Mock Factory:\n";
  pyxtp::Factory<Parent>::RegisterAll();
  keys = pyxtp::Calculators<Parent>().getKeys();
  printer(keys);
}
}  // namespace pyxtp

PYBIND11_MODULE(xtp_binds, module) {
  module.doc() =
      "VOTCA-XTP is a library which allows you to calculate the electronic "
      "properties of organic materials,"
      "https://votca.github.io";

  module.def("call_calculator", &pyxtp::call_calculator,
             R"pbdoc(
        Invoke a Votca XTP calculator

        Parameters
        ----------
        name
          Calculator's name

  )pbdoc");
}