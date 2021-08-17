#include "orb2mol.h"

#include "votca/xtp/molden.h"
#include <boost/format.hpp>

namespace votca {
namespace xtp {

void Orb2Mol::ParseOptions(const tools::Property&) {

  _moldenfile = _job_name + ".molden.input";
  _orbfile = _job_name + ".orb";
  _xyzfile = _job_name + ".xyz";
}

bool Orb2Mol::Run() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, _log) << "Loading data from " << _orbfile << std::flush;
  orbitals.ReadFromCpt(_orbfile);

  XTP_LOG(Log::error, _log) << "Start parsing" << std::flush;

  Molden writer(_log);
  writer.WriteFile(_moldenfile, orbitals);

  XTP_LOG(Log::error, _log) << "Done parsing \n" << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca
