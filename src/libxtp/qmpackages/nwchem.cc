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

#include "nwchem.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <votca/tools/filesystem.h>
#include <votca/xtp/ecpaobasis.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
using namespace std;

void NWChem::Initialize(tools::Property& options) {

  // NWChem file names
  string fileName = "system";

  _input_file_name = fileName + ".nw";
  _log_file_name = fileName + ".log";
  _shell_file_name = fileName + ".sh";
  _mo_file_name = fileName + ".movecs";

  ParseCommonOptions(options);

  if (_write_guess) {
    std::string::size_type iop_pos = _options.find("iterations 1");
    if (iop_pos != std::string::npos) {
      _options = _options + "\n iterations 1 ";
    }
  }
}

/* For QM/MM the molecules in the MM environment are represented by
 * their atomic partial charge distributions. Triggered by the option
 * keyword "set bq background" NWChem expects them in x,y,z,q format in the
 * backround.crg file.
 */

void NWChem::WriteChargeOption() {
  std::string::size_type iop_pos = _options.find("set bq background");
  if (iop_pos != std::string::npos) {
    _options = _options + "\n set bq background";
  }
}

long NWChem::WriteBackgroundCharges(ofstream& nw_file) {

  Index numberofcharges = 0;
  boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
  for (const std::unique_ptr<StaticSite>& site : _externalsites) {
    Eigen::Vector3d pos = site->getPos() * tools::conv::bohr2ang;
    string sitestring =
        boost::str(fmt % pos.x() % pos.y() % pos.z() % site->getCharge());
    if (site->getCharge() != 0.0) {
      nw_file << sitestring << endl;
      numberofcharges++;
    }
    std::vector<MinimalMMCharge> split_multipoles = SplitMultipoles(*site);
    for (const auto& mpoles : split_multipoles) {
      Eigen::Vector3d pos2 = mpoles._pos * tools::conv::bohr2ang;
      string multipole =
          boost::str(fmt % pos2.x() % pos2.y() % pos2.z() % mpoles._q);
      nw_file << multipole << endl;
      numberofcharges++;
    }
  }
  nw_file << endl;
  return numberofcharges;
}

bool NWChem::WriteGuess(const Orbitals& orbitals) {
  ofstream orb_file;
  // get name of temporary ascii file and open it
  std::string filebase = tools::filesystem::GetFileBase(_mo_file_name);
  std::string orb_file_name_ascii = _run_dir + "/" + filebase + ".mos";
  orb_file.open(orb_file_name_ascii);

  // header
  orb_file << "#generated by VOTCA\nbasisum\ngeomsum\n\nscf\nFri Sep 13 "
              "00:00:00 2013\nscf\n1\n\n8\nao basis\n1\n"
           << flush;
  Index size_of_basis = orbitals.getBasisSetSize();
  orb_file << size_of_basis << endl;
  orb_file << size_of_basis << endl;
  Eigen::MatrixXd MOs = ReorderMOsBack(orbitals);
  Index level = 1;
  Index ncolumns = 3;
  // write occupations as double in three columns
  // occupied levels
  Index column = 1;
  for (Index i = 0; i < orbitals.getNumberOfAlphaElectrons(); i++) {
    orb_file << FortranFormat(2.0);
    if (column == ncolumns) {
      orb_file << endl;
      column = 0;
    }
    column++;
  }
  // unoccupied levels
  for (Index i = orbitals.getNumberOfAlphaElectrons(); i < size_of_basis; i++) {
    orb_file << FortranFormat(0.0);
    if (column == ncolumns) {
      orb_file << endl;
      column = 0;
    }
    column++;
  }
  // extra endl
  if (column != 1) {
    orb_file << endl;
  }

  // write all energies in same format
  column = 1;
  for (Index i = 0; i < orbitals.MOs().eigenvalues().size(); ++i) {
    double energy = orbitals.MOs().eigenvalues()[i];
    orb_file << FortranFormat(energy);
    if (column == ncolumns) {
      orb_file << endl;
      column = 0;
    }
    column++;
  }
  if (column != 1) {
    orb_file << endl;
  }

  // write coefficients in same format
  for (Index i = 0; i < MOs.cols(); ++i) {
    Eigen::VectorXd mr = MOs.col(i);
    column = 1;
    for (Index j = 0; j < mr.size(); ++j) {
      orb_file << FortranFormat(mr[j]);
      if (column == ncolumns) {
        orb_file << endl;
        column = 0;
      }
      column++;
    }
    level++;
    if (column != 1) {
      orb_file << endl;
    }
  }
  orb_file << " 0.0000   0.0000" << endl;
  orb_file.close();
  // now convert this ascii file to binary
  std::string command = "cd " + _run_dir +
                        "; asc2mov 5000 system.mos system.movecs > convert.log";
  Index i = std::system(command.c_str());
  if (i == 0) {
    XTP_LOG(logDEBUG, *_pLog)
        << "Converted MO file from ascii to binary" << flush;
  } else {
    XTP_LOG(logERROR, *_pLog)
        << "Conversion of binary MO file to binary failed. " << flush;
    return false;
  }
  return true;
}

/**
 * Prepares the *.nw file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied
 */
bool NWChem::WriteInputFile(const Orbitals& orbitals) {

  std::string temp_suffix = "/id";
  std::string scratch_dir_backup = _scratch_dir;
  std::ofstream nw_file;
  std::ofstream crg_file;

  std::string nw_file_name_full = _run_dir + "/" + _input_file_name;
  std::string crg_file_name_full = _run_dir + "/background.crg";

  nw_file.open(nw_file_name_full);
  // header
  nw_file << "geometry noautoz noautosym" << endl;

  const QMMolecule& qmatoms = orbitals.QMAtoms();

  for (const QMAtom& atom : qmatoms) {
    Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
    nw_file << setw(3) << atom.getElement() << setw(12)
            << setiosflags(ios::fixed) << setprecision(5) << pos.x() << setw(12)
            << setiosflags(ios::fixed) << setprecision(5) << pos.y() << setw(12)
            << setiosflags(ios::fixed) << setprecision(5) << pos.z() << endl;
  }
  nw_file << "end\n";
  if (_write_charges) {
    // part for the MM charge coordinates
    crg_file.open(crg_file_name_full);
    Index numberofcharges = WriteBackgroundCharges(crg_file);
    crg_file << endl;
    crg_file.close();
    nw_file << endl;
    nw_file << "set bq:max_nbq " << numberofcharges << endl;
    nw_file << "bq background" << endl;
    nw_file << "load background.crg format 1 2 3 4" << endl;
    nw_file << "end\n" << endl;
  }

  if (_write_basis_set) {
    WriteBasisset(nw_file, qmatoms);
  }

  if (_write_pseudopotentials) {
    WriteECP(nw_file, qmatoms);
  }

  // write charge of the molecule
  nw_file << "\ncharge " << _charge << "\n";

  // writing scratch_dir info
  if (_scratch_dir != "") {
    std::string _temp("scratch_dir " + _scratch_dir + temp_suffix + "\n");
    nw_file << _temp;
  }
  if (_charge != 0) {
    std::string dft = "dft";
    if (_options.find(dft) != std::string::npos) {
      Index dftpos = Index(_options.find(dft));
      dftpos += Index(dft.size());
      std::string openshell =
          "\nodft\n" + (boost::format("mult %1%\n") % _spin).str();
      _options.insert(dftpos, openshell, 0, openshell.size());
    } else {
      throw runtime_error("NWCHEM: dft input data missing");
    }
  }

  nw_file << _options << "\n";
  if (_write_guess) {
    bool worked = WriteGuess(orbitals);
    if (!worked) {
      return false;
    }
  }

  nw_file << endl;
  nw_file.close();
  // and now generate a shell script to run both jobs, if neccessary
  XTP_LOG(logDEBUG, *_pLog)
      << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

  _scratch_dir = scratch_dir_backup + temp_suffix;

  WriteShellScript();
  _scratch_dir = scratch_dir_backup;

  return true;
}

bool NWChem::WriteShellScript() {
  ofstream shell_file;
  std::string shell_file_name_full = _run_dir + "/" + _shell_file_name;
  shell_file.open(shell_file_name_full);
  shell_file << "#!/bin/bash" << endl;
  shell_file << "mkdir -p " << _scratch_dir << endl;
  Index threads = OPENMP::getMaxThreads();
  if (threads == 1) {
    shell_file << _executable << " " << _input_file_name << " > "
               << _log_file_name << " 2> run.error" << endl;
  } else {
    shell_file << "mpirun -np " << boost::lexical_cast<std::string>(threads)
               << " " << _executable << " " << _input_file_name << " > "
               << _log_file_name << " 2> run.error" << endl;
  }
  shell_file.close();
  return true;
}

/**
 * Runs the NWChem job.
 */
bool NWChem::Run() {

  XTP_LOG(logDEBUG, *_pLog) << "Running NWChem job" << flush;

  if (std::system(nullptr)) {

    // NWChem overrides input information, if *.db and *.movecs files are
    // present better trash the old version
    std::string file_name = _run_dir + "/system.db";
    remove(file_name.c_str());
    file_name = _run_dir + "/" + _log_file_name;
    remove(file_name.c_str());

    std::string command = "cd " + _run_dir + "; sh " + _shell_file_name;

    Index check = std::system(command.c_str());
    if (check == -1) {
      XTP_LOG(logERROR, *_pLog)
          << _input_file_name << " failed to start" << flush;
      return false;
    }

    if (CheckLogFile()) {
      XTP_LOG(logDEBUG, *_pLog) << "Finished NWChem job" << flush;
      /* maybe we DO need to convert from fortran binary to ASCII first to avoid
    compiler-dependent issues */
      std::string command2;
      ascii_mo_file_name =
          tools::filesystem::GetFileBase(_mo_file_name) + ".mos";
      command2 = "cd " + _run_dir + "; mov2asc 10000 " + _mo_file_name + " " +
                 ascii_mo_file_name + "> convert.log";
      Index i = std::system(command2.c_str());
      if (i == 0) {
        XTP_LOG(logDEBUG, *_pLog)
            << "Converted MO file from binary to ascii" << flush;
      } else {
        XTP_LOG(logERROR, *_pLog)
            << "Conversion of binary MO file to ascii failed. " << flush;
        return false;
      }

      return true;
    } else {
      XTP_LOG(logDEBUG, *_pLog) << "NWChem job failed" << flush;
      return false;
    }
  } else {
    XTP_LOG(logERROR, *_pLog)
        << _input_file_name << " failed to start" << flush;
    return false;
  }

  return true;
}

/**
 * Cleans up after the NWChem job
 */
void NWChem::CleanUp() {

  // cleaning up the generated files
  if (_cleanup.size() != 0) {
    tools::Tokenizer tok_cleanup(_cleanup, ",");
    std::vector<std::string> cleanup_info;
    tok_cleanup.ToVector(cleanup_info);
    for (const std::string& substring : cleanup_info) {
      if (substring == "nw") {
        std::string file_name = _run_dir + "/" + _input_file_name;
        remove(file_name.c_str());
      }

      if (substring == "db") {
        std::string file_name = _run_dir + "/system.db";
        remove(file_name.c_str());
      }

      if (substring == "log") {
        std::string file_name = _run_dir + "/" + _log_file_name;
        remove(file_name.c_str());
      }

      if (substring == "movecs") {
        std::string file_name = _run_dir + "/" + _mo_file_name;
        remove(file_name.c_str());
        std::string file_name2 = _run_dir + "/" + ascii_mo_file_name;
        remove(file_name2.c_str());
      }

      if (substring == "gridpts") {
        std::string file_name = _run_dir + "/system.gridpts.*";
        remove(file_name.c_str());
      }
    }
  }
}

/**
 * Reads in the MO coefficients from a NWChem movecs file
 */
bool NWChem::ParseMOsFile(Orbitals& orbitals) {
  std::map<Index, std::vector<double> > coefficients;
  std::map<Index, double> energies;
  std::map<Index, double> occupancy;

  std::string line;
  Index levels = 0;
  // Index _level;
  Index basis_size = 0;
  Index number_of_electrons = 0;

  // opening the ascii MO file
  std::string orb_file_name_full = _run_dir + "/" + ascii_mo_file_name;
  std::ifstream input_file(orb_file_name_full);

  if (input_file.fail()) {
    XTP_LOG(logERROR, *_pLog)
        << "File " << orb_file_name_full
        << " with molecular orbitals is not found " << flush;
    return false;
  } else {
    XTP_LOG(logDEBUG, *_pLog)
        << "Reading MOs from " << orb_file_name_full << flush;
  }

  // the first 12 lines are garbage info
  for (Index i = 0; i < 12; i++) {
    getline(input_file, line);
  }
  // next line has basis set size
  input_file >> basis_size;
  XTP_LOG(logDEBUG, *_pLog) << "Basis set size: " << basis_size << flush;

  // next line has number of stored MOs
  input_file >> levels;
  XTP_LOG(logDEBUG, *_pLog) << "Energy levels: " << levels << flush;

  /* next lines contain information about occupation of the MOs
   *  - each line has 3 numbers
   *  - from _occ we can determine the number of electrons/2 */
  Index n_lines = ((levels - 1) / 3);
  Index n_rest = levels - 3 * n_lines;
  // read in the data
  Index imo = 0;
  for (Index i = 0; i < n_lines; i++) {
    for (Index j = 0; j < 3; j++) {
      input_file >> occupancy[imo];
      if (occupancy[imo] == 2.0) {
        number_of_electrons++;
      }
      imo++;
    }
  }
  if (n_rest != 0) {
    for (Index i = 0; i < n_rest; i++) {
      input_file >> occupancy[imo];
      imo++;
    }
  }

  XTP_LOG(logDEBUG, *_pLog)
      << "Alpha electrons: " << number_of_electrons << flush;

  Index occupied_levels = number_of_electrons;
  Index unoccupied_levels = levels - occupied_levels;
  XTP_LOG(logDEBUG, *_pLog) << "Occupied levels: " << occupied_levels << flush;
  XTP_LOG(logDEBUG, *_pLog)
      << "Unoccupied levels: " << unoccupied_levels << flush;

  // reset index and read MO energies the same way
  imo = 0;
  for (Index i = 0; i < n_lines; i++) {
    for (Index j = 0; j < 3; j++) {
      input_file >> energies[imo];
      imo++;
    }
  }
  if (n_rest != 0) {
    for (Index i = 0; i < n_rest; i++) {
      input_file >> energies[imo];
      imo++;
    }
  }

  // Now, the same for the coefficients
  double coef;
  for (Index imo2 = 0; imo2 < levels; imo2++) {
    for (Index i = 0; i < n_lines; i++) {
      for (Index j = 0; j < 3; j++) {
        input_file >> coef;
        coefficients[imo2].push_back(coef);
      }
    }
    if (n_rest != 0) {
      for (Index i = 0; i < n_rest; i++) {
        input_file >> coef;
        coefficients[imo2].push_back(coef);
      }
    }
  }

  // copying information to the orbitals object
  orbitals.setBasisSetSize(basis_size);
  orbitals.setNumberOfAlphaElectrons(number_of_electrons);
  orbitals.setNumberOfOccupiedLevels(occupied_levels);
  // copying energies to a matrix
  orbitals.MOs().eigenvalues().resize(levels);
  //_level = 1;
  for (Index i = 0; i < levels; i++) {
    orbitals.MOs().eigenvalues()[i] = energies[i];
  }

  // copying orbitals to the matrix
  orbitals.MOs().eigenvectors().resize(basis_size, levels);
  for (Index i = 0; i < levels; i++) {
    for (Index j = 0; j < basis_size; j++) {
      orbitals.MOs().eigenvectors()(j, i) = coefficients[i][j];
    }
  }

  ReorderOutput(orbitals);
  XTP_LOG(logDEBUG, *_pLog) << "Done reading MOs" << flush;

  return true;
}

StaticSegment NWChem::GetCharges() const {
  StaticSegment result("charges", 0);
  if (!CheckLogFile()) {
    throw std::runtime_error("logfile not correctly formatted");
  }
  std::string line;
  ifstream input_file((_run_dir + "/" + _log_file_name));
  bool has_charges = false;
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);

    std::string::size_type charge_pos = line.find("ESP");
    if (charge_pos != std::string::npos) {
      XTP_LOG(logDEBUG, *_pLog) << "Getting charges" << flush;
      has_charges = true;
      // two empty lines
      getline(input_file, line);
      getline(input_file, line);

      // now starts the data in format
      // _id type x y z q

      std::vector<std::string> row = GetLineAndSplit(input_file, "\t ");
      Index nfields = Index(row.size());
      while (nfields == 6) {
        Index atom_id = boost::lexical_cast<Index>(row.at(0)) - 1;
        std::string atom_type = row.at(1);
        double x = boost::lexical_cast<double>(row.at(2));
        double y = boost::lexical_cast<double>(row.at(3));
        double z = boost::lexical_cast<double>(row.at(4));
        Eigen::Vector3d pos = {x, y, z};
        pos *= tools::conv::ang2bohr;
        double atom_charge = boost::lexical_cast<double>(row.at(5));
        row = GetLineAndSplit(input_file, "\t ");
        nfields = row.size();

        StaticSite temp = StaticSite(atom_id, atom_type, pos);
        temp.setCharge(atom_charge);
        result.push_back(temp);
      }
    }
  }
  if (!has_charges) {
    throw std::runtime_error("Could not find charges in logfile");
  }
  return result;
}

Eigen::Matrix3d NWChem::GetPolarizability() const {
  if (!CheckLogFile()) {
    throw std::runtime_error("logfile not correctly formatted");
  }
  std::string line;
  ifstream input_file((_run_dir + "/" + _log_file_name));
  bool has_pol = false;

  Eigen::Matrix3d pol = Eigen::Matrix3d::Zero();
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);

    std::string::size_type pol_pos =
        line.find("DFT Linear Response polarizability / au");
    if (pol_pos != std::string::npos) {
      XTP_LOG(logDEBUG, *_pLog) << "Getting polarizability" << flush;
      getline(input_file, line);
      tools::Tokenizer tok(line, " ");
      std::vector<std::string> line_split = tok.ToVector();
      double frequency = std::stod(line_split[2]);
      if (std::abs(frequency) > 1e-9) {
        XTP_LOG(logDEBUG, *_pLog)
            << "Warning: Polarizability was calculated for frequency "
            << frequency << " normally f=0 for static polarizability" << flush;
      }
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);
      for (Index i = 0; i < 3; i++) {
        getline(input_file, line);
        tools::Tokenizer tok2(line, " ");
        std::vector<std::string> values = tok2.ToVector();
        if (values.size() != 4) {
          throw std::runtime_error("Polarisation line " + line +
                                   " cannot be parsed");
        }
        Eigen::Vector3d row;
        row << std::stod(values[1]), std::stod(values[2]), std::stod(values[3]);
        pol.row(i) = row;
      }

      has_pol = true;
    }
  }
  if (!has_pol) {
    throw std::runtime_error("Could not find polarisation in logfile");
  }
  return pol;
}

bool NWChem::CheckLogFile() const {

  // check if the log file exists

  ifstream input_file((_run_dir + "/" + _log_file_name));

  if (input_file.fail()) {
    XTP_LOG(logERROR, *_pLog) << "NWChem LOG is not found" << flush;
    return false;
  };

  if (input_file.peek() == std::ifstream::traits_type::eof()) {
    XTP_LOG(logERROR, *_pLog)
        << "NWChem run failed. Check OpenMPI version!" << flush;
    return false;
  };

  /* Checking the log file is a pain in the *** since NWChem throws an error
   * for our 'iterations 1'  runs (even though it outputs the required data
   * correctly. The only way that works for both scf and noscf runs is to
   * check for "Total DFT energy" near the end of the log file.
   */

  input_file.seekg(0, ios_base::end);  // go to the EOF
  char ch;
  std::string::size_type total_energy_pos = std::string::npos;
  std::string::size_type diis_pos = std::string::npos;
  do {
    // get the beginning of the line
    do {
      input_file.seekg(-2, ios_base::cur);
      input_file.get(ch);
      // cout << "\nNext Char: " << ch << " TELL G " <<
      // (Index)_input_file.tellg()
      // << endl;
    } while (ch != '\n' || (Index)input_file.tellg() == -1);

    std::string line;
    getline(input_file, line);  // Read the current line
    total_energy_pos = line.find("Total DFT energy");
    diis_pos = line.find("diis");
    // whatever is found first, determines the completeness of the file
    if (total_energy_pos != std::string::npos) {
      return true;
    } else if (diis_pos != std::string::npos) {
      XTP_LOG(logERROR, *_pLog) << "NWChem LOG is incomplete" << flush;
      return false;
    } else {
      // go to previous line
      //_input_file.get(ch);
      do {
        input_file.seekg(-2, ios_base::cur);
        input_file.get(ch);
        // cout << "\nChar: " << ch << endl;
      } while (ch != '\n' || (Index)input_file.tellg() == -1);
    }
  } while (total_energy_pos == std::string::npos &&
           diis_pos == std::string::npos);

  input_file.close();
  return true;
}

/**
 * Parses the Gaussian Log file and stores data in the Orbitals object
 */
bool NWChem::ParseLogFile(Orbitals& orbitals) {

  std::string line;
  std::vector<std::string> results;

  Index basis_set_size = 0;

  XTP_LOG(logDEBUG, *_pLog) << "Parsing " << _log_file_name << flush;
  std::string log_file_name_full = _run_dir + "/" + _log_file_name;
  // check if LOG file is complete
  if (!CheckLogFile()) {
    return false;
  }

  // save qmpackage name
  orbitals.setQMpackage(getPackageName());
  orbitals.setDFTbasisName(_basisset_name);
  if (_write_pseudopotentials) {
    orbitals.setECPName(_ecp_name);
  }

  // Start parsing the file line by line
  ifstream input_file(log_file_name_full);
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);
    /*
     * basis set size (is only required for overlap matrix reading, rest is
     * in orbitals file and could be skipped
     */
    std::string::size_type basis_pos = line.find("number of functions");
    if (basis_pos != std::string::npos) {
      tools::Tokenizer tok(line, ":");
      results = tok.ToVector();
      std::string bf = results.back();
      boost::trim(bf);
      basis_set_size = boost::lexical_cast<Index>(bf);
      orbitals.setBasisSetSize(basis_set_size);
      XTP_LOG(logDEBUG, *_pLog)
          << "Basis functions: " << basis_set_size << flush;
    }

    std::string::size_type energy_pos = line.find("Total DFT energy");
    if (energy_pos != std::string::npos) {

      tools::Tokenizer tok(line, "=");
      results = tok.ToVector();
      std::string energy = results.back();
      boost::trim(energy);
      orbitals.setQMEnergy(boost::lexical_cast<double>(energy));
      XTP_LOG(logDEBUG, *_pLog) << (boost::format("QM energy[Hrt]: %4.6f ") %
                                    orbitals.getDFTTotalEnergy())
                                       .str()
                                << flush;
    }

    std::string::size_type coordinates_pos = line.find("Output coordinates");
    if (coordinates_pos != std::string::npos) {
      XTP_LOG(logDEBUG, *_pLog) << "Getting the coordinates" << flush;
      //_has_coordinates = true;
      bool has_QMAtoms = orbitals.hasQMAtoms();
      // three garbage lines
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);
      // now starts the data in format
      // _id type Qnuc x y z
      std::vector<std::string> row = GetLineAndSplit(input_file, "\t ");
      Index nfields = Index(row.size());

      while (nfields == 6) {
        Index atom_id = boost::lexical_cast<Index>(row.at(0)) - 1;
        std::string atom_type = row.at(1);
        double x = boost::lexical_cast<double>(row.at(3));
        double y = boost::lexical_cast<double>(row.at(4));
        double z = boost::lexical_cast<double>(row.at(5));
        Eigen::Vector3d pos(x, y, z);
        pos *= tools::conv::ang2bohr;
        if (has_QMAtoms == false) {
          orbitals.QMAtoms().push_back(QMAtom(atom_id, atom_type, pos));
        } else {
          QMAtom& pAtom = orbitals.QMAtoms().at(atom_id);
          pAtom.setPos(pos);
        }
        atom_id++;
        row = GetLineAndSplit(input_file, "\t ");
        nfields = Index(row.size());
      }
    }

    // Check for ScaHFX = factor of HF exchange included in functional
    std::string::size_type HFX_pos = line.find("Hartree-Fock (Exact) Exchange");
    if (HFX_pos != std::string::npos) {

      tools::Tokenizer tok(line, "\t ");
      results = tok.ToVector();
      double ScaHFX = boost::lexical_cast<double>(results.back());
      orbitals.setScaHFX(ScaHFX);
      XTP_LOG(logDEBUG, *_pLog)
          << "DFT with " << ScaHFX << " of HF exchange!" << flush;
    }

  }  // end of reading the file line-by-line

  XTP_LOG(logDEBUG, *_pLog) << "Done parsing" << flush;
  return true;
}

void NWChem::WriteBasisset(ofstream& nw_file, const QMMolecule& qmatoms) {

  std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();
  BasisSet bs;
  bs.Load(_basisset_name);
  XTP_LOG(logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;
  nw_file << "basis spherical" << endl;
  for (const std::string& element_name : UniqueElements) {
    const Element& element = bs.getElement(element_name);
    for (const Shell& shell : element) {
      // nwchem can only use S,P,SP,D,F,G shells so we split them up if not SP
      if (!shell.isCombined()) {
        // shell type, number primitives, scale factor
        nw_file << element_name << " "
                << boost::algorithm::to_lower_copy(shell.getType()) << endl;
        for (const GaussianPrimitive& gaussian : shell) {
          for (Index icontr = 0; icontr < Index(gaussian.Contractions().size());
               icontr++) {
            if (gaussian.Contractions()[icontr] != 0.0) {
              nw_file << FortranFormat(gaussian.decay()) << " "
                      << FortranFormat(gaussian.Contractions()[icontr]) << endl;
            }
          }
        }
      } else {
        for (const char& subtype : shell.getType()) {
          nw_file << element_name << " "
                  << boost::algorithm::to_lower_copy(std::string(1, subtype))
                  << endl;

          for (const GaussianPrimitive& gaussian : shell) {
            nw_file << FortranFormat(gaussian.decay()) << " "
                    << FortranFormat(gaussian.Contractions()[FindLmax(
                           std::string(1, subtype))])
                    << endl;
          }
        }
      }
    }
  }
  nw_file << "end\n";
  nw_file << endl;

  return;
}

void NWChem::WriteECP(ofstream& nw_file, const QMMolecule& qmatoms) {

  std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();

  ECPBasisSet ecp;
  ecp.Load(_ecp_name);

  XTP_LOG(logDEBUG, *_pLog) << "Loaded Pseudopotentials " << _ecp_name << flush;

  for (const std::string& element_name : UniqueElements) {
    try {
      ecp.getElement(element_name);
    } catch (std::runtime_error& error) {
      XTP_LOG(logDEBUG, *_pLog)
          << "No pseudopotential for " << element_name << " available" << flush;
      continue;
    }
    const ECPElement& element = ecp.getElement(element_name);
    // element name, [possibly indeces of centers], zero to indicate the end
    nw_file << element_name << " nelec " << element.getNcore() << endl;
    for (const ECPShell& shell : element) {
      string shelltype = shell.getType();
      if (shell.getL() == element.getLmax()) {
        shelltype = "ul";
      }
      nw_file << element_name << " " << shelltype << endl;
      for (const ECPGaussianPrimitive& gaussian : shell) {
        nw_file << "    " << gaussian._power << " "
                << FortranFormat(gaussian._decay) << " "
                << FortranFormat(gaussian._contraction) << endl;
      }
    }
  }
  nw_file << "end\n";
  nw_file << endl;
  return;
}

std::string NWChem::FortranFormat(const double& number) {
  std::stringstream ssnumber;
  if (number >= 0) {
    ssnumber << "    ";
  } else {
    ssnumber << "   ";
  }
  ssnumber << setiosflags(ios::fixed) << setprecision(15) << std::scientific
           << number;
  std::string snumber = ssnumber.str();
  return snumber;
}

}  // namespace xtp
}  // namespace votca
