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
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/bse.h"
#include "votca/xtp/ecpbasisset.h"
#include "votca/xtp/gwbse.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"

using boost::format;
using namespace boost::filesystem;
using std::flush;
namespace votca {
namespace xtp {

Index GWBSE::CountCoreLevels() {
  Index ignored_corelevels = 0;
  if (!orbitals_.hasECPName()) {
    ECPBasisSet basis;
    basis.Load("corelevels");
    Index coreElectrons = 0;
    for (const auto& atom : orbitals_.QMAtoms()) {
      coreElectrons += basis.getElement(atom.getElement()).getNcore();
    }
    ignored_corelevels = coreElectrons / 2;
  }
  return ignored_corelevels;
}

void GWBSE::Initialize(tools::Property& options) {

  std::string key = Identify();

  // getting level ranges
  Index rpamax = 0;
  Index rpamin = 0;  // never changes
  Index qpmin = 0;
  Index qpmax = 0;
  Index bse_vmin = 0;
  Index bse_cmax = 0;

  Index homo = orbitals_.getHomo();  // indexed from 0
  Index num_of_levels = orbitals_.getBasisSetSize();
  Index num_of_occlevels = orbitals_.getNumberOfAlphaElectrons();

  std::string ranges = options.get(key + ".ranges").as<std::string>();

  // now check validity, and get rpa, qp, and bse level ranges accordingly

  if (ranges == "factor") {

    double rpamaxfactor = options.get(key + ".rpamax").as<double>();
    rpamax = Index(rpamaxfactor * double(num_of_levels)) -
             1;  // total number of levels

    double qpminfactor = options.get(key + ".qpmin").as<double>();
    qpmin =
        num_of_occlevels - Index(qpminfactor * double(num_of_occlevels)) - 1;

    double qpmaxfactor = options.get(key + ".qpmax").as<double>();
    qpmax =
        num_of_occlevels + Index(qpmaxfactor * double(num_of_occlevels)) - 1;

    double bseminfactor = options.get(key + ".bsemin").as<double>();
    bse_vmin =
        num_of_occlevels - Index(bseminfactor * double(num_of_occlevels)) - 1;

    double bsemaxfactor = options.get(key + ".bsemax").as<double>();
    bse_cmax =
        num_of_occlevels + Index(bsemaxfactor * double(num_of_occlevels)) - 1;

  } else if (ranges == "explicit") {
    // get explicit numbers
    rpamax = options.get(key + ".rpamax").as<Index>();
    qpmin = options.get(key + ".qpmin").as<Index>();
    qpmax = options.get(key + ".qpmax").as<Index>();
    bse_vmin = options.get(key + ".bsemin").as<Index>();
    bse_cmax = options.get(key + ".bsemax").as<Index>();
  } else if (ranges == "default") {
    rpamax = num_of_levels - 1;
    qpmin = 0;
    qpmax = 3 * homo + 1;
    bse_vmin = 0;
    bse_cmax = 3 * homo + 1;
  } else if (ranges == "full") {
    rpamax = num_of_levels - 1;
    qpmin = 0;
    qpmax = num_of_levels - 1;
    bse_vmin = 0;
    bse_cmax = num_of_levels - 1;
  }
  std::string ignore_corelevels =
      options.get(key + ".ignore_corelevels").as<std::string>();

  if (ignore_corelevels == "RPA" || ignore_corelevels == "GW" ||
      ignore_corelevels == "BSE") {
    Index ignored_corelevels = CountCoreLevels();
    if (ignore_corelevels == "RPA") {
      rpamin = ignored_corelevels;
    }
    if (ignore_corelevels == "GW" || ignore_corelevels == "RPA") {
      if (qpmin < ignored_corelevels) {
        qpmin = ignored_corelevels;
      }
    }
    if (ignore_corelevels == "GW" || ignore_corelevels == "RPA" ||
        ignore_corelevels == "BSE") {
      if (bse_vmin < ignored_corelevels) {
        bse_vmin = ignored_corelevels;
      }
    }

    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Ignoring " << ignored_corelevels
        << " core levels for " << ignore_corelevels << " and beyond." << flush;
  }

  // check maximum and minimum sizes
  if (rpamax > num_of_levels) {
    rpamax = num_of_levels - 1;
  }
  if (qpmax > num_of_levels) {
    qpmax = num_of_levels - 1;
  }
  if (bse_cmax > num_of_levels) {
    bse_cmax = num_of_levels - 1;
  }
  if (bse_vmin < 0) {
    bse_vmin = 0;
  }
  if (qpmin < 0) {
    qpmin = 0;
  }

  gwopt_.homo = homo;
  gwopt_.qpmin = qpmin;
  gwopt_.qpmax = qpmax;
  gwopt_.rpamin = rpamin;
  gwopt_.rpamax = rpamax;

  bseopt_.vmin = bse_vmin;
  bseopt_.cmax = bse_cmax;
  bseopt_.homo = homo;
  bseopt_.qpmin = qpmin;
  bseopt_.qpmax = qpmax;
  bseopt_.rpamin = rpamin;
  bseopt_.rpamax = rpamax;

  orbitals_.setRPAindices(rpamin, rpamax);
  orbitals_.setGWindices(qpmin, qpmax);
  orbitals_.setBSEindices(bse_vmin, bse_cmax);
  orbitals_.SetFlagUseHqpOffdiag(bseopt_.use_Hqp_offdiag);

  Index bse_vmax = homo;
  Index bse_cmin = homo + 1;
  Index bse_vtotal = bse_vmax - bse_vmin + 1;
  Index bse_ctotal = bse_cmax - bse_cmin + 1;
  Index bse_size = bse_vtotal * bse_ctotal;

  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " RPA level range [" << rpamin
                              << ":" << rpamax << "]" << flush;
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " GW  level range [" << qpmin
                              << ":" << qpmax << "]" << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " BSE level range occ[" << bse_vmin << ":" << bse_vmax
      << "]  virt[" << bse_cmin << ":" << bse_cmax << "]" << flush;

  gwopt_.reset_3c = options.get(key + ".rebuild_threecenter_freq").as<Index>();

  bseopt_.nmax = options.get(key + ".exctotal").as<Index>();
  if (bseopt_.nmax > bse_size || bseopt_.nmax < 0) {
    bseopt_.nmax = bse_size;
  }

  bseopt_.davidson_correction =
      options.get(key + ".eigensolver.davidson_correction").as<std::string>();

  bseopt_.davidson_tolerance =
      options.get(key + ".eigensolver.davidson_tolerance").as<std::string>();

  bseopt_.davidson_update =
      options.get(key + ".eigensolver.davidson_update").as<std::string>();

  bseopt_.davidson_maxiter =
      options.get(key + ".eigensolver.davidson_maxiter").as<Index>();

  bseopt_.useTDA = options.get(key + ".useTDA").as<bool>();
  orbitals_.setTDAApprox(bseopt_.useTDA);
  if (!bseopt_.useTDA) {
    XTP_LOG(Log::error, *pLog_) << " BSE type: full" << flush;
  } else {
    XTP_LOG(Log::error, *pLog_) << " BSE type: TDA" << flush;
  }

  Index full_bse_size = (bseopt_.useTDA) ? bse_size : 2 * bse_size;
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " BSE Hamiltonian has size "
                              << full_bse_size << "x" << full_bse_size << flush;

  bseopt_.use_Hqp_offdiag = options.get(key + ".use_Hqp_offdiag").as<bool>();

  if (!bseopt_.use_Hqp_offdiag) {
    XTP_LOG(Log::error, *pLog_)
        << " BSE without Hqp offdiagonal elements" << flush;
  } else {
    XTP_LOG(Log::error, *pLog_)
        << " BSE with Hqp offdiagonal elements" << flush;
  }

  bseopt_.max_dyn_iter =
      options.get(key + ".dynamical_screening_max_iterations").as<Index>();
  bseopt_.dyn_tolerance =
      options.get(key + ".dynamical_screening_tolerance").as<double>();
  if (bseopt_.max_dyn_iter > 0) {
    do_dynamical_screening_bse_ = true;
  }

  functional_ = options.get(key + ".vxc.functional").as<std::string>();
  grid_ = options.get(key + ".vxc.grid").as<std::string>();

  auxbasis_name_ = options.get(key + ".auxbasisset").as<std::string>();
  dftbasis_name_ = options.get(key + ".basisset").as<std::string>();
  if (dftbasis_name_ != orbitals_.getDFTbasisName()) {
    throw std::runtime_error(
        "Name of the Basisset from .orb file: " + orbitals_.getDFTbasisName() +
        " and from GWBSE optionfile " + dftbasis_name_ + " do not agree.");
  }

  std::string mode = options.get(key + ".mode").as<std::string>();
  if (mode == "G0W0") {
    gwopt_.gw_sc_max_iterations = 1;
  } else if (mode == "evGW") {
    gwopt_.g_sc_limit = 0.1 * gwopt_.gw_sc_limit;
    gwopt_.eta = 0.1;
  }

  // Check if the exact exchange of the functional matches the exact exchange in
  // orb file. If not display a warning.
  double ScaHFX_temp = Vxc_Potential<Vxc_Grid>::getExactExchange(functional_);
  if (ScaHFX_temp != orbitals_.getScaHFX()) {
    XTP_LOG(Log::error, *pLog_)
        << (boost::format(
                "WARNING: GWBSE exact exchange a=%s differs from qmpackage "
                "exact exchange a=%s, \n probably your functionals are "
                "inconsistent \n or dft orbitals were loaded from a package"
                " other than votca or orca. \n The GWBSE exact exchange will "
                "be used.") %
            ScaHFX_temp % orbitals_.getScaHFX())
               .str()
        << std::flush;
    orbitals_.setScaHFX(ScaHFX_temp);
  }

  XTP_LOG(Log::error, *pLog_) << " Running GW as: " << mode << flush;
  gwopt_.ScaHFX = orbitals_.getScaHFX();

  gwopt_.shift = options.get(key + ".scissor_shift").as<double>();
  gwopt_.g_sc_limit =
      options.get(key + ".g_sc_limit").as<double>();  // convergence criteria
                                                      // for qp iteration
                                                      // [Hartree]]
  gwopt_.g_sc_max_iterations =
      options.get(key + ".g_sc_max_iterations").as<Index>();  // convergence
                                                              // criteria for qp
                                                              // iteration
                                                              // [Hartree]]

  if (mode == "evGW") {
    gwopt_.gw_sc_max_iterations =
        options.get(key + ".gw_sc_max_iterations").as<Index>();
  }

  gwopt_.gw_sc_limit =
      options.get(key + ".gw_sc_limit").as<double>();  // convergence criteria
                                                       // for shift it
  XTP_LOG(Log::error, *pLog_)
      << " g_sc_limit [Hartree]: " << gwopt_.g_sc_limit << flush;
  if (gwopt_.gw_sc_max_iterations > 1) {
    XTP_LOG(Log::error, *pLog_)
        << " gw_sc_limit [Hartree]: " << gwopt_.gw_sc_limit << flush;
  }
  bseopt_.min_print_weight =
      options.get(key + ".bse_print_weight").as<double>();
  // print exciton WF composition weight larger than minimum

  // possible tasks
  std::string tasks_string = options.get(key + ".tasks").as<std::string>();
  boost::algorithm::to_lower(tasks_string);
  if (tasks_string.find("all") != std::string::npos) {
    do_gw_ = true;
    do_bse_singlets_ = true;
    do_bse_triplets_ = true;
  }
  if (tasks_string.find("gw") != std::string::npos) {
    do_gw_ = true;
  }
  if (tasks_string.find("singlets") != std::string::npos) {
    do_bse_singlets_ = true;
  }
  if (tasks_string.find("triplets") != std::string::npos) {
    do_bse_triplets_ = true;
  }

  XTP_LOG(Log::error, *pLog_) << " Tasks: " << flush;
  if (do_gw_) {
    XTP_LOG(Log::error, *pLog_) << " GW " << flush;
  }
  if (do_bse_singlets_) {
    XTP_LOG(Log::error, *pLog_) << " singlets " << flush;
  }
  if (do_bse_triplets_) {
    XTP_LOG(Log::error, *pLog_) << " triplets " << flush;
  }
  XTP_LOG(Log::error, *pLog_) << " Store: " << flush;
  if (do_gw_) {
    XTP_LOG(Log::error, *pLog_) << " GW " << flush;
  }

  if (options.exists(key + ".fragments")) {
    std::vector<tools::Property*> prop_region =
        options.Select(key + ".fragments.fragment");
    Index index = 0;
    for (tools::Property* prop : prop_region) {
      std::string indices =
          prop->ifExistsReturnElseThrowRuntimeError<std::string>("indices");
      fragments_.push_back(QMFragment<BSE_Population>(index, indices));
      index++;
    }
  }

  gwopt_.sigma_integration =
      options.get(key + ".sigma_integrator").as<std::string>();
  XTP_LOG(Log::error, *pLog_)
      << " Sigma integration: " << gwopt_.sigma_integration << flush;
  gwopt_.eta = options.get(key + ".eta").as<double>();
  XTP_LOG(Log::error, *pLog_) << " eta: " << gwopt_.eta << flush;
  if (gwopt_.sigma_integration == "exact") {
    XTP_LOG(Log::error, *pLog_)
        << " RPA Hamiltonian size: " << (homo + 1 - rpamin) * (rpamax - homo)
        << flush;
  }
  gwopt_.order = options.get(key + ".quadrature_order").as<Index>();
  XTP_LOG(Log::error, *pLog_)
      << " Quadrature integration order : " << gwopt_.order << flush;
  gwopt_.quadrature_scheme =
      options.get(key + ".quadrature_scheme").as<std::string>();
  XTP_LOG(Log::error, *pLog_)
      << " Quadrature integration scheme : " << gwopt_.quadrature_scheme
      << flush;
  gwopt_.alpha = options.get(key + ".alpha").as<double>();
  XTP_LOG(Log::error, *pLog_)
      << " Alpha smoothing parameter : " << gwopt_.alpha << flush;

  gwopt_.qp_solver = options.get(key + ".qp_solver").as<std::string>();
  gwopt_.qp_grid_steps = options.get(key + ".qp_grid_steps").as<Index>();
  gwopt_.qp_grid_spacing = options.get(key + ".qp_grid_spacing").as<double>();
  XTP_LOG(Log::error, *pLog_) << " QP solver: " << gwopt_.qp_solver << flush;
  if (gwopt_.qp_solver == "grid") {
    XTP_LOG(Log::error, *pLog_)
        << " QP grid steps: " << gwopt_.qp_grid_steps << flush;
    XTP_LOG(Log::error, *pLog_)
        << " QP grid spacing: " << gwopt_.qp_grid_spacing << flush;
  }
  gwopt_.gw_mixing_order =
      options.get(key + ".gw_mixing_order").as<Index>();  // max history in
                                                          // mixing (0: plain,
                                                          // 1: linear, >1
                                                          // Anderson)

  gwopt_.gw_mixing_alpha = options.get(key + ".gw_mixing_alpha").as<double>();

  if (mode == "evGW") {
    if (gwopt_.gw_mixing_order == 0) {
      XTP_LOG(Log::error, *pLog_) << " evGW with plain update " << std::flush;
    } else if (gwopt_.gw_mixing_order == 1) {
      XTP_LOG(Log::error, *pLog_) << " evGW with linear update using alpha "
                                  << gwopt_.gw_mixing_alpha << std::flush;
    } else {
      XTP_LOG(Log::error, *pLog_) << " evGW with Anderson update with history "
                                  << gwopt_.gw_mixing_order << " using alpha "
                                  << gwopt_.gw_mixing_alpha << std::flush;
    }
  }

  sigma_plot_states_ =
      options.get(key + ".sigma_plot_states").as<std::string>();
  sigma_plot_steps_ = options.get(key + ".sigma_plot_steps").as<Index>();
  sigma_plot_spacing_ = options.get(key + ".sigma_plot_spacing").as<double>();
  sigma_plot_filename_ =
      options.get(key + ".sigma_plot_filename").as<std::string>();
  if (!sigma_plot_states_.empty()) {
    XTP_LOG(Log::error, *pLog_)
        << " Sigma plot states: " << sigma_plot_states_ << flush;
    XTP_LOG(Log::error, *pLog_)
        << " Sigma plot steps: " << sigma_plot_steps_ << flush;
    XTP_LOG(Log::error, *pLog_)
        << " Sigma plot spacing: " << sigma_plot_spacing_ << flush;
    XTP_LOG(Log::error, *pLog_)
        << " Sigma plot filename: " << sigma_plot_filename_ << flush;
  }
}

void GWBSE::addoutput(tools::Property& summary) {

  const double hrt2ev = tools::conv::hrt2ev;
  tools::Property& gwbse_summary = summary.add("GWBSE", "");
  if (do_gw_) {
    gwbse_summary.setAttribute("units", "eV");
    gwbse_summary.setAttribute(
        "DFTEnergy",
        (format("%1$+1.6f ") % (orbitals_.getDFTTotalEnergy() * hrt2ev)).str());

    tools::Property& dft_summary = gwbse_summary.add("dft", "");
    dft_summary.setAttribute("HOMO", gwopt_.homo);
    dft_summary.setAttribute("LUMO", gwopt_.homo + 1);

    for (Index state = 0; state < gwopt_.qpmax + 1 - gwopt_.qpmin; state++) {
      tools::Property& level_summary = dft_summary.add("level", "");
      level_summary.setAttribute("number", state + gwopt_.qpmin);
      level_summary.add(
          "dft_energy",
          (format("%1$+1.6f ") %
           (orbitals_.MOs().eigenvalues()(state + gwopt_.qpmin) * hrt2ev))
              .str());
      level_summary.add(
          "gw_energy",
          (format("%1$+1.6f ") % (orbitals_.QPpertEnergies()(state) * hrt2ev))
              .str());

      level_summary.add("qp_energy",
                        (format("%1$+1.6f ") %
                         (orbitals_.QPdiag().eigenvalues()(state) * hrt2ev))
                            .str());
    }
  }
  if (do_bse_singlets_) {
    tools::Property& singlet_summary = gwbse_summary.add("singlets", "");
    for (Index state = 0; state < bseopt_.nmax; ++state) {
      tools::Property& level_summary = singlet_summary.add("level", "");
      level_summary.setAttribute("number", state + 1);
      level_summary.add(
          "omega", (format("%1$+1.6f ") %
                    (orbitals_.BSESinglets().eigenvalues()(state) * hrt2ev))
                       .str());
      if (orbitals_.hasTransitionDipoles()) {

        const Eigen::Vector3d& dipoles = (orbitals_.TransitionDipoles())[state];
        double f = 2 * dipoles.squaredNorm() *
                   orbitals_.BSESinglets().eigenvalues()(state) / 3.0;

        level_summary.add("f", (format("%1$+1.6f ") % f).str());
        tools::Property& dipol_summary = level_summary.add(
            "Trdipole", (format("%1$+1.4f %2$+1.4f %3$+1.4f") % dipoles.x() %
                         dipoles.y() % dipoles.z())
                            .str());
        dipol_summary.setAttribute("unit", "e*bohr");
        dipol_summary.setAttribute("gauge", "length");
      }
    }
  }
  if (do_bse_triplets_) {
    tools::Property& triplet_summary = gwbse_summary.add("triplets", "");
    for (Index state = 0; state < bseopt_.nmax; ++state) {
      tools::Property& level_summary = triplet_summary.add("level", "");
      level_summary.setAttribute("number", state + 1);
      level_summary.add(
          "omega", (format("%1$+1.6f ") %
                    (orbitals_.BSETriplets().eigenvalues()(state) * hrt2ev))
                       .str());
    }
  }
  return;
}

/*
 *    Many-body Green's fuctions theory implementation
 *
 *  data required from orbitals file
 *  - atomic coordinates
 *  - DFT molecular orbitals (energies and coeffcients)
 *  - DFT exchange-correlation potential matrix in atomic orbitals
 *  - number of electrons, number of levels
 */

Eigen::MatrixXd GWBSE::CalculateVXC(const AOBasis& dftbasis) {
  if (orbitals_.getXCFunctionalName().empty()) {
    orbitals_.setXCFunctionalName(functional_);
  } else {
    if (!(functional_ == orbitals_.getXCFunctionalName())) {
      throw std::runtime_error("Functionals from DFT " +
                               orbitals_.getXCFunctionalName() + " GWBSE " +
                               functional_ + " differ!");
    }
  }

  Vxc_Grid grid;
  grid.GridSetup(grid_, orbitals_.QMAtoms(), dftbasis);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Setup grid for integration with gridsize: " << grid_
      << " with " << grid.getGridSize() << " points, divided into "
      << grid.getBoxesSize() << " boxes" << flush;
  Vxc_Potential<Vxc_Grid> vxcpotential(grid);
  vxcpotential.setXCfunctional(functional_);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Integrating Vxc with functional " << functional_
      << flush;
  Eigen::MatrixXd DMAT = orbitals_.DensityMatrixGroundState();
  Mat_p_Energy e_vxc_ao = vxcpotential.IntegrateVXC(DMAT);
  XTP_LOG(Log::info, *pLog_) << TimeStamp() << " Calculated Vxc" << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Set hybrid exchange factor: " << orbitals_.getScaHFX()
      << flush;
  Index qptotal = gwopt_.qpmax - gwopt_.qpmin + 1;
  Index basissize = Index(orbitals_.MOs().eigenvectors().rows());
  Eigen::MatrixXd mos =
      orbitals_.MOs().eigenvectors().block(0, gwopt_.qpmin, basissize, qptotal);

  Eigen::MatrixXd vxc = mos.transpose() * e_vxc_ao.matrix() * mos;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Calculated exchange-correlation expectation values "
      << flush;

  return vxc;
}

bool GWBSE::Evaluate() {

  // set the parallelization
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Using "
                              << OPENMP::getMaxThreads() << " threads" << flush;

  if (XTP_HAS_MKL_OVERLOAD()) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Using MKL overload for Eigen " << flush;
  } else {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Using native Eigen implementation, no BLAS overload " << flush;
  }
  Index nogpus = OpenMP_CUDA::UsingGPUs();
  if (nogpus > 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Using CUDA support for tensor multiplication with "
        << nogpus << " GPUs." << flush;
  }

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Molecule Coordinates [A] " << flush;
  for (QMAtom& atom : orbitals_.QMAtoms()) {
    std::string output = (boost::format("%5d"
                                        "%5s"
                                        "   %1.4f %1.4f %1.4f") %
                          atom.getId() % atom.getElement() %
                          (atom.getPos().x() * tools::conv::bohr2ang) %
                          (atom.getPos().y() * tools::conv::bohr2ang) %
                          (atom.getPos().z() * tools::conv::bohr2ang))
                             .str();

    XTP_LOG(Log::error, *pLog_) << output << flush;
  }

  std::string dft_package = orbitals_.getQMpackage();
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " DFT data was created by " << dft_package << flush;

  BasisSet dftbs;
  dftbs.Load(dftbasis_name_);

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Loaded DFT Basis Set " << dftbasis_name_ << flush;

  // fill DFT AO basis by going through all atoms
  AOBasis dftbasis;
  dftbasis.Fill(dftbs, orbitals_.QMAtoms());
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Filled DFT Basis of size "
                              << dftbasis.AOBasisSize() << flush;

  // load auxiliary basis set (element-wise information) from xml file
  BasisSet auxbs;
  auxbs.Load(auxbasis_name_);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Loaded Auxbasis Set " << auxbasis_name_ << flush;

  // fill auxiliary AO basis by going through all atoms
  AOBasis auxbasis;
  auxbasis.Fill(auxbs, orbitals_.QMAtoms());
  orbitals_.setAuxbasisName(auxbasis_name_);
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Filled Auxbasis of size "
                              << auxbasis.AOBasisSize() << flush;

  if ((do_bse_singlets_ || do_bse_triplets_) && fragments_.size() > 0) {
    for (const auto& frag : fragments_) {
      XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Fragment " << frag.getId()
                                  << " size:" << frag.size() << flush;
    }
  }

  if (!do_gw_ && !orbitals_.hasQPdiag()) {
    throw std::runtime_error(
        "You want no GW calculation but the orb file has no QPcoefficients for "
        "BSE");
  }
  TCMatrix_gwbse Mmn;
  // rpamin here, because RPA needs till rpamin
  Index max_3c = std::max(bseopt_.cmax, gwopt_.qpmax);
  Mmn.Initialize(auxbasis.AOBasisSize(), gwopt_.rpamin, max_3c, gwopt_.rpamin,
                 gwopt_.rpamax);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp()
      << " Calculating Mmn_beta (3-center-repulsion x orbitals)  " << flush;
  Mmn.Fill(auxbasis, dftbasis, orbitals_.MOs().eigenvectors());
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Removed " << Mmn.Removedfunctions()
      << " functions from Aux Coulomb matrix to avoid near linear dependencies"
      << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Calculated Mmn_beta (3-center-repulsion x orbitals)  "
      << flush;

  Eigen::MatrixXd Hqp;
  if (do_gw_) {

    std::chrono::time_point<std::chrono::system_clock> start =
        std::chrono::system_clock::now();
    Eigen::MatrixXd vxc = CalculateVXC(dftbasis);
    GW gw = GW(*pLog_, Mmn, vxc, orbitals_.MOs().eigenvalues());
    gw.configure(gwopt_);
    gw.CalculateGWPerturbation();

    if (!sigma_plot_states_.empty()) {
      gw.PlotSigma(sigma_plot_filename_, sigma_plot_steps_, sigma_plot_spacing_,
                   sigma_plot_states_);
    }

    // store perturbative QP energy data in orbitals object (DFT, S_x,S_c, V_xc,
    // E_qp)
    orbitals_.QPpertEnergies() = gw.getGWAResults();
    orbitals_.RPAInputEnergies() = gw.RPAInputEnergies();

    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Calculating offdiagonal part of Sigma  " << flush;
    gw.CalculateHQP();
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Calculated offdiagonal part of Sigma  " << flush;

    Hqp = gw.getHQP();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es =
        gw.DiagonalizeQPHamiltonian();
    if (es.info() == Eigen::ComputationInfo::Success) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Diagonalized QP Hamiltonian  " << flush;
    }

    orbitals_.QPdiag().eigenvectors() = es.eigenvectors();
    orbitals_.QPdiag().eigenvalues() = es.eigenvalues();
    std::chrono::duration<double> elapsed_time =
        std::chrono::system_clock::now() - start;
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " GW calculation took "
                                << elapsed_time.count() << " seconds." << flush;

  } else {
    if (orbitals_.getGWAmax() != gwopt_.qpmax ||
        orbitals_.getGWAmin() != gwopt_.qpmin ||
        orbitals_.getRPAmax() != gwopt_.rpamax ||
        orbitals_.getRPAmin() != gwopt_.rpamin) {
      throw std::runtime_error(
          "The ranges for GW and RPA do not agree with the ranges from the "
          ".orb file, rerun your GW calculation");
    }
    const Eigen::MatrixXd& qpcoeff = orbitals_.QPdiag().eigenvectors();

    Hqp = qpcoeff * orbitals_.QPdiag().eigenvalues().asDiagonal() *
          qpcoeff.transpose();
  }

  // proceed only if BSE requested
  if (do_bse_singlets_ || do_bse_triplets_) {

    std::chrono::time_point<std::chrono::system_clock> start =
        std::chrono::system_clock::now();

    BSE bse = BSE(*pLog_, Mmn);
    bse.configure(bseopt_, orbitals_.RPAInputEnergies(), Hqp);

    // store the direct contribution to the static BSE results
    Eigen::VectorXd Hd_static_contrib_triplet;
    Eigen::VectorXd Hd_static_contrib_singlet;

    if (do_bse_triplets_) {
      bse.Solve_triplets(orbitals_);
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Solved BSE for triplets " << flush;
      bse.Analyze_triplets(fragments_, orbitals_);
    }

    if (do_bse_singlets_) {
      bse.Solve_singlets(orbitals_);
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Solved BSE for singlets " << flush;
      bse.Analyze_singlets(fragments_, orbitals_);
    }

    // do perturbative dynamical screening in BSE
    if (do_dynamical_screening_bse_) {

      if (do_bse_triplets_) {
        bse.Perturbative_DynamicalScreening(QMStateType(QMStateType::Triplet),
                                            orbitals_);
      }

      if (do_bse_singlets_) {
        bse.Perturbative_DynamicalScreening(QMStateType(QMStateType::Singlet),
                                            orbitals_);
      }
    }

    std::chrono::duration<double> elapsed_time =
        std::chrono::system_clock::now() - start;
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " BSE calculation took "
                                << elapsed_time.count() << " seconds." << flush;
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " GWBSE calculation finished " << flush;
  return true;
}

}  // namespace xtp
}  // namespace votca
