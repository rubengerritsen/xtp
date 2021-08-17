/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

// Standard includes
#include <chrono>

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/gnode.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "kmcmultiple.h"

namespace votca {
namespace xtp {
void KMCMultiple::ParseSpecificOptions(const tools::Property& options) {

  _runtime = options.get(".runtime").as<double>();
  _field = options.get(".field").as<Eigen::Vector3d>();
  double mtobohr = 1E9 * tools::conv::nm2bohr;
  _field *=
      (tools::conv::ev2hrt / mtobohr);  // Converting from V/m to Hartree/bohr

  _outputtime = options.get(".outputtime").as<double>();
  _timefile = options.ifExistsReturnElseReturnDefault<std::string>(".timefile",
                                                                   _timefile);

  std::string carriertype = options.get(".carriertype").as<std::string>();
  _carriertype = QMStateType(carriertype);
  if (!_carriertype.isKMCState()) {
    throw std::runtime_error("KMC cannot be run for state:" +
                             _carriertype.ToLongString());
  }

  _log.setReportLevel(Log::current_level);
  _log.setCommonPreface("\n ...");
}

void KMCMultiple::PrintDiffandMu(const Eigen::Matrix3d& avgdiffusiontensor,
                                 double simtime, unsigned long step) {
  double absolute_field = _field.norm();

  if (absolute_field == 0) {
    unsigned long diffusionsteps = step / _diffusionresolution;
    Eigen::Matrix3d result =
        avgdiffusiontensor /
        (double(diffusionsteps) * 2.0 * simtime * double(_numberofcarriers));
    XTP_LOG(Log::error, _log)
        << "\nStep: " << step
        << " Diffusion tensor averaged over all carriers (nm^2/s):\n"
        << result * tools::conv::bohr2nm * tools::conv::bohr2nm << std::flush;
  } else {
    double average_mobility = 0;
    double bohr2Hrts_to_nm2Vs =
        tools::conv::bohr2nm * tools::conv::bohr2nm / tools::conv::hrt2ev;
    XTP_LOG(Log::error, _log) << "\nMobilities (nm^2/Vs): " << std::flush;
    for (Index i = 0; i < _numberofcarriers; i++) {
      Eigen::Vector3d velocity = _carriers[i].get_dRtravelled() / simtime;
      double mobility =
          velocity.dot(_field) / (absolute_field * absolute_field);
      XTP_LOG(Log::error, _log)
          << std::scientific << "    carrier " << i + 1
          << ": mu=" << mobility * bohr2Hrts_to_nm2Vs << std::flush;
      average_mobility +=
          velocity.dot(_field) / (absolute_field * absolute_field);
    }
    average_mobility /= double(_numberofcarriers);
    XTP_LOG(Log::error, _log)
        << std::scientific
        << "  Overall average mobility in field direction <mu>="
        << average_mobility * bohr2Hrts_to_nm2Vs << " nm^2/Vs  " << std::flush;
  }
}

void KMCMultiple::WriteToTrajectory(std::fstream& traj,
                                    std::vector<Eigen::Vector3d>& startposition,
                                    double simtime, unsigned long step) const {
  traj << simtime << "\t";
  traj << step << "\t";
  for (Index i = 0; i < _numberofcarriers; i++) {
    Eigen::Vector3d pos = startposition[i] + _carriers[i].get_dRtravelled();
    traj << pos.x() * tools::conv::bohr2nm << "\t";
    traj << pos.y() * tools::conv::bohr2nm << "\t";
    traj << pos.z() * tools::conv::bohr2nm;
    if (i < _numberofcarriers - 1) {
      traj << "\t";
    } else {
      traj << std::endl;
    }
  }
}

void KMCMultiple::WriteToEnergyFile(std::fstream& tfile, double simtime,
                                    unsigned long step) const {
  double absolute_field = _field.norm();
  double currentenergy = 0;
  double currentmobility = 0;
  Eigen::Vector3d dr_travelled_current = Eigen::Vector3d::Zero();
  double dr_travelled_field = 0.0;
  Eigen::Vector3d avgvelocity_current = Eigen::Vector3d::Zero();
  if (absolute_field != 0) {
    for (const auto& carrier : _carriers) {
      dr_travelled_current += carrier.get_dRtravelled();
      currentenergy += carrier.getCurrentEnergy();
    }
    dr_travelled_current /= double(_numberofcarriers);
    currentenergy /= double(_numberofcarriers);
    avgvelocity_current = dr_travelled_current / simtime;
    currentmobility =
        avgvelocity_current.dot(_field) / (absolute_field * absolute_field);
    dr_travelled_field = dr_travelled_current.dot(_field) / absolute_field;
  }
  double bohr2Hrts_to_nm2Vs =
      tools::conv::bohr2nm * tools::conv::bohr2nm / tools::conv::hrt2ev;
  tfile << simtime << "\t" << step << "\t"
        << currentenergy * tools::conv::hrt2ev << "\t"
        << currentmobility * bohr2Hrts_to_nm2Vs << "\t"
        << dr_travelled_field * tools::conv::bohr2nm << "\t"
        << dr_travelled_current.norm() * tools::conv::bohr2nm << "\t"
        << std::endl;
}

void KMCMultiple::PrintDiagDandMu(const Eigen::Matrix3d& avgdiffusiontensor,
                                  double simtime, unsigned long step) {
  unsigned long diffusionsteps = step / _diffusionresolution;
  Eigen::Matrix3d result =
      avgdiffusiontensor /
      (double(diffusionsteps) * 2.0 * simtime * double(_numberofcarriers));

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(result);
  double bohr2_nm2 = tools::conv::bohr2nm * tools::conv::bohr2nm;
  XTP_LOG(Log::error, _log) << "\nEigenvalues:\n " << std::flush;
  for (Index i = 0; i < 3; i++) {
    XTP_LOG(Log::error, _log)
        << "Eigenvalue: " << es.eigenvalues()(i) * bohr2_nm2 << std::flush
        << "Eigenvector: " << es.eigenvectors().col(i).x() << "   "
        << es.eigenvectors().col(i).y() << "   " << es.eigenvectors().col(i).z()
        << "\n"
        << std::flush;
  }
  double bohr2Hrts_to_nm2Vs =
      tools::conv::bohr2nm * tools::conv::bohr2nm / tools::conv::hrt2ev;
  // calculate average mobility from the Einstein relation
  if (_field.norm() == 0) {
    XTP_LOG(Log::error, _log)
        << "The following value is calculated using the Einstein relation "
           "and assuming an isotropic medium"
        << std::flush;
    double avgD = 1. / 3. * es.eigenvalues().sum();
    double average_mobility = std::abs(avgD / _temperature);
    XTP_LOG(Log::error, _log)
        << std::scientific << "  Overall average mobility <mu>="
        << average_mobility * bohr2Hrts_to_nm2Vs << " nm^2/Vs " << std::flush;
  }
}

void KMCMultiple::PrintChargeVelocity(double simtime) {
  Eigen::Vector3d avg_dr_travelled = Eigen::Vector3d::Zero();
  for (Index i = 0; i < _numberofcarriers; i++) {
    XTP_LOG(Log::error, _log)
        << std::scientific << "    carrier " << i + 1 << ": "
        << _carriers[i].get_dRtravelled().transpose() / simtime *
               tools::conv::bohr2nm
        << std::flush;
    avg_dr_travelled += _carriers[i].get_dRtravelled();
  }
  avg_dr_travelled /= double(_numberofcarriers);

  Eigen::Vector3d avgvelocity = avg_dr_travelled / simtime;
  XTP_LOG(Log::error, _log)
      << std::scientific << "  Overall average velocity (nm/s): "
      << avgvelocity.transpose() * tools::conv::bohr2nm << std::flush;
}

void KMCMultiple::RunVSSM() {

  std::chrono::time_point<std::chrono::system_clock> realtime_start =
      std::chrono::system_clock::now();
  XTP_LOG(Log::error, _log)
      << "\nAlgorithm: VSSM for Multiple Charges" << std::flush;
  XTP_LOG(Log::error, _log)
      << "number of carriers: " << _numberofcarriers << std::flush;
  XTP_LOG(Log::error, _log)
      << "number of nodes: " << _nodes.size() << std::flush;

  bool checkifoutput = (_outputtime != 0);
  double nexttrajoutput = 0;
  unsigned long maxsteps = boost::numeric_cast<unsigned long>(_runtime);
  unsigned long outputstep = boost::numeric_cast<unsigned long>(_outputtime);
  bool stopontime = false;

  if (_runtime > 100) {
    XTP_LOG(Log::error, _log)
        << "stop condition: " << maxsteps << " steps." << std::flush;

    if (checkifoutput) {
      XTP_LOG(Log::error, _log) << "output frequency: ";
      XTP_LOG(Log::error, _log)
          << "every " << outputstep << " steps." << std::flush;
    }
  } else {
    stopontime = true;
    XTP_LOG(Log::error, _log)
        << "stop condition: " << _runtime << " seconds runtime." << std::flush;

    if (checkifoutput) {
      XTP_LOG(Log::error, _log) << "output frequency:\n "
                                   "every "
                                << _outputtime << " seconds." << std::flush;
    }
  }
  XTP_LOG(Log::error, _log)
      << "(If you specify runtimes larger than 100 kmcmultiple assumes that "
         "you are specifying the number of steps for both runtime and "
         "outputtime.)"
      << std::flush;

  if (!stopontime && _outputtime != 0 && floor(_outputtime) != _outputtime) {
    throw std::runtime_error(
        "ERROR in kmcmultiple: runtime was specified in steps (>100) and "
        "outputtime in seconds (not an integer). Please use the same units for "
        "both input parameters.");
  }

  if (_numberofcarriers > Index(_nodes.size())) {
    throw std::runtime_error(
        "ERROR in kmcmultiple: specified number of carriers is greater than "
        "the "
        "number of nodes. This conflicts with single occupation.");
  }

  std::fstream traj;
  std::fstream tfile;

  if (checkifoutput) {

    XTP_LOG(Log::error, _log)
        << "Writing trajectory to " << _trajectoryfile << "." << std::flush;
    traj.open(_trajectoryfile, std::fstream::out);

    traj << "time[s]\tsteps\t";
    for (Index i = 0; i < _numberofcarriers; i++) {
      traj << "carrier" << i + 1 << "_x\t";
      traj << "carrier" << i + 1 << "_y\t";
      traj << "carrier" << i + 1 << "_z";
      if (i < _numberofcarriers - 1) {
        traj << "'\t";
      }
    }
    traj << std::endl;
    if (!_timefile.empty()) {
      XTP_LOG(Log::error, _log)
          << "Writing time dependence of energy and mobility to " << _timefile
          << "." << std::flush;
      tfile.open(_timefile, std::fstream::out);
      tfile << "time[s]\t "
               "steps\tenergy_per_carrier[eV]\tmobility[nm**2/"
               "Vs]\tdistance_fielddirection[nm]\tdistance_absolute[nm]"
            << std::endl;
    }
  }
  RandomlyCreateCharges();
  std::vector<Eigen::Vector3d> startposition(_numberofcarriers,
                                             Eigen::Vector3d::Zero());
  for (Index i = 0; i < _numberofcarriers; i++) {
    startposition[i] = _carriers[i].getCurrentPosition();
  }

  traj << 0 << "\t";
  traj << 0 << "\t";
  for (Index i = 0; i < _numberofcarriers; i++) {
    traj << startposition[i].x() * tools::conv::bohr2nm << "\t";
    traj << startposition[i].y() * tools::conv::bohr2nm << "\t";
    traj << startposition[i].z() * tools::conv::bohr2nm;
    if (i < _numberofcarriers - 1) {
      traj << "\t";
    } else {
      traj << std::endl;
    }
  }

  std::vector<GNode*> forbiddennodes;
  std::vector<GNode*> forbiddendests;

  Eigen::Matrix3d avgdiffusiontensor = Eigen::Matrix3d::Zero();

  double simtime = 0.0;
  unsigned long step = 0;

  while (((stopontime && simtime < _runtime) ||
          (!stopontime && step < maxsteps))) {

    std::chrono::duration<double> elapsed_time =
        std::chrono::system_clock::now() - realtime_start;
    if (elapsed_time.count() > (_maxrealtime * 60. * 60.)) {
      XTP_LOG(Log::error, _log)
          << "\nReal time limit of " << _maxrealtime << " hours ("
          << Index(_maxrealtime * 60 * 60 + 0.5)
          << " seconds) has been reached. Stopping here.\n"
          << std::flush;
      break;
    }

    double cumulated_rate = 0;
    for (const auto& carrier : _carriers) {
      cumulated_rate += carrier.getCurrentEscapeRate();
    }
    if (cumulated_rate <= 0) {  // this should not happen: no possible jumps
                                // defined for a node
      throw std::runtime_error(
          "ERROR in kmcmultiple: Incorrect rates. All "
          "the escape rates for the current setting are 0.");
    }

    double dt = Promotetime(cumulated_rate);

    simtime += dt;
    step++;

    for (auto& carrier : _carriers) {
      carrier.updateOccupationtime(dt);
    }

    ResetForbiddenlist(forbiddennodes);
    bool level1step = true;
    while (level1step) {

      // determine which electron will escape
      GNode* newnode = nullptr;
      Chargecarrier* affectedcarrier = ChooseAffectedCarrier(cumulated_rate);

      if (CheckForbidden(affectedcarrier->getCurrentNode(), forbiddennodes)) {
        continue;
      }
      ResetForbiddenlist(forbiddendests);
      while (true) {
        // LEVEL 2

        const GLink& event =
            ChooseHoppingDest(affectedcarrier->getCurrentNode());
        newnode = event.getDestination();

        if (newnode == nullptr) {
          AddtoForbiddenlist(affectedcarrier->getCurrentNode(), forbiddennodes);
          break;  // select new escape node (ends level 2 but without setting
                  // level1step to 1)
        }

        // check after the event if this was allowed
        if (CheckForbidden(*newnode, forbiddendests)) {
          continue;
        }

        // if the new segment is unoccupied: jump; if not: add to forbidden
        // list and choose new hopping destination
        if (newnode->isOccupied()) {
          if (CheckSurrounded(affectedcarrier->getCurrentNode(),
                              forbiddendests)) {
            AddtoForbiddenlist(affectedcarrier->getCurrentNode(),
                               forbiddennodes);
            break;  // select new escape node (ends level 2 but without
                    // setting level1step to 1)
          }
          AddtoForbiddenlist(*newnode, forbiddendests);
          continue;  // select new destination
        } else {
          affectedcarrier->jumpAccordingEvent(event);
          level1step = false;
          break;  // this ends LEVEL 2 , so that the time is updated and the
                  // next MC step started
        }
        // END LEVEL 2
      }
      // END LEVEL 1
    }

    if (step % _diffusionresolution == 0) {
      for (const auto& carrier : _carriers) {
        avgdiffusiontensor += (carrier.get_dRtravelled()) *
                              (carrier.get_dRtravelled()).transpose();
      }
    }

    if (step != 0 && step % _intermediateoutput_frequency == 0) {
      PrintDiffandMu(avgdiffusiontensor, simtime, step);
    }

    if (checkifoutput) {
      bool outputsteps = (!stopontime && step % outputstep == 0);
      bool outputtime = (stopontime && simtime > nexttrajoutput);
      if (outputsteps || outputtime) {
        // write to trajectory file
        nexttrajoutput = simtime + _outputtime;
        WriteToTrajectory(traj, startposition, simtime, step);
        if (!_timefile.empty()) {
          WriteToEnergyFile(tfile, simtime, step);
        }
      }
    }
  }  // KMC

  if (checkifoutput) {
    traj.close();
    if (!_timefile.empty()) {
      tfile.close();
    }
  }

  WriteOccupationtoFile(simtime, _occfile);

  XTP_LOG(Log::error, _log) << "\nfinished KMC simulation after " << step
                            << " steps.\n"
                               "simulated time "
                            << simtime << " seconds.\n"
                            << std::flush;

  PrintChargeVelocity(simtime);

  XTP_LOG(Log::error, _log) << "\nDistances travelled (nm): " << std::flush;
  for (Index i = 0; i < _numberofcarriers; i++) {
    XTP_LOG(Log::error, _log)
        << std::scientific << "    carrier " << i + 1 << ": "
        << _carriers[i].get_dRtravelled().transpose() * tools::conv::bohr2nm
        << std::flush;
  }

  PrintDiffandMu(avgdiffusiontensor, simtime, step);
  PrintDiagDandMu(avgdiffusiontensor, simtime, step);

  return;
}

bool KMCMultiple::Evaluate(Topology& top) {

  XTP_LOG(Log::error, _log) << "\n-----------------------------------"
                               "\n      KMC FOR MULTIPLE CHARGES"
                               "\n-----------------------------------\n"
                            << std::flush;

  XTP_LOG(Log::info, _log) << "\nInitialising random number generator"
                           << std::flush;
  _RandomVariable.init(_seed);

  LoadGraph(top);
  RunVSSM();
  std::cout << _log;
  return true;
}

}  // namespace xtp
}  // namespace votca
