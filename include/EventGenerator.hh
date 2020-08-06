// EventGenerator - Produces simulated high-energy particle collisions
//                  using the Pythia8 and FastJet libraries.
// Copyright (C) 2017-2020 Patrick T. Komiske III
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef EventGenerator_HH
#define EventGenerator_HH

// C++ standard library
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Pythia
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"

// FastJet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "JetMatcher.hh"

namespace mc {

// these modes determine the status of the jet
const unsigned char M_FULLY_MATCHED = 1;
const unsigned char M_HADRON_LEVEL  = 2;
const unsigned char M_PARTON_LEVEL  = 4;
const unsigned char M_HADRON_JET    = 8;
const unsigned char M_HADRON_HARD   = 16;
const unsigned char M_PARTON_JET    = 32;
const unsigned char M_PARTON_HARD   = 64;
const unsigned char M_PH_SAME_HARDS = 128;
const unsigned char M_PARTIAL_MATCH = ~M_FULLY_MATCHED;

// commonly used types
using fastjet::PseudoJet;

// essentially map massless particles masses to zero
inline std::string m_round(double x) {
  std::ostringstream ss;
  if (x < 1e-5) x = 0;
  ss << std::setprecision(6) << x;
  return ss.str();
}

// deal with wrapping around issues
// ensures phi is within pi of ref_phi
const double PI = 3.14159265358979323846;
inline double phi_fix(double phi, double ref_phi) {
  double diff = phi - ref_phi;
  if (diff > PI) phi -= 2*PI;
  else if (diff < -PI) phi += 2*PI;
  return phi; 
}

// holds the indices of a matched set of jets/hardprocess particles
struct MatchedUnit {
  double phi;
  int hard_parton_ind, parton_jet_ind, 
      hard_hadron_ind, hadron_jet_ind;
  unsigned char status;
};

///////////////////////////////////////////////////////////////////////////////
// EventGenerator - Primary class for doing event generation
///////////////////////////////////////////////////////////////////////////////
class EventGenerator {

///////////////////////////////////////////////////////////////////////////////
// PUBLIC PROPERTIES
///////////////////////////////////////////////////////////////////////////////
public:

  // event properties
  unsigned ntot, nmaxjets;
  bool fully_matched_only, hardprocess, partonlevel, hadronlevel;

  // jet properties
  double jet_R;
  std::vector<MatchedUnit> matched_units;

  // event containers
  std::vector<PseudoJet> origs, partons, hadrons, parton_jets, hadron_jets;
  std::vector<std::vector<PseudoJet>> decays;

  // vectors to store info about particle id
  std::vector<int> parton_pids, hadron_pids;

  // other options
  bool verbose, exit;
  unsigned print_every;
  std::string user_comments;

///////////////////////////////////////////////////////////////////////////////
// PRIVATE PROPERTIES
///////////////////////////////////////////////////////////////////////////////
private:

  // command line variables
  int argc_;
  char** argv_;

  // counters
  std::chrono::steady_clock::time_point start_time_;
  unsigned force_hadron_level_failed_, write_counter_;

  // filepath to output file
  std::string outfilepath_;
  std::ofstream outfile_;

  // random seed
  int seed_;

  // event properties
  std::string process_, description_;
  double ecm_, pthatmin_, pthatmax_, biaspow_;
  bool no_ue_, no_isr_, no_fsr_, 
       no_hardprocess_, no_partonlevel_, no_hadronlevel_,
       output_within_R_;

  // jet properties
  double ptjetmin_, ptjetmax_, absrapmax_;
  std::string jet_alg_name_;
  fastjet::JetAlgorithm jet_alg_;
  unsigned char base_status_;
  std::vector<double> hadron_jet_features_, parton_jet_features_;

  // printing options
  bool print_pythia_;
  unsigned num_print_event_, num_print_proc_, outprecision_;

  // setup options
  std::vector<unsigned> hardproc_ids_;
  std::vector<std::string> pythia_commands_;

  // Pythia instance
  Pythia8::Pythia pythia_;

  // FastJet objects
  fastjet::JetDefinition jet_def_;
  fastjet::Selector jet_sel_;
  std::unique_ptr<fastjet::ClusterSequence> parton_cs_, hadron_cs_;

  // jet matching
  JetMatcher parton_hadron_matcher, hard_parton_matcher;
  std::vector<int> parton_hard_inds_, hadron_hard_inds_, 
                   hadron_parton_inds_, parton_hadron_inds_;

///////////////////////////////////////////////////////////////////////////////
// PUBLIC METHODS
///////////////////////////////////////////////////////////////////////////////
public:

  // constructor from command line and optional comments to add to description
  EventGenerator(int argc, char** argv, std::string comments = "") :
    user_comments(comments),
    start_time_(std::chrono::steady_clock::now()),
    force_hadron_level_failed_(0),
    write_counter_(0),
    pythia_("../share/Pythia8/xmldoc", false)
  { parse(argc, argv); }

  // parses command line input into internal vairables
  void parse(int argc, char** argv);

  // method to advance event generator, returns how many units were obtained
  unsigned next();

  // determines if unit j was fully matched
  bool unit_fully_matched(unsigned j) const {
    assert(j < matched_units.size());
    return matched_units[j].status & M_FULLY_MATCHED;
  }

  // accesses the hadron jet in unit j
  std::pair<PseudoJet,bool> get_hadron_jet(unsigned j) const {
    int hadron_jet_ind(matched_units[j].hadron_jet_ind);
    if (hadron_jet_ind == -1) return std::make_pair(PseudoJet(false), false);
    return std::make_pair(hadron_jets[hadron_jet_ind], true);
  }

  // accesses the parton jet in unit j
  std::pair<PseudoJet,bool> get_parton_jet(unsigned j) const {
    int parton_jet_ind(matched_units[j].parton_jet_ind);
    if (parton_jet_ind == -1) return std::make_pair(PseudoJet(false), false);
    return std::make_pair(parton_jets[parton_jet_ind], true);
  }

  // sets per-jet features determined by the user, cleared after each call to write_unit
  void set_hadron_jet_features(const std::vector<double> & features) { 
    hadron_jet_features_ = features;
  }
  void set_parton_jet_features(const std::vector<double> & features) {
    parton_jet_features_ = features;
  }

  // writes specified unit to output file
  void write_unit(unsigned j);

  // prints update to std::cout
  void print(unsigned i) const {
    std::ostringstream oss;
    oss << std::setprecision(6)
        << process_ 
        << ", seed " << seed_ 
        << ": Produced " << i 
        << " events, " << double(i)/ntot*100
        << "% done in " << std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time_).count()
        << "s";
    std::cout << oss.str() << std::endl;
  }

  // writes a summary of the event generation, to be called at the end of the run
  void summarize() {
    summarize(outfile_);
    if (verbose && !outfilepath_.empty()) summarize(std::cout);
  }
  void summarize(std::ostream &);

///////////////////////////////////////////////////////////////////////////////
// PRIVATE METHODS
///////////////////////////////////////////////////////////////////////////////
private:

  // initialization functions
  void setup_process();
  void setup_jet_alg();
  void setup_jet_matchers();
  void setup_pythia();

  // print descrption of options
  void describe(std::ostream &) const;

  // store particles for fastjet from a pythia event record
  void store_hardprocess();
  void p2fj(std::vector<PseudoJet> & pjs, std::vector<int> & pids);

  // accesses vectors of matched indices, preserving -1 values
  int hadron_hard_ind(int j) const {
    return j != -1 && unsigned(j) < hadron_hard_inds_.size() ? hadron_hard_inds_[j] : -1;
  }
  int hadron_parton_ind(int j) const {
    return j != -1 && unsigned(j) < hadron_parton_inds_.size() ? hadron_parton_inds_[j] : -1;
  }
  int parton_hadron_ind(int j) const {
    return j != -1 && unsigned(j) < parton_hadron_inds_.size() ? parton_hadron_inds_[j] : -1;
  }
  int parton_hard_ind(int j) const {
    return j != -1 && unsigned(j) < parton_hard_inds_.size() ? parton_hard_inds_[j] : -1;
  }

  // writes hard process particle to file
  void write_hard(int hard_ind, double ref_phi, const char * id);

  // writes jet to file
  void write_jet(const PseudoJet & jet, double ref_phi, char id, const std::vector<double> & features);

  // writes constituents to file (ignores within-R info)
  void write_constituents(const PseudoJet & jet, const std::vector<int> & pids, double ref_phi, char id);

  // writes jet constituents and particles within R of other object to file
  void write_constituents_and_within_R(const PseudoJet & primary_jet, const PseudoJet & secondary_jet,
                                       const std::vector<PseudoJet> & primary_particles,
                                       const std::vector<int> & pids, double ref_phi, char id);

  // write particles that are within R of given jet, but only those not previously written
  void write_within_R(const PseudoJet & jet,
                      const std::vector<PseudoJet> & particles,
                      const std::vector<int> & pids,
                      const std::vector<char> & writes,
                      double ref_phi, char id);

  // write all particles that are within R of given jet
  void write_within_R(const PseudoJet & jet,
                      const std::vector<PseudoJet> & particles,
                      const std::vector<int> & pids,
                      double ref_phi, char id);

  // write various kinematics to file 
  void write_ptyphimid(const PseudoJet & pj, double ref_phi);
  void write_ptyphim(const PseudoJet & pj, double ref_phi, bool round_mass = false);
  void write_ptyphipid(const PseudoJet & pj, double ref_phi, int pid);
  void write_ptyphi(const PseudoJet & pj, double ref_phi);

}; // EventGenerator

} // namespace mc

#endif // EventGenerator_HH
