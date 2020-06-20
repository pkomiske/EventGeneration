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

// C++ standard library
#include <exception>
#include <limits>
#include <stdexcept>
#include <tuple>

// boost
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
namespace po = boost::program_options;

#include "EventGenerator.hh"
#include "fastjet/version.hh"

namespace mc {

// namespace shortenings
using fastjet::ClusterSequence;

///////////////////////////////////////////////////////////////////////////////
// EventGenerator - Public methods 
///////////////////////////////////////////////////////////////////////////////

// parses command line input into internal vairables
void EventGenerator::parse(int argc, char** argv) {

  // reset internals
  argc_ = argc;
  argv_ = argv;
  exit = false;

  // to catch bad input
  try {

    po::variables_map vm;
    po::options_description od("Options for event generation with Pythia8.");
    od.add_options()

    // important options (will also be positional)
    ("outfilepath,o", po::value<std::string>(&outfilepath_)->default_value(""),
                        "Path to the output file.")
    ("seed,s",        po::value<int>(&seed_)->default_value(-1),
                        "Random seed to use.")
    ("process,p",     po::value<std::string>(&process_)->default_value("QCD"),
                        "Name of process to generate: [Zg, Zq, Zjet, QCD, W, Top].")
    ("ntot,n",        po::value<uint> (&ntot)->default_value(1),
                        "Number of events to generate.")

    // event generation options
    ("ecm",            po::value<double>(&ecm_)->default_value(14000),
                         "Center of mass energy in GeV.")
    ("pthatmin",       po::value<double>(&pthatmin_)->default_value(-1),
                         "Pythia pthatmin parameter.")
    ("pthatmax",       po::value<double>(&pthatmax_)->default_value(-1),
                         "Pythia pthatmax parameter.")
    ("biaspow,b",      po::value<double>(&biaspow_)->default_value(-1),
                         "Pythia bias power.")
    ("no-ue",          po::bool_switch(&no_ue_),
                         "Turns off underlying event.")
    ("no_isr",         po::bool_switch(&no_isr_),
                         "Turns off initial-state radiation.")
    ("no_fsr",         po::bool_switch(&no_fsr_),
                         "Turns off final-state radiation.")
    ("no-hardprocess", po::bool_switch(&no_hardprocess_),
                         "Turns off outputting the hard process.")
    ("no-partonlevel", po::bool_switch(&no_partonlevel_),
                         "Turns off outputing parton-level events.")
    ("no-hadronlevel", po::bool_switch(&no_hadronlevel_),
                         "Turns off outputing hadron-level events.")

    // jet options
    ("jet-R,R",     po::value<double>(&jet_R)->default_value(1.0),
                      "Jet radius to use.")
    ("ptjetmin",    po::value<double>(&ptjetmin_)->default_value(500),
                      "Minimum jet pt.")
    ("ptjetmax",    po::value<double>(&ptjetmax_)->default_value(-1),
                      "Maximum jet pt.")
    ("absrapmax,y", po::value<double>(&absrapmax_)->default_value(5),
                      "Maximum absolute jet rapidity.")
    ("jet-alg,a",   po::value<std::string>(&jet_alg_name_)->default_value("antikt"),
                      "Jet algorithm to use: [akt, kt, ca]")

    // jet matching options
    ("hard-parton-match",   po::value<double>(&hard_parton_matcher.deltaRMax())->default_value(1),
                              "Matching factor between hard process particle and parton-level jet.")
    ("parton-hadron-match", po::value<double>(&parton_hadron_matcher.deltaRMax())->default_value(1),
                              "Matching factor between parton-level and hadron-level jets.")

    // printing/outputting options
    ("print-pythia",              po::bool_switch(&print_pythia_),
                                    "Prints Pythia output.")
    ("print-every,e",             po::value<uint>(&print_every)->default_value(1000),
                                    "Printing frequency.")
    ("num-print-event",           po::value<uint>(&num_print_event_)->default_value(0),
                                    "Number of Pythia events to print.")
    ("num-print-proc",            po::value<uint>(&num_print_proc_)->default_value(0),
                                    "Number of Pythia hard processses to print.")
    ("output-fully-matched-only", po::bool_switch(&fully_matched_only),
                                    "Output only fully matched events.")
    ("output-within-R",           po::bool_switch(&output_within_R_),
                                    "Outputs all partons/hadrons within R of matched object.")
    ("outprecision", po::value<uint>(&outprecision_)->default_value(12),
                                    "Number of digits to print for floats.")

    // generic options
    ("help,h",      "Produce help message.")
    ("verbose,v", po::bool_switch(&verbose),
                    "Turns on other verbose printing in the program.")
    ; // ends add_options()

    // allow important arguments to be provided as positional arguments
    po::positional_options_description pod;
    pod.add("outfilepath", 1);
    pod.add("seed", 1);
    pod.add("process", 1);
    pod.add("ntot", 1);
    
    // actually parse command line
    po::store(po::command_line_parser(argc_, argv_).options(od).positional(pod).run(), vm);

    // print help if requested
    if (vm.count("help")) {
      exit = true;
      std::cout << "Usage: " << argv_[0];
      for (uint i = 0; i < pod.max_total_count(); i++)
        std::cout << ' ' << pod.name_for_position(i);
      std::cout << "\n\n" << od << '\n';
      return;
    }

    // check things
    po::notify(vm);

    // setup event generation
    setup_process();
    setup_jet_alg();
    setup_jet_matchers();
    setup_pythia();

    // setup output file (will be std::cout if no filepath provided)
    if (outfilepath_.empty()) {
      outfile_.std::ostream::rdbuf(std::cout.rdbuf());
      describe(outfile_);
    }
    else {
      outfile_.open(outfilepath_);
      describe(outfile_);
      if (verbose) describe(std::cout);
    }

    // set output precision
    outfile_ << std::setprecision(outprecision_);
  }

  // catch any exceptions (from po::notify)
  catch (std::exception & e) {
    std::cerr << e.what() << '\n';
    throw e;
  }
}

// method to advance event generator, returns how many units were obtained
uint EventGenerator::next() {

  // iterate until we get an event with some jets
  while (true) {

    // generate process with pythia
    if (!pythia_.next()) continue;

    // get hard process info
    store_hardprocess();

    // store partons
    bool have_any_jets(false);
    if (partonlevel) {

      // examines event record and stores partons
      p2fj(partons, parton_pids);

      // findJunctions = false in call to forceHadronLevel since Pythia will have set that up already
      if (hadronlevel && !pythia_.forceHadronLevel(false)) {
        force_hadron_level_failed_++;
        continue;
      }

      // get parton jets (the unique pointer deletes the last cluster sequence)
      parton_cs_ = std::unique_ptr<ClusterSequence>(new ClusterSequence(partons, jet_def_));
      parton_jets = jet_sel_(parton_cs_->inclusive_jets(pthatmin_));
      if (parton_jets.size()) {
        have_any_jets = true;

        // match to the hard process
        if (hardprocess)
          std::tie(parton_hard_inds_, std::ignore) = hard_parton_matcher(origs, parton_jets);
      }
    }

    // store hadrons
    if (hadronlevel) {

      // examines event record and stores hadrons
      p2fj(hadrons, hadron_pids);

      // get hadron jets
      hadron_cs_ = std::unique_ptr<ClusterSequence>(new ClusterSequence(hadrons, jet_def_));
      hadron_jets = jet_sel_(hadron_cs_->inclusive_jets(pthatmin_));
      if (hadron_jets.size()) {
        have_any_jets = true;

        // match to hard process or parton level if applicable
        if (hardprocess)
          std::tie(hadron_hard_inds_, std::ignore) = hard_parton_matcher(origs, hadron_jets);
        if (partonlevel)
          std::tie(hadron_parton_inds_, parton_hadron_inds_) = parton_hadron_matcher(parton_jets, hadron_jets);
      }
    }

    if (have_any_jets) break;
  }

  // clear units
  matched_units.clear();

  // iterate over hadron jets
  for (uint j = 0; j < hadron_jets.size(); j++) {

    // setup new matched unit
    matched_units.emplace_back();
    MatchedUnit & unit(matched_units.back());

    unit.phi = hadron_jets[j].phi();
    unit.hadron_jet_ind = j;
    unit.parton_jet_ind = hadron_parton_ind(j);
    unit.hard_hadron_ind = hadron_hard_ind(j);
    unit.hard_parton_ind = parton_hard_ind(unit.parton_jet_ind);

    // determine status of this set of inds
    unsigned char status(base_status_ | M_HADRON_JET);
    if (partonlevel) { // we are considering parton-level info
      if (unit.parton_jet_ind != -1) { // we expect a parton-level jet to have matched
        status |= M_PARTON_JET; // we have a parton-level jet
        if (hardprocess){ // we are considering hardprocess info
          if (unit.hard_parton_ind != -1) status |= M_PARTON_HARD; // hard process matched parton jet
          else status &= M_PARTIAL_MATCH; // incomplete because hard process did not match parton jet
        }
      }
      else status &= M_PARTIAL_MATCH; // incomplete because no parton-level jet
    }

    if (hardprocess) { // we are considering hard process info
      if (unit.hard_hadron_ind != -1) status |= M_HADRON_HARD; // hard process matched hadron jet
      else status &= M_PARTIAL_MATCH; // incomplete because hard process did not match hadron jet
    }

    // check if parton and hadron jets matched the same hard process particle
    if ((status & M_PARTON_HARD) && (status & M_HADRON_HARD)) {
      if (unit.hard_parton_ind != unit.hard_hadron_ind) status &= M_PARTIAL_MATCH;
      else status |= M_PH_SAME_HARDS;
    }

    unit.status = status;
  }

  // iterate over parton jets
  for (uint j = 0; j < parton_jets.size(); j++) {

    // check that we didn't already output this along with a hadron jet
    if (parton_hadron_ind(j) != -1) continue;

    // setup new matched unit
    matched_units.emplace_back();
    MatchedUnit & unit(matched_units.back());

    unit.phi = parton_jets[j].phi();
    unit.parton_jet_ind = j;
    unit.hard_parton_ind = parton_hard_ind(j);
    unit.hadron_jet_ind = parton_hadron_ind(j); // should always be -1 here
    unit.hard_hadron_ind = hadron_hard_ind(unit.hadron_jet_ind); // should always be -1 here
    assert(unit.hadron_jet_ind == -1 && unit.hard_hadron_ind == -1);

    unsigned char status(base_status_ | M_PARTON_JET);
    if (hadronlevel) status &= M_PARTIAL_MATCH; // here we have no associated hadron jet, so incomplete if expecting that
    if (hardprocess) { // we are considering hardprocess info
      if (unit.hard_parton_ind != -1) status |= M_PARTON_HARD; // hard process matched parton jet
      else status &= M_PARTIAL_MATCH; // incomplete because hard process did not match parton jet
    }

    unit.status = status;
  }

  return matched_units.size();
}

// writes unit to output file
void EventGenerator::write_unit(uint j) {

  // check for appropriate input and acquire unit
  assert(j < matched_units.size());
  const MatchedUnit & unit(matched_units[j]);
  const unsigned char status(unit.status);

  // output unit index if provided, and status of unit
  outfile_ << "i " << write_counter_++      << '\n'
           << "S " << int(status)           << '\n'
           << "W " << pythia_.info.weight() << '\n';

  // output single hard process if joint, else just the parton-matched one
  if (status & M_PH_SAME_HARDS)
    write_hard(unit.hard_hadron_ind, unit.phi, "O");
  else {
    if (status & M_HADRON_HARD)
      write_hard(unit.hard_hadron_ind, unit.phi, "OH");
    if (status & M_PARTON_HARD)
      write_hard(unit.hard_parton_ind, unit.phi, "OP");
  }

  // write hadron jet if we have it
  if (status & M_HADRON_JET)
    write_jet(hadron_jets[unit.hadron_jet_ind], unit.phi, 'H', hadron_jet_features_);

  // write parton jet if we have it
  if (status & M_PARTON_JET)
    write_jet(parton_jets[unit.parton_jet_ind], unit.phi, 'P', parton_jet_features_);

  // write hadrons
  if (status & M_HADRON_LEVEL) { // we have hadrons at all
    if (status & M_HADRON_JET) { // we have a hadron jet
      if (status & M_PARTON_JET) // we also have a parton jet
        write_constituents_and_within_R(hadron_jets[unit.hadron_jet_ind], parton_jets[unit.parton_jet_ind], 
                                        hadrons, hadron_pids, unit.phi, 'H');
      else // no parton jet, just output hadron constituents
        write_constituents(hadron_jets[unit.hadron_jet_ind], hadron_pids, unit.phi, 'H');
    }
    else { // no hadron jet, output hadrons within R of parton jet
      assert(status & M_PARTON_JET); // should always be true since we should have either a parton or hadron jet
      write_within_R(parton_jets[unit.parton_jet_ind], hadrons, hadron_pids, unit.phi, 'H');
    } 
  }

  // write partons
  if (status & M_PARTON_LEVEL) { // we have partons at all
    if (status & M_PARTON_JET) { // we have a parton jet
      if (status & M_HADRON_JET) // we also have a hadron jet
        write_constituents_and_within_R(parton_jets[unit.parton_jet_ind], hadron_jets[unit.hadron_jet_ind], 
                                        partons, parton_pids, unit.phi, 'P');
      else // no hadron jet, just output parton constituents
        write_constituents(parton_jets[unit.parton_jet_ind], parton_pids, unit.phi, 'P');
    }
    else { // no parton jet, output partons within R of hadron jet
      assert(status & M_HADRON_JET); // should always be true since we should have either a parton or hadron jet
      write_within_R(hadron_jets[unit.hadron_jet_ind], partons, parton_pids, unit.phi, 'P');
    } 
  }

  // newline between events
  outfile_ << '\n';

  // clear features
  hadron_jet_features_.clear();
  parton_jet_features_.clear();
}

// writes a summary of the event generation, to be called at the end of the run
void EventGenerator::summarize(std::ostream & os) {

  // summarize jet matching
  if (hardprocess) {
    hard_parton_matcher.summarize(os, "# ");
    os << "#\n";
  }
  if (partonlevel && hadronlevel) {
    parton_hadron_matcher.summarize(os, "# ");
    os << "#\n";
  }

  // summarize forceHadronLevel
  if (hadronlevel)
    os << "# number of forceHadronLevel fails: " << force_hadron_level_failed_ << "\n#\n";

  // output cross section info
  auto time_diff(std::chrono::steady_clock::now() - start_time_);
  os << "# CrossSection: " << pythia_.info.sigmaGen()*1000000 << " +- " 
                           << pythia_.info.sigmaErr()*1000000 << " nb \n"
     << "# TotalWeight: "  << pythia_.info.weightSum() << '\n'
     << "# Duration: "     << std::chrono::duration<double>(time_diff).count() << "s\n";
}


///////////////////////////////////////////////////////////////////////////////
// EventGenerator - Private methods 
///////////////////////////////////////////////////////////////////////////////

// translates options to pythia settings
void EventGenerator::setup_process() {

  if (no_partonlevel_ && no_hadronlevel_)
    throw std::runtime_error("parton-level and hadron-level cannot both be off");
  hardprocess = !no_hardprocess_;
  partonlevel = !no_partonlevel_;
  hadronlevel = !no_hadronlevel_;

  boost::to_lower(process_);
  if (process_ == "zq") {
    description_ = "Z(->nunubar)q";
    hardproc_ids_ = {6};
    pythia_commands_ = {"WeakBosonAndParton:qg2gmZq = on",
                        "WeakZ0:gmZmode = 2",
                        "23:onMode = off",
                        "23:onIfAny = 12 14 16"};
  }
  else if (process_ == "zg") {
    description_ = "Z(->nunubar)g";
    hardproc_ids_ = {6};
    pythia_commands_ = {"WeakBosonAndParton:qqbar2gmZg = on",
                        "WeakZ0:gmZmode = 2",
                        "23:onMode = off",
                        "23:onIfAny = 12 14 16"};
  }
  else if (process_ == "zjet") {
    description_ = "Z(->nunubar)jet";
    hardproc_ids_ = {6};
    pythia_commands_ = {"WeakBosonAndParton:qg2gmZq = on",
                        "WeakBosonAndParton:qqbar2gmZg = on",
                        "WeakZ0:gmZmode = 2",
                        "23:onMode = off",
                        "23:onIfAny = 12 14 16"};
  }
  else if (process_ == "qcd") {
    description_ = "HardQCD:all";
    hardproc_ids_ = {5, 6};
    pythia_commands_ = {"HardQCD:all = on"};
  }
  else if (process_ == "w") {
    description_ = "W+W-->qqbar,qqbar";
    hardproc_ids_ = {5, 6};
    pythia_commands_ = {"WeakDoubleBoson:ffbar2WW = on",
                        "24:onMode = off",
                        "24:onIfAny = 1 2 3 4 5"};
  }
  else if (process_ == "top") {
    description_ = "ttbar->W+(->qqbar)b,W-(->qqbar)bbar";
    hardproc_ids_ = {5, 6};
    pythia_commands_ = {"Top:gg2ttbar = on",
                        "Top:qqbar2ttbar = on",
                        "24:onMode = off",
                        "24:onIfAny = 1 2 3 4 5"};
  }
  else throw std::runtime_error("Unknown process " + process_);

  nmaxjets = hardproc_ids_.size();
  hadron_hard_inds_ = parton_hard_inds_ = hadron_parton_inds_ = std::vector<int>(nmaxjets, -1);

  // the base status for a jet with the specified options
  base_status_ = M_FULLY_MATCHED;
  if (hadronlevel) base_status_ |= M_HADRON_LEVEL;
  if (partonlevel) base_status_ |= M_PARTON_LEVEL;
}

// translates options to jet selection
void EventGenerator::setup_jet_alg() {

  boost::to_lower(jet_alg_name_);
  if (jet_alg_name_ == "antikt" || jet_alg_name_ == "akt" || jet_alg_name_ == "anti-kt")
    jet_alg_ = fastjet::antikt_algorithm;
  else if (jet_alg_name_ == "kt")
    jet_alg_ = fastjet::kt_algorithm;
  else if (jet_alg_name_ == "ca")
    jet_alg_ = fastjet::cambridge_algorithm;
  else throw std::runtime_error("Unrecognized jet algorithm " + jet_alg_name_);

  double actual_ptjetmax(ptjetmax_ == -1 ? std::numeric_limits<double>::infinity() : ptjetmax_);
  jet_def_ = fastjet::JetDefinition(jet_alg_, jet_R);
  jet_sel_ = fastjet::SelectorNHardest(nmaxjets) &&
             fastjet::SelectorAbsRapMax(absrapmax_) &&
             fastjet::SelectorPtRange(ptjetmin_, actual_ptjetmax);
  ClusterSequence::set_fastjet_banner_stream(new std::ostringstream());
}

// convolves options with the JetMatcher objects that have been set up
void EventGenerator::setup_jet_matchers() {
  hard_parton_matcher.name() = "hard-parton-" + std::to_string(seed_);
  parton_hadron_matcher.name() = "parton-hadron-" + std::to_string(seed_);
  hard_parton_matcher.deltaRMax() *= jet_R;
  parton_hadron_matcher.deltaRMax() *= jet_R;
}

// feeds commands to Pythia8 object
void EventGenerator::setup_pythia() {
  pythia_.settings.mode("Beams:idA", 2212);
  pythia_.settings.mode("Beams:idB", 2212);
  pythia_.settings.mode("Beams:eCM", ecm_);

  pythia_.settings.flag("Random:setSeed", true);
  pythia_.settings.mode("Random:seed", seed_);

  pythia_.settings.flag("PartonLevel:MPI", !no_ue_);
  pythia_.settings.flag("PartonLevel:ISR", !no_isr_);
  pythia_.settings.flag("PartonLevel:FSR", !no_fsr_);
  pythia_.settings.flag("HadronLevel:all", no_partonlevel_);

  // take 25% wide pthat[min/max] if not provided
  if (pthatmin_ == -1) pthatmin_ = 0.75*ptjetmin_;
  if (pthatmax_ == -1 && ptjetmax_ != -1) pthatmax_ = 1.25*ptjetmax_;
  pythia_.settings.parm("PhaseSpace:pTHatMin", pthatmin_);
  pythia_.settings.parm("PhaseSpace:pTHatMax", pthatmax_);
  pythia_.settings.flag("PhaseSpace:bias2Selection", biaspow_ > 0);
  pythia_.settings.parm("PhaseSpace:bias2electionPow", biaspow_);

  pythia_.settings.flag("Print:quiet", !print_pythia_);
  pythia_.settings.mode("Next:numberShowEvent", num_print_event_);
  pythia_.settings.mode("Next:numberShowProcess", num_print_proc_);

  // iterate over process-level commands
  for (std::vector<std::string>::iterator it = pythia_commands_.begin();
       it != pythia_commands_.end(); ++it)
    pythia_.readString(*it);

  pythia_.init();
}

// writes a nicely formatted description of all options to the provided stream
void EventGenerator::describe(std::ostream & os) const {

  os << "# invocation:";
  for (int i = 0; i < argc_; i++) os << ' ' << argv_[i];

  time_t t(time(NULL));
  os << '\n'
     << "# date: "    << asctime(localtime(&t))
     << "# pythia: "  << PYTHIA_VERSION           << '\n'
     << "# fastjet: " << fastjet::fastjet_version << '\n'
     << "#\n"
     << "# outfilepath: " << outfilepath_ << '\n'
     << "# seed: "        << seed_        << '\n'
     << "# print_every: " << print_every  << '\n'
     << "#\n"
     << "# process: "     << process_     << '\n'
     << "# description: " << description_ << '\n'
     << "# ntot: "        << ntot         << '\n'
     << "# ecm: "         << ecm_         << '\n'
     << "# pthatmin: "    << pthatmin_    << '\n'
     << "# pthatmax: "    << pthatmax_    << '\n'
     << "# biaspow: "     << (biaspow_ > 0 ? std::to_string(biaspow_) : "off") << '\n'
     << "#\n"
     << "# ue: "  << (no_ue_ ? "off" : "on")  << '\n'
     << "# isr: " << (no_isr_ ? "off" : "on") << '\n'
     << "# fsr: " << (no_fsr_ ? "off" : "on") << '\n'
     << "# partonlevel: " << (no_partonlevel_ ? "no" : "yes") << '\n'
     << "# hadronlevel: " << (no_hadronlevel_ ? "no" : "yes") << '\n'
     << "#\n"
     << "# jet_alg: "   << jet_alg_name_ << '\n'
     << "# jet_R: "     << jet_R         << '\n'
     << "# ptjetmin_: " << ptjetmin_     << '\n'
     << "# ptjetmax_: " << ptjetmax_     << '\n'
     << "# absrapmax: " << absrapmax_    << '\n'
     << "#\n"
     << "# hard_parton_deltaRMax: "   << hard_parton_matcher.deltaRMax()   << '\n'
     << "# parton_hadron_deltaRMax: " << parton_hadron_matcher.deltaRMax() << '\n'
     << "# fully_matched_only: "           << (fully_matched_only ? "yes" : "no")    << '\n'
     << "# output_within_R: "         << (output_within_R_ ? "yes" : "no") << '\n'
     << "#\n"
     << "# status_bits:\n"
     << "#   M_FULLY_MATCHED : " << int(M_FULLY_MATCHED) << '\n'
     << "#   M_HADRON_LEVEL  : " << int(M_HADRON_LEVEL)  << '\n'
     << "#   M_PARTON_LEVEL  : " << int(M_PARTON_LEVEL)  << '\n'
     << "#   M_HADRON_JET    : " << int(M_HADRON_JET)    << '\n'
     << "#   M_HADRON_HARD   : " << int(M_HADRON_HARD)   << '\n'
     << "#   M_PARTON_JET    : " << int(M_PARTON_JET)    << '\n'
     << "#   M_PARTON_HARD   : " << int(M_PARTON_HARD)   << '\n'
     << "#   M_PH_SAME_HARDS : " << int(M_PH_SAME_HARDS) << '\n'
     << "#\n"
     << "# format of matched units:\n"
     << "#   i       unit index within this file\n"
     << "#   S       status, see status_bits above\n"
     << "#   W       event-generator weight\n"
     << "#   O[H/P]  pt y phi m pid : hard-process particle, [H/P] indicates specific matching\n"
     << "#   D       pt y phi m pid : decay product from the above hard-process particle\n"
     << "#   J(H/P)  pt y phi m [user_values] : jet, hadron-level or parton-level, user provided features\n"
     << "#   H[C][R] pt y phi pid : hadron, jet constituent (C), within jet_R of parton jet (R)\n"
     << "#   P[C][R] pt y phi pid : parton, jet constituent (C), within jet_R of hadron jet (R)\n";

  if (!user_comments.empty())
    os << "#\n" << user_comments;

  os << '\n'; // newline to indicate start of events
}

// stores outgoing particles and any decays of theirs from the hard process
inline void EventGenerator::store_hardprocess() {
  origs.clear();
  decays.clear();
  for (uint i : hardproc_ids_) {
    origs.emplace_back(pythia_.process[i].p());
    origs.back().set_user_index(pythia_.process[i].id());

    decays.push_back(std::vector<PseudoJet>());
    for (int j : pythia_.process[i].daughterListRecursive()) {
      decays.back().emplace_back(pythia_.process[j].p());
      decays.back().back().set_user_index(pythia_.process[j].id());
    }
  }
}

// store particles for fastjet from a pythia event record
inline void EventGenerator::p2fj(std::vector<PseudoJet> & pjs, std::vector<int> & pids) {
  pjs.clear();
  pids.clear();
  for (int i = 0; i < pythia_.event.size(); i++) {
    const Pythia8::Particle & particle(pythia_.event[i]);
    
    // always ignore neutrinos
    int idabs  = particle.idAbs();
    if (particle.isFinal() && idabs != 12 && idabs != 14 && idabs != 16) {
      pjs.emplace_back(particle.p());

      // user index will be its position in the ids array
      pjs.back().set_user_index(pids.size());
      pids.push_back(particle.id());
    }
  }
}

// writes hard process particle to file
inline void EventGenerator::write_hard(int hard_ind, 
                                       double ref_phi,
                                       const char * id) {
  outfile_ << id << ' ';
  write_ptyphimid(origs[hard_ind], ref_phi);
  outfile_ << '\n';
  for (const PseudoJet & decay : decays[hard_ind]) {
    outfile_ << "D ";
    write_ptyphimid(decay, ref_phi);
    outfile_ << '\n';
  }
}

// writes jet to file
inline void EventGenerator::write_jet(const PseudoJet & jet,
                                      double ref_phi, char id, 
                                      const std::vector<double> & features) {
  outfile_ << 'J' << id << ' ';
  write_ptyphim(jet, ref_phi);
  if (!features.empty())
    for (double feature : features)
      outfile_ << ' ' << feature;
  outfile_ << '\n';
}

// writes constituents to file (ignores within-R info)
inline void EventGenerator::write_constituents(const PseudoJet & jet,
                                               const std::vector<int> & pids,
                                               double ref_phi, char id) {
  for (const PseudoJet & pj : jet.constituents()) {
    outfile_ << id << "C ";
    write_ptyphipid(pj, ref_phi, pids[pj.user_index()]);
    outfile_ << '\n';
  }
}

// writes jet constituents and particles within R of other object to file
inline void
EventGenerator::write_constituents_and_within_R(const PseudoJet & primary_jet,
                                                const PseudoJet & secondary_jet,
                                                const std::vector<PseudoJet> & primary_particles,
                                                const std::vector<int> & pids,
                                                double ref_phi, char id) {

  // vector of bools to keep track of which particles have been output
  std::vector<char> writes(pids.size(), false);

  // output constituents from primary jet
  for (const PseudoJet & pj : primary_jet.constituents()) {
    outfile_ << id << 'C';

    // indicate that we've looked through this one
    uint index(pj.user_index());
    writes[index] = true;

    // test for const within R of seconday jet
    if (output_within_R_ && secondary_jet.delta_R(pj) <= jet_R)
      outfile_ << 'R';

    outfile_ << ' ';
    write_ptyphipid(pj, ref_phi, pids[index]);
    outfile_ << '\n';
  }

  // output particles that are within R of seconday jet
  write_within_R(secondary_jet, primary_particles, pids, writes, ref_phi, id);
}

// write particles that are within R of given jet, but only those not previously written
inline void EventGenerator::write_within_R(const PseudoJet & jet,
                                           const std::vector<PseudoJet> & particles,
                                           const std::vector<int> & pids,
                                           const std::vector<char> & writes,
                                           double ref_phi, char id) {

  if (!output_within_R_) return;

  for (const PseudoJet & pj : particles) {
    uint index(pj.user_index());
    if (!writes[index] && jet.delta_R(pj) <= jet_R) {
      outfile_ << id << "R ";
      write_ptyphipid(pj, ref_phi, pids[index]);
      outfile_ << '\n';
    }
  }
}

// write all particles that are within R of given jet
inline void EventGenerator::write_within_R(const PseudoJet & jet,
                                           const std::vector<PseudoJet> & particles,
                                           const std::vector<int> & pids,
                                           double ref_phi, char id) {

  if (!output_within_R_) return;

  for (const PseudoJet & pj : particles) {
    if (jet.delta_R(pj) <= jet_R) {
      outfile_ << id << "R ";
      write_ptyphipid(pj, ref_phi, pids[pj.user_index()]);
      outfile_ << '\n';
    }
  }
}

inline void EventGenerator::write_ptyphimid(const PseudoJet & pj, double ref_phi) {
  write_ptyphim(pj, ref_phi, true);
  outfile_ << ' ' << pj.user_index();
}

inline void EventGenerator::write_ptyphim(const PseudoJet & pj,
                                          double ref_phi, bool round_mass) {
  write_ptyphi(pj, ref_phi);
  outfile_ << ' ';
  if (round_mass) outfile_ << m_round(pj.m());
  else outfile_ << pj.m();
}

inline void EventGenerator::write_ptyphipid(const PseudoJet & pj,
                                            double ref_phi, int pid) {
  write_ptyphi(pj, ref_phi);
  outfile_ << ' ' << pid;
}

inline void EventGenerator::write_ptyphi(const PseudoJet & pj, double ref_phi) {
  outfile_ << pj.pt() << ' '
           << pj.rap() << ' '
           << phi_fix(pj.phi(), ref_phi);
}

} // namespace mc