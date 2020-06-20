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

#ifndef JetMatcher_HH
#define JetMatcher_HH

// C++ standard library
#include <string>
#include <utility>
#include <vector>

// FastJet
#include "fastjet/PseudoJet.hh"

namespace mc {

typedef unsigned int uint;
using fastjet::PseudoJet;

///////////////////////////////////////////////////////////////////////////////
// JetMatcher - Associates jets from two related sources
///////////////////////////////////////////////////////////////////////////////
class JetMatcher {
private:

  double deltaRMax_;
  std::string name_;

  // counters
  uint numNewMatchedMultipleTimes, numRefMatchedMultipleTimes,
       numHardestNewFailMatch, numHardestRefFailMatch,
       numHardestRefMatchedNotHardestNew, nTotal;

public:

  // constructors
  JetMatcher() : JetMatcher(0, "") {}
  JetMatcher(double deltaRMax, std::string name) :
    deltaRMax_(deltaRMax), name_(name),

    // counters
    numNewMatchedMultipleTimes(0), numRefMatchedMultipleTimes(0),
    numHardestNewFailMatch(0), numHardestRefFailMatch(0),
    numHardestRefMatchedNotHardestNew(0), nTotal(0)
  {}

  // access parameters
  double & deltaRMax() { return deltaRMax_; }
  const double & deltaRMax() const { return deltaRMax_; }
  std::string & name() { return name_; }

  // function that actual determines if two jets match
  // this can be updated by derived classes
  virtual bool match_specific_jets(const PseudoJet & jet1, const PseudoJet & jet2) const {
    
    // delta R matching
    return jet1.delta_R(jet2) <= deltaRMax_;
  }

  // match collections of jets
  std::pair<std::vector<int>, std::vector<int>> 
  operator()(const std::vector<PseudoJet> & newJets, const std::vector<PseudoJet> & refJets);

  // summarize the statistics of the jet matcher
  void summarize(std::ostream &, const char * = "") const;

}; // JetMatcher

} // namespace mc

#endif // JetMatcher_HH
