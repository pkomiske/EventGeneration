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
#include <algorithm>
#include <iomanip>
#include <sstream>

#include "JetMatcher.hh"

namespace mc {

std::pair<std::vector<int>, std::vector<int>> 
JetMatcher::operator()(const std::vector<PseudoJet> & newJets,
                       const std::vector<PseudoJet> & refJets) {

  const uint sizeNew(newJets.size()), sizeRef(refJets.size());

  // containers of indices with the mapping between new and old jets
  std::vector<int> new2ref(sizeNew, -1), ref2new(sizeRef, -1);

  // if we have no new or no ref jets, return
  if (sizeNew == 0 || sizeRef == 0) return std::make_pair(ref2new, new2ref);

  for (uint j = 0; j < sizeRef; j++) {
    const PseudoJet & refJet(refJets[j]);

    for (uint i = 0; i < sizeNew; i++) {
      const PseudoJet & newJet(newJets[i]);

      // if new matches ref and ref has not previously matched
      if (match_specific_jets(newJet, refJet)) {
        if (ref2new[j] == -1) {
          if (new2ref[i] == -1) {
            new2ref[i] = j;
            ref2new[j] = i;
            break;
          }
          else numNewMatchedMultipleTimes++;
        }
        else numRefMatchedMultipleTimes++;
      }
    }
  }

  // check that highest two new jets matched
  for (uint i = 0, m = std::min<uint>(2, sizeNew); i < m; i++) {
    if (new2ref[i] == -1) numHardestNewFailMatch++;
  }

  // check that highest two ref jets matched and get matching new jets
  for (uint i = 0, m = std::min<uint>(2, sizeRef); i < m; i++) {
    int new_i(ref2new[i]);

    if (new_i == -1) numHardestRefFailMatch++;
    else if (new_i >= 2) numHardestRefMatchedNotHardestNew++;
  }
  nTotal += sizeRef;

  return std::make_pair(ref2new, new2ref);
}

// summarize
void JetMatcher::summarize(std::ostream & os, const char * pre) const {
  std::ostringstream oss;
  oss << std::setprecision(6);
  oss << pre << std::setw(name_.size()) << ' '
      << "    hRefNoMatch    hNewNoMatch     hRef!~hNew  nRefMatchMany  nNewMatchMany\n"
      << pre << name_ 
      << std::setw(15) << numHardestRefFailMatch/double(nTotal)
      << std::setw(15) << numHardestNewFailMatch/double(nTotal)
      << std::setw(15) << numHardestRefMatchedNotHardestNew/double(nTotal)
      << std::setw(15) << numRefMatchedMultipleTimes/double(nTotal)
      << std::setw(15) << numNewMatchedMultipleTimes/double(nTotal) << '\n';
  os << oss.str();
}

} // namespace mc