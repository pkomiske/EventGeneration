/** 
 *  Patrick Komiske, MIT, Copyright (C) 2020.
 */

// Events generation with Pythia and jet finding with FastJet
#include "EventGenerator.hh"

int main(int argc, char** argv) {

  // user comments
  std::ostringstream oss;
  oss << "# user_values: [SoftDrop(zcut=0.1, beta=0).mass, SoftDrop(zcut=0.1, beta=1).mass,\n"
      << "#               SafeDrop(zcut=0.1, f=0.0).mass, SafeDrop(zcut=0.1, f=0.25).mass,\n"
      << "#               SafeDrop(zcut=0.1, f=0.5).mass, SafeDrop(zcut=0.1, f=0.75).mass,\n"
      << "#               SafeDrop(zcut=0.1, f=1.0).mass, SafeDrop(zcut=0.2, f=0.0).mass,\n"
      << "#               SafeDrop(zcut=0.2, f=0.25).mass, SafeDrop(zcut=0.2, f=0.5).mass,\n"
      << "#               SafeDrop(zcut=0.2, f=0.75).mass, SafeDrop(zcut=0.2, f=1.0).mass,\n"
      << "#               SafeDrop(zcut=0.4, f=0.0).mass, SafeDrop(zcut=0.4, f=0.25).mass,\n"
      << "#               SafeDrop(zcut=0.4, f=0.5).mass, SafeDrop(zcut=0.4, f=0.75).mass,\n"
      << "#               SafeDrop(zcut=0.4, f=1.0).mass,\n"
      << "#               IVGCR(Area, 0.1).(mass, emd), IVGCR(Area, 0.2).(mass, emd),\n"
      << "#               IVGCR(Area, 0.3).(mass, emd), IVGCR(Area, 0.4).(mass, emd),\n"
      << "#               IVGCR(EMD, 0.1).(mass, emd), IVGCR(EMD, 0.2).(mass, emd),\n"
      << "#               IVGCR(EMD, 0.2).(mass, emd), IVGCR(EMD, 0.3).(mass, emd)]\n";

  mc::EventGenerator evgen(argc, argv, oss.str());
  if (evgen.exit) return 1;

  uint i(0);
  double analysis_time(0);
  while(i < evgen.ntot) {

    // this advances pythia until we have a valid event with some jets
    // returns how many sets of jets we have to consider
    uint nunits(evgen.next());

    // iterate over units
    for (uint j = 0; j < nunits; j++) {

      // check that unit is fully_matched (only increment i on fully_matched units anyway)
      bool fully_matched(evgen.unit_fully_matched(j));
      if (!fully_matched && evgen.fully_matched_only) continue;

      // get specific jets 
      std::pair<fastjet::PseudoJet,bool> hadron_jet_pair(evgen.get_hadron_jet(j)), parton_jet_pair(evgen.get_parton_jet(j));

      // analyze hadron jet
      if (hadron_jet_pair.second) {
        std::chrono::steady_clock::time_point start(std::chrono::steady_clock::now());
        const fastjet::PseudoJet & hadron_jet(hadron_jet_pair.first);

        // determine hadron jet features
        std::vector<fastjet::PseudoJet> hadron_consts(hadron_jet.constituents());
        std::vector<double> hadron_features{double(hadron_consts.size()), 0};

        // count number of pions in the jet (just as an example)
        // (the user index indexes into the pid vector)
        for (const fastjet::PseudoJet & pj : hadron_consts) {
          if (abs(evgen.hadron_pids[pj.user_index()]) == 211)
            hadron_features.back() += 1;
        }
        
        // set features to be output
        evgen.set_hadron_jet_features(hadron_features);

        analysis_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      }

      // analyze parton jet
      if (parton_jet_pair.second) {
        std::chrono::steady_clock::time_point start(std::chrono::steady_clock::now());
        const fastjet::PseudoJet & parton_jet(parton_jet_pair.first);

        // determine parton jet features
        std::vector<fastjet::PseudoJet> parton_consts(parton_jet.constituents());
        std::vector<double> parton_features{double(parton_consts.size()), 0};

        // count number of pions in the jet (just as an example)
        // (the user index indexes into the pid vector)
        for (const fastjet::PseudoJet & pj : parton_consts) {
          if (abs(evgen.parton_pids[pj.user_index()]) == 211)
            parton_features.back() += 1;
        }
        
        // set features to be output
        evgen.set_parton_jet_features(parton_features);

        analysis_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      } 

      // output unit
      evgen.write_unit(j);

      // update i
      if (fully_matched) i++;

      // print update
      if (i % evgen.print_every == 0) evgen.print(i);

      // check we have enough units already
      if (i == evgen.ntot) break;
    }
  }

  // summarize
  evgen.summarize();

  // output time
  std::cout << std::setprecision(12) << "Total analysis time: " << analysis_time << '\n';

  return 0;
}
