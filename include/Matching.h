#include <iostream>
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/common/include/CommonModules.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/HOTVR.hh"
#include "fastjet/contrib/HOTVRinfo.hh"

/** This class has two important purposes.
* 1. convert selected genparticles into pseudojets (those to be clustered)
* 2. matching a parton jet containing a top to the correponding particle jet (run_matching)
*/

class Matching{
 private:
  std::vector<fastjet::PseudoJet> _stableParticles;
  std::vector<fastjet::PseudoJet> _partons_to_cluster;
  std::vector<fastjet::PseudoJet> _Wb_partons_to_cluster;

  std::vector<TopJet> _matched_jets;
  std::vector<TopJet> _matched_parton_jets;
  std::vector<std::pair<TopJet, TopJet>> _matched_pairs;

  std::vector<fastjet::PseudoJet> _matched_pseudojets;
  std::vector<std::pair<fastjet::PseudoJet, std::vector<fastjet::PseudoJet>>> _matched_pseudojets_and_constituents;

  fastjet::PseudoJet convert_particle(GenParticle* genparticle);
  bool IsNeutrino(GenParticle* p);
  bool IsTop(GenParticle* p);
  bool IsW(GenParticle* p);
  bool IsZ(GenParticle* p);
  bool IsH(GenParticle* p);
  bool Isb(GenParticle* p);
  bool IsParton(GenParticle* p);
  bool IsStableParticle(GenParticle* p);
  bool BeforeTopDecay(GenParticle* p, std::vector<GenParticle>* genparticles);
  bool BeforeWZDecay(GenParticle* p, std::vector<GenParticle>* genparticles);
  bool BeforeHDecay(GenParticle* p, std::vector<GenParticle>* genparticles);
  bool FinalStateParton(GenParticle* p, std::vector<GenParticle>* genparticles);
  GenParticle* GetDaughter(GenParticle* p, std::vector<GenParticle>* genparticles, int n);

  bool IsHadronicTop(GenParticle* p,  std::vector<GenParticle>* genparticles);
  bool IsHadronicWZH(GenParticle* p,  std::vector<GenParticle>* genparticles);
  void UpdateSkipList(GenParticle* top, std::vector<GenParticle>* genparticles, std::vector<bool>& skip);

  bool IsMatched(TopJet jet, double matching_distance, TopJet refjet);
  bool IsMatched(fastjet::PseudoJet jet, double matching_distance, fastjet::PseudoJet refjet);

 public:
  Matching();
  void run_matching(std::vector<TopJet> particle_jets, std::vector<TopJet> parton_jets);
  void run_matching(std::vector<fastjet::PseudoJet> jets, std::vector<fastjet::PseudoJet> denominator_jets);
  void run_matching(std::vector<fastjet::PseudoJet> jets, std::vector<std::vector<fastjet::PseudoJet>> constituents, std::vector<fastjet::PseudoJet> denominator_jets);
  int  PrimaryEventPDGID(std::vector<GenParticle> *genparticles);
  bool IsFullyHadronicEvent(std::vector<GenParticle> *genparticles);
  bool IsGravitonEvent(std::vector<GenParticle> *genparticles);
  void set_stable_particles(std::vector<GenParticle>* genparticles);
  void set_partons(std::vector<GenParticle>* genparticles, int target=6);
  void set_Wb_partons_from_top(std::vector<GenParticle>* genparticles);

  uint i_closest_particlejet(TopJet refjet, std::vector<TopJet> jets);
  TopJet get_closest_jet(TopJet jet, std::vector<TopJet> jets);

  void Reset();

  std::vector<fastjet::PseudoJet> get_stableParticles() { return _stableParticles; };
  std::vector<fastjet::PseudoJet> get_partons() { return _partons_to_cluster; };
  std::vector<fastjet::PseudoJet> get_Wb_partons() { return _Wb_partons_to_cluster; };

  std::vector<TopJet> get_matched_jets() { return _matched_jets; };
  std::vector<TopJet> get_matched_parton_jets() { return _matched_parton_jets; };
  std::vector<std::pair<TopJet, TopJet>> get_matched_pairs() { return _matched_pairs; };

  std::vector<fastjet::PseudoJet> get_matched_pseudojets() { return _matched_pseudojets; };
  std::vector<std::pair<fastjet::PseudoJet, std::vector<fastjet::PseudoJet>>> get_matched_jets_and_constituents() { return _matched_pseudojets_and_constituents; };
};
