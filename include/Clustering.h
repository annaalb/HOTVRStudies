#include <iostream>
#include <memory>

#include "UHH2/common/include/CommonModules.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
//#include "fastjet/contrib/IteratedSoftDrop.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/HOTVRinfo.hh"
#include "fastjet/contrib/HOTVR.hh"
#include "fastjet/contrib/VariableRPlugin.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

#include <math.h>
#include <random>


/** Class for the clustering. The class uses fastjet plugins to cluster jets from input pseudojets.
* The clustered jets can be obtained via getter methods.
*/
class Clustering{
private:
  fastjet::ClusterSequence* _clust_seq;
  bool _settings_not_shown=true;

  double _mu, _theta, _max_r, _min_r, _rho, _subjetptmin, _jetptmin, _alpha;
  double _z_cut, _beta, _pt_threshold;
  string _clustering_algorithmus;
  double _a, _b, _c;

  std::vector<fastjet::PseudoJet> _hotvr_jets;
  std::vector<fastjet::PseudoJet> _vr_jets;
  std::vector<fastjet::PseudoJet> _vr_jets_SD;
  std::vector<fastjet::PseudoJet> _vr_jets_ISD;  // Variable R plus iterative SD
  std::vector<fastjet::PseudoJet> _parton_fatjets;
  std::vector<fastjet::PseudoJet> _W_parton_fatjets;
  std::vector<fastjet::PseudoJet> _b_parton_fatjets;

// save jet constituents
  std::vector<std::vector<fastjet::PseudoJet>> _vr_jet_constituents;
  std::vector<std::vector<fastjet::PseudoJet>> _hotvr_jet_constituents;
  std::vector<std::vector<fastjet::PseudoJet>> _rejected_cluster_constituents;
  std::vector<std::vector<fastjet::PseudoJet>> _soft_cluster_constituents;
  std::vector<std::vector<fastjet::PseudoJet>> _rejected_subjets_constituents;

  std::vector<std::vector<fastjet::PseudoJet>> _jet0_subjets_constituents;
  std::vector<std::vector<fastjet::PseudoJet>> _jet1_subjets_constituents;
  std::vector<std::vector<fastjet::PseudoJet>> _jet2_subjets_constituents;
  std::vector<std::vector<fastjet::PseudoJet>> _jet3_subjets_constituents;

// jets rejected via grooming in the clustering process
  std::vector<fastjet::PseudoJet> _rejected_cluster;
  std::vector<fastjet::PseudoJet> _soft_cluster;
  std::vector<fastjet::PseudoJet> _rejected_subjets;

  std::vector<TopJet> _top_hotvr_jets;
  std::vector<TopJet> _top_vr_jets;
  std::vector<TopJet> _top_parton_jets;
  std::vector<TopJet> _W_parton_jets;
  std::vector<TopJet> _b_parton_jets;

  std::vector<TopJet> _top_rejected_cluster;
  std::vector<TopJet> _top_soft_cluster;
  std::vector<Jet> _top_rejected_subjets;

  TopJet convert_jet(fastjet::PseudoJet pseudojet, double tau1, double tau2, double tau3);
  TopJet convert_jet(fastjet::PseudoJet pseudojet);
  std::vector<Jet> convert_subjets(std::vector<fastjet::PseudoJet> subjets);
  double get_value(std::string word);
  void show_settings();
  std::vector<fastjet::PseudoJet> add_ghosts(std::vector<fastjet::PseudoJet> gen_in);
  std::vector<std::vector<fastjet::PseudoJet>> save_constituents(std::vector<fastjet::PseudoJet> jet_in);

public:
// Constructor
  Clustering();
  Clustering(std::string clustering);
  void Reset();
// clustering method called in the HOTVRStudiesModule
  void cluster_jets(std::vector<fastjet::PseudoJet>, int nevent=0);
// specific clustering methods for HOTVR; HOTVR with SD; VR SD; AK10
  void cluster_HOTVR_jets(std::vector<fastjet::PseudoJet>);
  void cluster_HOTVR_SD_jets(std::vector<fastjet::PseudoJet>, int nevent=0);
  void cluster_VR_SD_jets(std::vector<fastjet::PseudoJet>);
  void cluster_VR_ISD_jets(std::vector<fastjet::PseudoJet>);
  void cluster_parton_jets(std::vector<fastjet::PseudoJet>, bool);
  void cluster_W_parton_jets(std::vector<fastjet::PseudoJet>);
  void cluster_b_parton_jets(std::vector<fastjet::PseudoJet>);

// Getter for pseudojets
std::vector<fastjet::PseudoJet> get_clustered_jets(); // depending on clustering mode in config file
std::vector<std::vector<fastjet::PseudoJet>> get_clustered_jet_constituents();
//HOTVR clustering
  std::vector<fastjet::PseudoJet> get_hotvr_jets(){return _hotvr_jets;};
//Variable R plus grooming
  std::vector<fastjet::PseudoJet> get_vr_jets(){return _vr_jets;};
  std::vector<fastjet::PseudoJet> get_vr_jets_SD(){return _vr_jets_SD;};
  std::vector<fastjet::PseudoJet> get_vr_jets_ISD(){return _vr_jets_ISD;};
// clustering on parton level
  std::vector<fastjet::PseudoJet> get_parton_jets(){return _parton_fatjets;};
  std::vector<fastjet::PseudoJet> get_W_parton_jets(){return _W_parton_fatjets;};
  std::vector<fastjet::PseudoJet> get_b_parton_jets(){return _b_parton_fatjets;};


// Getter for jet constituents
std::vector<std::vector<fastjet::PseudoJet>> get_vr_jet_constituents(){return _vr_jet_constituents;};
std::vector<std::vector<fastjet::PseudoJet>> get_hotvr_jet_constituents(){return _hotvr_jet_constituents;};

// Getter for TopJets
  std::vector<TopJet> get_top_hotvr_jets(){return _top_hotvr_jets;};
  std::vector<TopJet> get_top_vr_jets(){return _top_vr_jets;};
  std::vector<TopJet> get_top_parton_jets(){return _top_parton_jets;};
  std::vector<TopJet> get_top_W_parton_jets(){return _W_parton_jets;};
  std::vector<TopJet> get_top_b_parton_jets(){return _b_parton_jets;};

// Getter for rejected jets or subjets or soft jets
  std::vector<fastjet::PseudoJet> get_rejected_cluster(){return _rejected_cluster;};
  std::vector<fastjet::PseudoJet> get_soft_cluster(){return _soft_cluster;};
  std::vector<fastjet::PseudoJet> get_rejected_subjets(){return _rejected_subjets;};

  std::vector<std::vector<fastjet::PseudoJet>> get_rejected_cluster_constituents(){return _rejected_cluster_constituents;};
  std::vector<std::vector<fastjet::PseudoJet>> get_soft_cluster_constituents(){return _soft_cluster_constituents;};
  std::vector<std::vector<fastjet::PseudoJet>> get_rejected_subjets_constituents(){return _rejected_subjets_constituents;};

  std::vector<std::vector<fastjet::PseudoJet>> get_jet0_subjets_constituents(){return _jet0_subjets_constituents;};
  std::vector<std::vector<fastjet::PseudoJet>> get_jet1_subjets_constituents(){return _jet1_subjets_constituents;};
  std::vector<std::vector<fastjet::PseudoJet>> get_jet2_subjets_constituents(){return _jet2_subjets_constituents;};
  std::vector<std::vector<fastjet::PseudoJet>> get_jet3_subjets_constituents(){return _jet3_subjets_constituents;};

};
