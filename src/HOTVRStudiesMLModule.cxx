#include <iostream>
#include <memory>
// include general UHH2 classes
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
// include everything from the HOTVR Studies
#include "UHH2/HOTVRStudies/include/Matching.h"
#include "UHH2/HOTVRStudies/include/Clustering.h"
#include "UHH2/HOTVRStudies/include/TopTagger.h"
#include "UHH2/HOTVRStudies/include/HOTVRJetsHists.h"
#include "UHH2/HOTVRStudies/include/VRJetsHists.h"
#include "UHH2/HOTVRStudies/include/Makefiles.h"
// include from fastjet
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/HOTVRinfo.hh"
#include "fastjet/contrib/HOTVR.hh"

using namespace std;
using namespace uhh2;
using namespace fastjet;
using namespace contrib;

/** AnalysisModule for HOTVR Studies.
*
* This is the central class which calls other AnalysisModules, Hists or Selection classes.
* This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
*/

class HOTVRStudiesMLModule: public AnalysisModule {
public:

  explicit HOTVRStudiesMLModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
//initialize hist classes
  std::unique_ptr<VRJetsHists> hist_jets;
  std::unique_ptr<VRJetsHists> hist_jets_200;
  std::unique_ptr<VRJetsHists> hist_jets_350;
  std::unique_ptr<VRJetsHists> hist_jets_500;

// initialize event handle
  Event::Handle<vector<vector<double>>> h_jet_info_px;
  Event::Handle<vector<vector<double>>> h_jet_info_py;
  Event::Handle<vector<vector<double>>> h_jet_info_pz;
  Event::Handle<vector<vector<double>>> h_jet_info_m;

// initialize classes
  string m_clustering;
  Matching* matching;
  Clustering* clustering;

// initialize vectors of jets
// ... containing pseudojets
  vector<fastjet::PseudoJet> pseudojets;
  vector<fastjet::PseudoJet> jets;
  vector<vector<fastjet::PseudoJet>> jets_constituents;
  vector<fastjet::PseudoJet> parton_pseudojets;
  vector<fastjet::PseudoJet> parton_jets;

  string dataset_version;

  bool isTTbar;
  bool is_mc, is_qcd;

  int nevent=0;
};

/*
 ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
 ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

HOTVRStudiesMLModule::HOTVRStudiesMLModule(Context & ctx){
  cout << "Starting HOTVRStudiesMLModule!" << endl;
// get info from xml
  m_clustering = ctx.get("Clustering");
// check for the dataset version (ttbar or QCD)
  dataset_version = ctx.get("dataset_version");
  isTTbar = dataset_version.find("ttbar") == 0;
  is_mc = ctx.get("dataset_type") == "MC";
  is_qcd = (dataset_version.find("QCD") == 0);

  ctx.undeclare_all_event_output();

//declare event handle
  //for ML studies
  h_jet_info_px = ctx.declare_event_output<vector<vector<double>>>("jet_info_px");
  h_jet_info_py = ctx.declare_event_output<vector<vector<double>>>("jet_info_py");
  h_jet_info_pz = ctx.declare_event_output<vector<vector<double>>>("jet_info_pz");
  h_jet_info_m = ctx.declare_event_output<vector<vector<double>>>("jet_info_m");

// Set up Hists classes:
  hist_jets.reset(new VRJetsHists(ctx, "VRJetsHists_jets", is_qcd));
  hist_jets_200.reset(new VRJetsHists(ctx, "VRJetsHists_jets_200", is_qcd));
  hist_jets_350.reset(new VRJetsHists(ctx, "VRJetsHists_jets_350", is_qcd));
  hist_jets_500.reset(new VRJetsHists(ctx, "VRJetsHists_jets_500", is_qcd));

}
/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/
bool HOTVRStudiesMLModule::process(Event & event) {
// here define genparticles and read from event
  auto genparticles = event.genparticles;

//MATCHING -> class to filter the stable particles from the genparticles and create a pseudojet from them
  matching = new Matching();
  matching->set_stable_particles(genparticles);
  pseudojets = matching->get_stableParticles();
  matching->set_partons(genparticles);
  parton_pseudojets = matching->get_partons();
//CLUSTERING -> clusters jets depending on the algorithm that is chosen
  clustering = new Clustering(m_clustering);
  clustering->cluster_jets(pseudojets); // cluster the pseudojets, possible modes defined in hotvr.config: "HOTVR, HOTVR_SD, VR"
  jets = clustering->get_clustered_jets();
  jets_constituents = clustering->get_clustered_jet_constituents();
  // cluster and get the clustered parton jets (AK10 jets)
  clustering->cluster_parton_jets(parton_pseudojets, isTTbar);
  parton_jets = clustering->get_parton_pseudojets();
  matching->run_matching(jets, jets_constituents, parton_jets); // match the clustered jets to parton level jets
  vector<pair<PseudoJet, vector<PseudoJet>>> matched_jets_and_constituents = matching->get_matched_jets_and_constituents();
  vector<vector<double>> jets_constituents_px;
  vector<vector<double>> jets_constituents_py;
  vector<vector<double>> jets_constituents_pz;
  vector<vector<double>> jets_constituents_m;

  vector<double> constituents_px;
  vector<double> constituents_py;
  vector<double> constituents_pz;
  vector<double> constituents_m;
  for (size_t j = 0; j < matched_jets_and_constituents.size(); j++) { // loop over jets to fill hists
    auto jet = matched_jets_and_constituents[j].first;
    auto jet_constituents = matched_jets_and_constituents[j].second;
    hist_jets->fill_pseudojet(jet); // fill basic hists (no specific HOTVR info like subjets etc.)
    if(jet.pt()>200 && jet.pt()<350)  hist_jets_200->fill_pseudojet(jet);
    if(jet.pt()>350 && jet.pt()<500)  hist_jets_350->fill_pseudojet(jet);
    if(jet.pt()>500)  hist_jets_500->fill_pseudojet(jet);
  //  hist_jets->fill_pseudojet_constituents(jet_constituents);
    for (size_t n = 0; n < jet_constituents.size(); n++) { // loop over constituents
        auto constituent = jet_constituents[n];
        constituents_px.push_back(constituent.px()); // save double in constituents vector
        constituents_py.push_back(constituent.py());
        constituents_pz.push_back(constituent.pz());
        constituents_m.push_back(constituent.m());
      }
    jets_constituents_px.push_back(constituents_px); // add constituents to jet list
    jets_constituents_py.push_back(constituents_py);
    jets_constituents_pz.push_back(constituents_pz);
    jets_constituents_m.push_back(constituents_m);
  constituents_px.clear();
  constituents_py.clear();
  constituents_pz.clear();
  constituents_m.clear();
  }
//set the event handle
  // event.set(h_jet_info_px, jets_constituents_px);
  // event.set(h_jet_info_py, jets_constituents_py);
  // event.set(h_jet_info_pz, jets_constituents_pz);
  // event.set(h_jet_info_m, jets_constituents_m);
  // delete clustering infos
  clustering->Reset();
// decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the HOTVRStudiesMLModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(HOTVRStudiesMLModule)
