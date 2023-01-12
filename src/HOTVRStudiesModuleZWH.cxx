#include <iostream>
#include <memory>
// include general UHH2 classes
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/PrintingModules.h"

// include everything from the HOTVR Studies
#include "UHH2/HOTVRStudies/include/Matching.h"
#include "UHH2/HOTVRStudies/include/Clustering.h"
#include "UHH2/HOTVRStudies/include/TopTagger.h"
#include "UHH2/HOTVRStudies/include/TopTagPerformanceHists.h"
#include "UHH2/HOTVRStudies/include/HOTVRJetsHists.h"
#include "UHH2/HOTVRStudies/include/VRJetsHists.h"
#include "UHH2/HOTVRStudies/include/SoftClusterHists.h"
#include "UHH2/HOTVRStudies/include/Makefiles.h"
#include "UHH2/HOTVRStudies/include/EfficiencyHists.h"
// include from HOTVR for the performance hists
#include "UHH2/HOTVR/include/HOTVRIds.h"
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

class HOTVRStudiesModuleZWH: public AnalysisModule {
public:

  explicit HOTVRStudiesModuleZWH(Context & ctx);
  virtual bool process(Event & event) override;

private:
//initialize hist classes
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_200;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_400;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_600;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_800;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_1000;

  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub2;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub2_200;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub2_400;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub2_600;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub2_800;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub2_1000;

  std::unique_ptr<HOTVRJetsHists> hist_matched_jets;
  std::unique_ptr<HOTVRJetsHists> hist_matched_jets_200;
  std::unique_ptr<HOTVRJetsHists> hist_matched_jets_400;
  std::unique_ptr<HOTVRJetsHists> hist_matched_jets_600;
  std::unique_ptr<HOTVRJetsHists> hist_matched_jets_800;
  std::unique_ptr<HOTVRJetsHists> hist_matched_jets_1000;

  std::unique_ptr<HOTVRJetsHists> hist_tagged_jets;
  std::unique_ptr<HOTVRJetsHists> hist_tagged_jets_200;
  std::unique_ptr<HOTVRJetsHists> hist_tagged_jets_400;
  std::unique_ptr<HOTVRJetsHists> hist_tagged_jets_600;
  std::unique_ptr<HOTVRJetsHists> hist_tagged_jets_800;
  std::unique_ptr<HOTVRJetsHists> hist_tagged_jets_1000;

  std::unique_ptr<HOTVRJetsHists> hist_njets_matched;

// initialize event handle
  Event::Handle<vector<TopJet>> h_HOTVR_jets;
  Event::Handle<vector<TopJet>> h_parton_jets;

  Event::Handle<vector<TopJet>> h_matched_jets;
  Event::Handle<vector<TopJet>> h_matched_parton_jets;
  Event::Handle<vector<pair<TopJet, TopJet>>> h_matched_pairs; //contains pair of matched jets

  int n_points = 100;

// initialize classes
  string m_clustering;
  Matching* matching;
  Clustering* clustering;
  TopTagger* toptagger;

// initialize vectors of jets
// ... containing pseudojets
  vector<fastjet::PseudoJet> pseudojets;
  vector<fastjet::PseudoJet> hotvr_jets;
  vector<vector<fastjet::PseudoJet>> hotvr_jets_constituents;
  vector<fastjet::PseudoJet> parton_pseudojets;
  vector<fastjet::PseudoJet> parton_jets;

  vector<Jet> _rejected_subjets; // rejected subjets with ptsub
  vector<TopJet> _rejected_cluster; // rejected cluster (jets without subjets)
  vector<TopJet> _soft_cluster; // soft cluster rejected via softdrop / massjump condition

  string dataset_version;

  bool isTTbar, isW, isZ, isH;
  bool is_mc, is_qcd;

  int nevent=0;
  bool debug = false;
};

/*
 ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
 ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

HOTVRStudiesModuleZWH::HOTVRStudiesModuleZWH(Context & ctx){
  cout << "Starting HOTVRStudiesModuleZWH!" << endl;
  if(debug){cout << "HOTVRStudiesModuleZWH: Debugging mode :) " << '\n';}
// get info from xml
  m_clustering = ctx.get("Clustering");
// check for the dataset version (ttbar or QCD)
  dataset_version = ctx.get("dataset_version");
  isTTbar = dataset_version.find("ttbar") == 0;
  isW = dataset_version.find("WW") == 0;
  isH = dataset_version.find("HH") == 0;
  isZ = dataset_version.find("ZZ") == 0;

  is_mc = ctx.get("dataset_type") == "MC";
  is_qcd = (dataset_version.find("QCD") == 0);

  h_HOTVR_jets = ctx.get_handle<vector<TopJet>>("HOTVR_jets");
  h_parton_jets = ctx.get_handle<vector<TopJet>>("parton_jets");

//declare event handle
  h_matched_jets = ctx.get_handle<vector<TopJet>>("matched_jets");
  h_matched_parton_jets = ctx.get_handle<vector<TopJet>>("matched_parton_jets");
  h_matched_pairs = ctx.get_handle<vector<pair<TopJet, TopJet>>>("matched_pairs");

// Set up Hists classes:
  hist_hotvr_jets.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets", is_qcd));
  hist_hotvr_jets_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_200", is_qcd));
  hist_hotvr_jets_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_400", is_qcd));
  hist_hotvr_jets_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_600", is_qcd));
  hist_hotvr_jets_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_800", is_qcd));
  hist_hotvr_jets_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_1000", is_qcd));

  hist_hotvr_jets_Nsub2.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub2", is_qcd));
  hist_hotvr_jets_Nsub2_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub2_200", is_qcd));
  hist_hotvr_jets_Nsub2_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub2_400", is_qcd));
  hist_hotvr_jets_Nsub2_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub2_600", is_qcd));
  hist_hotvr_jets_Nsub2_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub2_800", is_qcd));
  hist_hotvr_jets_Nsub2_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub2_1000", is_qcd));

  hist_matched_jets.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_matched_jets", is_qcd));
  hist_matched_jets_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_matched_jets_200", is_qcd));
  hist_matched_jets_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_matched_jets_400", is_qcd));
  hist_matched_jets_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_matched_jets_600", is_qcd));
  hist_matched_jets_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_matched_jets_800", is_qcd));
  hist_matched_jets_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_matched_jets_1000", is_qcd));

  hist_tagged_jets.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_tagged_jets", is_qcd));
  hist_tagged_jets_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_tagged_jets_200", is_qcd));
  hist_tagged_jets_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_tagged_jets_400", is_qcd));
  hist_tagged_jets_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_tagged_jets_600", is_qcd));
  hist_tagged_jets_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_tagged_jets_800", is_qcd));
  hist_tagged_jets_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_tagged_jets_1000", is_qcd));

  hist_njets_matched.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_njets_matched", is_qcd));

}
/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/
bool HOTVRStudiesModuleZWH::process(Event & event) {
  if(debug){std::cout << "Begin process..." << '\n';}

  //get the topjets from the Clustering Module
  vector<TopJet> _top_hotvr_jets = event.get(h_HOTVR_jets);
  vector<TopJet> _top_parton_jets = event.get(h_parton_jets);

//fill hists with hotvr jets and corresponding rejected subjets
  hist_hotvr_jets->fill_n_jets(event, _top_hotvr_jets);
  for(uint j=0; j<_top_hotvr_jets.size(); ++j){ // loop over hotvr jets
    //std::cout << "Jet pt "<<jet.pt() << '\n';
    TopJet jet = _top_hotvr_jets[j];
    hist_hotvr_jets->fill_topjet(event, jet);
    if(jet.pt()>200 &&jet.pt()<400)  hist_hotvr_jets_200->fill_topjet(event, jet);
    if(jet.pt()>400 &&jet.pt()<600)  hist_hotvr_jets_400->fill_topjet(event, jet);
    if(jet.pt()>600 &&jet.pt()<800)  hist_hotvr_jets_600->fill_topjet(event, jet);
    if(jet.pt()>800 &&jet.pt()<1000)  hist_hotvr_jets_800->fill_topjet(event, jet);
    if(jet.pt()>1000 &&jet.pt()<1200)  hist_hotvr_jets_1000->fill_topjet(event, jet);

// fill hists after Nsub2 cut
if (jet.subjets().size()>1) {
    hist_hotvr_jets_Nsub2->fill_topjet(event, jet);
    if(jet.pt()>200 &&jet.pt()<400)  hist_hotvr_jets_Nsub2_200->fill_topjet(event, jet);
    if(jet.pt()>400 &&jet.pt()<600)  hist_hotvr_jets_Nsub2_400->fill_topjet(event, jet);
    if(jet.pt()>600 &&jet.pt()<800)  hist_hotvr_jets_Nsub2_600->fill_topjet(event, jet);
    if(jet.pt()>800 &&jet.pt()<1000)  hist_hotvr_jets_Nsub2_800->fill_topjet(event, jet);
    if(jet.pt()>1000 &&jet.pt()<1200)  hist_hotvr_jets_Nsub2_1000->fill_topjet(event, jet);
  }
} // end loop over hotvr jets

// ------MATCHING--------
matching = new Matching();
//run_matching: loop over the parton jets and match them to the hotvr jets
  vector<TopJet> matched_jets;
  vector<TopJet> matched_parton_jets;
  matching->run_matching(_top_hotvr_jets, _top_parton_jets);
  matched_jets = matching->get_matched_jets();
  matched_parton_jets = matching->get_matched_parton_jets();
  vector<pair<TopJet, TopJet>> matched_pair = matching->get_matched_pairs();

  hist_njets_matched->fill_n_jets(event, matched_jets);
//TopTagger for matched hotvr jets
  toptagger = new TopTagger();
  vector<TopJet> matched_jets_tagged;
  vector<pair<TopJet, TopJet>> matched_pair_tagged;
  for(uint j=0; j<matched_pair.size(); ++j){ // loop over matched jets
    TopJet parton_jet=matched_pair[j].second;
    TopJet matched_jet=matched_pair[j].first;

    //fill hists with matched jets
    hist_matched_jets->fill_topjet(event, matched_jet);
    if(parton_jet.pt()>200 && parton_jet.pt()<400)  hist_matched_jets_200->fill_topjet(event, matched_jet);
    if(parton_jet.pt()>400 && parton_jet.pt()<600)  hist_matched_jets_400->fill_topjet(event, matched_jet);
    if(parton_jet.pt()>600 && parton_jet.pt()<800)  hist_matched_jets_600->fill_topjet(event, matched_jet);
    if(parton_jet.pt()>800 && parton_jet.pt()<1000)  hist_matched_jets_800->fill_topjet(event, matched_jet);
    if(parton_jet.pt()>1000 && parton_jet.pt()<1200)  hist_matched_jets_1000->fill_topjet(event, matched_jet);

// --------apply Tagger-----------------
    string tagger;
    if (isW) tagger = "W";
    if (isZ) tagger = "Z";
    if (isH) tagger = "H";

    if(toptagger->Is_tagged(tagger, matched_jets[j])){
      matched_jets_tagged.push_back(matched_jet);
      matched_pair_tagged.push_back(matched_pair[j]);
      // fill hists with tagged jets
      hist_tagged_jets->fill_topjet(event, matched_jet);

      if(parton_jet.pt()>200 && parton_jet.pt()<400){
        hist_tagged_jets_200->fill_topjet(event, matched_jet);
        }

      if(parton_jet.pt()>400 && parton_jet.pt()<600){
        hist_tagged_jets_400->fill_topjet(event, matched_jet);
      }

      if(parton_jet.pt()>600 && parton_jet.pt()<800){
        hist_tagged_jets_600->fill_topjet(event, matched_jet);
      }

      if(parton_jet.pt()>800 && parton_jet.pt()<1000){
        hist_tagged_jets_800->fill_topjet(event, matched_jet);
      }

      if(parton_jet.pt()>1000 && parton_jet.pt()<1200){
        hist_tagged_jets_1000->fill_topjet(event, matched_jet);
      }

  } // end if top tag
  } // end loop over matched jets

// decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the HOTVRStudiesModuleZWH is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(HOTVRStudiesModuleZWH)
