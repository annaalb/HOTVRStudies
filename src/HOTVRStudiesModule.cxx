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

#include <UHH2/HOTVRStudies/include/ModuleBASE.h>

using namespace std;
using namespace uhh2;
using namespace fastjet;
using namespace contrib;

/** AnalysisModule for HOTVR Studies.
*
* This is the central class which calls other AnalysisModules, Hists or Selection classes.
* This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
*/

class HOTVRStudiesModule: public ModuleBASE {
public:

  explicit HOTVRStudiesModule(Context & ctx);
  virtual bool process(Event & event) override;
  void book_histograms(uhh2::Context&, string, vector<string>, bool, string);
  void fill_histograms(uhh2::Event&, TopJet, TopJet, string, string);


private:
  //initialize hist classes
  vector<string> pTs;

  std::unique_ptr<HOTVRJetsHists> hist_njets_matched;

  // initialize event handle
  Event::Handle<vector<TopJet>> h_HOTVR_jets;
  Event::Handle<vector<TopJet>> h_parton_jets;

  Event::Handle<vector<TopJet>> h_matched_jets;
  Event::Handle<vector<TopJet>> h_matched_parton_jets;
  Event::Handle<vector<pair<TopJet, TopJet>>> h_matched_pairs; //contains pair of matched jets
  // denominator for the ROC curves first: ttbar, second:qcd

  // for performance study
  std::unique_ptr<AnalysisModule> cl_topjet_hotvr;
  std::unique_ptr<Selection> sel_ntop_hotvr;
  std::unique_ptr<Selection> sel_toptag_hotvr, sel_toptag_softdrop;

  // --- Scans --- //
  // std::unique_ptr<Hists> hist_hotvr_pre, hist_softdrop_pre;
  // std::vector<std::unique_ptr<Selection>> sel_hotvr, sel_softdrop;
  // std::vector<std::unique_ptr<Hists>> hist_hotvr, hist_softdrop;
  int n_points = 100;
  //----EfficiencyHists----
  //  std::unique_ptr<EfficiencyHists> hist_eff_hotvr;

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
  vector<fastjet::PseudoJet> vr_jets;
  vector<vector<fastjet::PseudoJet>> vr_jets_constituents;
  vector<fastjet::PseudoJet> vr_jets_SD;
  vector<fastjet::PseudoJet> vr_jets_ISD;
  vector<fastjet::PseudoJet> parton_pseudojets;
  vector<fastjet::PseudoJet> parton_jets;

  vector<Jet> _rejected_subjets; // rejected subjets with ptsub
  vector<TopJet> _rejected_cluster; // rejected cluster (jets without subjets)
  vector<TopJet> _soft_cluster; // soft cluster rejected via softdrop / massjump condition

  string dataset_version;

  bool isTTbar;
  bool is_mc, is_qcd;

  int nevent=0;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  Event::Handle<TTbarGen> h_ttbargen;
  bool debug = false;
};

void HOTVRStudiesModule::book_histograms(uhh2::Context& ctx, string jet_coll, vector<string> pTs, bool is_qcd, string selection=""){
  for (auto pt:pTs){
    book_HFolder("hist_"+jet_coll+"_jets"+selection+"_"+pt, new HOTVRJetsHists(ctx,"HOTVRJetsHists_"+jet_coll+"_jets"+selection+"_"+pt, is_qcd));
  }
  book_HFolder("hist_"+jet_coll+"_jets"+selection, new HOTVRJetsHists(ctx,"HOTVRJetsHists_"+jet_coll+"_jets"+selection, is_qcd));
}

void HOTVRStudiesModule::fill_histograms(uhh2::Event& event,TopJet jet, TopJet pt_jet, string jet_coll, string selection=""){

  HFolder("hist_"+jet_coll+"_jets"+selection)->fill_topjet(event,jet); 

  if(pt_jet.pt()>200 &&pt_jet.pt()<400) HFolder("hist_"+jet_coll+"_jets"+selection+"_200")->fill_topjet(event,jet); 
  else if(pt_jet.pt()>400 &&pt_jet.pt()<600)HFolder("hist_"+jet_coll+"_jets"+selection+"_400")->fill_topjet(event,jet); 
  else if(pt_jet.pt()>600 &&pt_jet.pt()<800) HFolder("hist_"+jet_coll+"_jets"+selection+"_600")->fill_topjet(event,jet); 
  else if(pt_jet.pt()>800 &&pt_jet.pt()<1000) HFolder("hist_"+jet_coll+"_jets"+selection+"_800")->fill_topjet(event,jet); 
  else if(pt_jet.pt()>1000 &&pt_jet.pt()<1200) HFolder("hist_"+jet_coll+"_jets"+selection+"_1000")->fill_topjet(event,jet); 
}


/*
 ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
 ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

HOTVRStudiesModule::HOTVRStudiesModule(Context & ctx){
  cout << "Starting HOTVRStudiesModule!" << endl;
  if(debug){cout << "HOTVRStudiesModule: Debugging mode :) " << '\n';}
  // get info from xml
  m_clustering = ctx.get("Clustering");
  // check for the dataset version (ttbar or QCD)
  dataset_version = ctx.get("dataset_version");
  isTTbar = dataset_version.find("ttbar") == 0;
  is_mc = ctx.get("dataset_type") == "MC";
  is_qcd = (dataset_version.find("QCD") == 0);

  h_HOTVR_jets = ctx.get_handle<vector<TopJet>>("HOTVR_jets");
  h_parton_jets = ctx.get_handle<vector<TopJet>>("parton_jets");

  //declare event handle
  h_matched_jets = ctx.get_handle<vector<TopJet>>("matched_jets");
  h_matched_parton_jets = ctx.get_handle<vector<TopJet>>("matched_parton_jets");
  h_matched_pairs = ctx.get_handle<vector<pair<TopJet, TopJet>>>("matched_pairs");

  // Set up Hists classes:
  pTs = {"200","400", "600","800","1000"};
  book_histograms(ctx, "hotvr", pTs, is_qcd, "");

  //  book_histograms(ctx, "hotvr", pTs, is_qcd, "_mGr130");
  //  book_histograms(ctx, "hotvr", pTs, is_qcd, "_mSm130");

  book_histograms(ctx, "hotvr", pTs, is_qcd, "_Nsub3");

  // book_histograms(ctx, "hotvr", pTs, is_qcd, "_Nsubeq2");
  // book_histograms(ctx, "hotvr", pTs, is_qcd, "_Nsubeq3");
  // book_histograms(ctx, "hotvr", pTs, is_qcd, "_Nsub4");

  book_histograms(ctx, "hotvr", pTs, is_qcd, "_fpt");
  book_histograms(ctx, "hotvr", pTs, is_qcd, "_mass");

  book_histograms(ctx, "matched", pTs, is_qcd, "");

  // book_histograms(ctx, "matched", pTs, is_qcd, "_Nsub3");
  // book_histograms(ctx, "matched", pTs, is_qcd, "_fpt");

  book_histograms(ctx, "tagged", pTs, is_qcd, "");
  book_histograms(ctx, "tagged", pTs, is_qcd, "_tau32");

  hist_njets_matched.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_njets_matched", is_qcd));

  //for ROC curves
  double top_pt_min = 0.0;
  double top_eta_max = 100000;
  TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);
  TopJetId id_hotvr = AndId<TopJet>(id_topjet, Tau32Groomed(0.47));  // definition of hotvr id (HOTVR SD WP 0.47)

  // for (int i = 0; i < n_points; ++i) //loop over working points and create hists
  // {
  //   std::unique_ptr<Selection> sel;
  //   std::unique_ptr<Hists> hist;
  //   TopJetId id;
  //   double tau32 = 1 - (i * 1./n_points);
  //   // - HOTVR - //
  //   id = AndId<TopJet>(id_topjet, Tau32Groomed(tau32)); // define id on specific tau
  //   sel.reset(new NTopJetSelection(1, -1, id, h_matched_jets)); // selection on specific tau
  //   sel_hotvr.push_back(std::move(sel)); // fill vector with selection bool
  //
  //   hist.reset(new TopTagPerformanceHists(ctx, "HOTVR_Performance_" + to_string(i), is_qcd,  id, h_matched_pairs)); // reserve hist for sp. tau
  //   hist_hotvr.push_back(std::move(hist)); // fill vector with hists
  // }
  //
  // hist_hotvr_pre.reset(new TopTagPerformanceHists(ctx, "HOTVR_Pre", is_qcd, id_topjet, h_matched_pairs)); // reserve preselection hist
  //
  // //efficiency hists -> give hotvr id
  // hist_eff_hotvr.reset(new EfficiencyHists(ctx, "HOTVR_Eff", is_qcd, id_hotvr, h_matched_pairs));

}
/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/
bool HOTVRStudiesModule::process(Event & event) {
  if(debug){std::cout << "Begin process..." << '\n';}

  //get the topjets from the Clustering Module
  vector<TopJet> _top_hotvr_jets = event.get(h_HOTVR_jets);
  vector<TopJet> _top_parton_jets = event.get(h_parton_jets);

  //fill hists with hotvr jets and corresponding rejected subjets
  HFolder("hist_hotvr_jets")->fill_n_jets(event, _top_hotvr_jets);

  for(uint j=0; j<_top_hotvr_jets.size(); ++j){ // loop over hotvr jets
    //std::cout << "Jet pt "<<jet.pt() << '\n';
    TopJet jet = _top_hotvr_jets[j];
    fill_histograms(event, jet, jet, "hotvr",""); 
  
    // fill hists for m > 130 GeV and m < 130 GeV
    // if (jet.v4().M()>130) { fill_histograms(event, jet, jet, "hotvr","mGr130"); }
    // else if (jet.v4().M()<130) { fill_histograms(event, jet, jet, "hotvr","mSm130"); }

    // fill hists after Nsub3 cut
    if (jet.subjets().size()>2) { fill_histograms(event, jet, jet, "hotvr","Nsub3"); }

    // fill hists for jets with == 2 subjets
    // if (jet.subjets().size()==2) { fill_histograms(event, jet, jet, "hotvr","Nsubeq2"); } 

    // fill hists after Nsub3 cut
    //   if (jet.subjets().size()==3) { fill_histograms(event, jet, jet, "hotvr","Nsubeq3");}

    // fill hists after Nsub3 cut
    //     if (jet.subjets().size()>3) { fill_histograms(event, jet, jet, "hotvr","Nsub3");}

    // fill hists after fpt cut
    if (jet.subjets().size()>2 &&jet.hotvr_fpt1()<0.8) { fill_histograms(event, jet, jet, "hotvr","fpt"); }

    // fill hist after mass cut
    if (jet.subjets().size()>2 && 140<jet.v4().M() &&jet.v4().M()<220) { fill_histograms(event, jet, jet, "hotvr","mass"); }

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

  HFolder("hist_njets_matched")->fill_n_jets(event, matched_jets);

  //TopTagger for matched hotvr jets
  toptagger = new TopTagger();
  vector<TopJet> matched_jets_tagged;
  vector<pair<TopJet, TopJet>> matched_pair_tagged;
  for(uint j=0; j<matched_pair.size(); ++j){ // loop over matched jets
    TopJet parton_jet=matched_pair[j].second;
    TopJet matched_jet=matched_pair[j].first;

    //fill hists with matched jets
    fill_histograms(event, matched_jet, parton_jet, "matched",""); 

    // if(matched_jet.subjets().size()>2){ fill_histograms(event, matched_jet, parton_jet, "matched","_Nsub3");}
    // if(matched_jet.hotvr_fpt1()<0.8){ fill_histograms(event, matched_jet, parton_jet, "matched","_fpt");}

    // --------apply Top Tag-----------------
    if(toptagger->Is_tagged("sd", matched_jets[j])){
      matched_jets_tagged.push_back(matched_jet);
      matched_pair_tagged.push_back(matched_pair[j]);
    
      // fill hists with tagged jets
      double tau32 = matched_jet.tau3_groomed()/matched_jet.tau2_groomed();
      fill_histograms(event, matched_jet, parton_jet, "tagged",""); 
      if(tau32<0.4) fill_histograms(event, matched_jet, parton_jet, "tagged","tau32"); 

    } // end if top tag
  } // end loop over matched jets

  //set the event handle  these jets are the matched TopJets
  //  event.set(h_matched_parton_jets, matched_parton_jets);
  //  event.set(h_matched_jets, matched_jets);
  
  //for ROC curves
  //set the event handle for the denominator
  // event.set(h_matched_pairs, matched_pair);
  // hist_hotvr_pre->fill(event);   // fill hists before tau selection (all parton jets)
  // hist_eff_hotvr->fill(event);   // fill hists for efficiency plots
  // //set the event handle for the numerator
  // event.set(h_matched_pairs, matched_pair_tagged);
  // for (int i = 0; i < n_points; ++i){ //loop over working points
  //   if (sel_hotvr.at(i)->passes(event)){
  //     hist_hotvr.at(i)->fill(event); //fill hist if at least one jet of the event passes the sel
  //   }
  // }

  // decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the HOTVRStudiesModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(HOTVRStudiesModule)
