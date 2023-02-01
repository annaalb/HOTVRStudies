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
#include "UHH2/HOTVRStudies/include/HOTVRWtopHists.h"
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

class HOTVRWTopMatchingModule: public AnalysisModule {
public:

  explicit HOTVRWTopMatchingModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
//initialize hist classes
std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets;
std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_200;
std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_400;
std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_600;
std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_800;
std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_1000;

std::unique_ptr<HOTVRWtopHists> hist_W_top_parton_jets;
std::unique_ptr<HOTVRWtopHists> hist_W_top_parton_jets_200;
std::unique_ptr<HOTVRWtopHists> hist_W_top_parton_jets_400;
std::unique_ptr<HOTVRWtopHists> hist_W_top_parton_jets_600;
std::unique_ptr<HOTVRWtopHists> hist_W_top_parton_jets_800;
std::unique_ptr<HOTVRWtopHists> hist_W_top_parton_jets_1000;

std::unique_ptr<HOTVRWtopHists> hist_W_b_matched_parton_jets;
std::unique_ptr<HOTVRWtopHists> hist_W_matched_parton_jets;

std::unique_ptr<HOTVRJetsHists> hist_W_parton_jets;
std::unique_ptr<HOTVRJetsHists> hist_W_parton_jets_200;
std::unique_ptr<HOTVRJetsHists> hist_W_parton_jets_400;
std::unique_ptr<HOTVRJetsHists> hist_W_parton_jets_600;
std::unique_ptr<HOTVRJetsHists> hist_W_parton_jets_800;
std::unique_ptr<HOTVRJetsHists> hist_W_parton_jets_1000;

std::unique_ptr<HOTVRJetsHists> hist_b_parton_jets;
std::unique_ptr<HOTVRJetsHists> hist_b_parton_jets_200;
std::unique_ptr<HOTVRJetsHists> hist_b_parton_jets_400;
std::unique_ptr<HOTVRJetsHists> hist_b_parton_jets_600;
std::unique_ptr<HOTVRJetsHists> hist_b_parton_jets_800;
std::unique_ptr<HOTVRJetsHists> hist_b_parton_jets_1000;

std::unique_ptr<HOTVRJetsHists> hist_top_parton_jets;
std::unique_ptr<HOTVRJetsHists> hist_top_parton_jets_200;
std::unique_ptr<HOTVRJetsHists> hist_top_parton_jets_400;
std::unique_ptr<HOTVRJetsHists> hist_top_parton_jets_600;
std::unique_ptr<HOTVRJetsHists> hist_top_parton_jets_800;
std::unique_ptr<HOTVRJetsHists> hist_top_parton_jets_1000;

  std::unique_ptr<HOTVRWtopHists> hist_matched_jets;
  std::unique_ptr<HOTVRWtopHists> hist_matched_jets_200;
  std::unique_ptr<HOTVRWtopHists> hist_matched_jets_400;
  std::unique_ptr<HOTVRWtopHists> hist_matched_jets_600;
  std::unique_ptr<HOTVRWtopHists> hist_matched_jets_800;
  std::unique_ptr<HOTVRWtopHists> hist_matched_jets_1000;

// initialize event handle
  Event::Handle<vector<TopJet>> h_HOTVR_jets;
  Event::Handle<vector<TopJet>> h_parton_jets;
  Event::Handle<vector<TopJet>> h_W_parton_jets;

  Event::Handle<vector<TopJet>> h_matched_jets;
  Event::Handle<vector<TopJet>> h_matched_parton_jets;
  Event::Handle<vector<pair<TopJet, TopJet>>> h_matched_pairs; //contains pair of matched jets
    // denominator for the ROC curves first: ttbar, second:qcd

// initialize classes
  string m_clustering;
  Matching* matching;
  Clustering* clustering;
  TopTagger* toptagger;

// initialize vectors of jets
// ... containing pseudojets
  vector<fastjet::PseudoJet> pseudojets;
  vector<fastjet::PseudoJet> hotvr_jets;
  vector<fastjet::PseudoJet> parton_pseudojets;
  vector<fastjet::PseudoJet> Wb_parton_pseudojets;

  vector<fastjet::PseudoJet> parton_jets;
  vector<fastjet::PseudoJet> W_parton_jets;
  vector<fastjet::PseudoJet> b_parton_jets;

  vector<TopJet> _top_hotvr_jets;
  vector<TopJet> _top_parton_jets;
  vector<TopJet> _top_parton_jets_W;
  vector<TopJet> _top_parton_jets_b;

  string dataset_version;

  bool isTTbar;
  bool is_mc, is_qcd;

  int nevent=0;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  Event::Handle<TTbarGen> h_ttbargen;
  bool debug = false;
};

/*
 ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
 ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

HOTVRWTopMatchingModule::HOTVRWTopMatchingModule(Context & ctx){
  cout << "Starting HOTVRWTopMatchingModule!" << endl;
  if(debug){cout << "HOTVRWTopMatchingModule: Debugging mode :) " << '\n';}
// get info from xml
  m_clustering = ctx.get("Clustering");
// check for the dataset version (ttbar or QCD)
  dataset_version = ctx.get("dataset_version");
  isTTbar = dataset_version.find("ttbar") == 0;
  is_mc = ctx.get("dataset_type") == "MC";
  is_qcd = (dataset_version.find("QCD") == 0);
  if(debug){std::cout << "Got info from config file" << '\n';}

  h_HOTVR_jets = ctx.get_handle<vector<TopJet>>("HOTVR_jets");
  h_parton_jets = ctx.get_handle<vector<TopJet>>("parton_jets");
  h_W_parton_jets = ctx.get_handle<vector<TopJet>>("W_parton_jets");

//declare event handle
  h_matched_jets = ctx.get_handle<vector<TopJet>>("matched_jets");
  h_matched_parton_jets = ctx.get_handle<vector<TopJet>>("matched_parton_jets");
  h_matched_pairs = ctx.get_handle<vector<pair<TopJet, TopJet>>>("matched_pairs");

  if(debug){std::cout << "Set handle" << '\n';}

  if (isTTbar) {
    if(debug){std::cout << "HOTVRClusteringModule: Event is ttbar!" << '\n';}
    const std::string ttbar_gen_label("ttbargen");
    ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, true));
    h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
  }
// Set up Hists classes:
hist_hotvr_jets.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets", is_qcd));
hist_hotvr_jets_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_200", is_qcd));
hist_hotvr_jets_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_400", is_qcd));
hist_hotvr_jets_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_600", is_qcd));
hist_hotvr_jets_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_800", is_qcd));
hist_hotvr_jets_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_1000", is_qcd));

hist_W_top_parton_jets.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_W_top_parton_jets", is_qcd));
hist_W_top_parton_jets_200.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_W_top_parton_jets_200", is_qcd));
hist_W_top_parton_jets_400.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_W_top_parton_jets_400", is_qcd));
hist_W_top_parton_jets_600.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_W_top_parton_jets_600", is_qcd));
hist_W_top_parton_jets_800.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_W_top_parton_jets_800", is_qcd));
hist_W_top_parton_jets_1000.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_W_top_parton_jets_1000", is_qcd));

hist_W_b_matched_parton_jets.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_W_b_matched_parton_jets", is_qcd));
hist_W_matched_parton_jets.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_W_matched_parton_jets", is_qcd));

hist_W_parton_jets.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_W_parton_jets", is_qcd));
hist_W_parton_jets_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_W_parton_jets_200", is_qcd));
hist_W_parton_jets_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_W_parton_jets_400", is_qcd));
hist_W_parton_jets_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_W_parton_jets_600", is_qcd));
hist_W_parton_jets_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_W_parton_jets_800", is_qcd));
hist_W_parton_jets_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_W_parton_jets_1000", is_qcd));

hist_b_parton_jets.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_b_parton_jets", is_qcd));
hist_b_parton_jets_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_b_parton_jets_200", is_qcd));
hist_b_parton_jets_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_b_parton_jets_400", is_qcd));
hist_b_parton_jets_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_b_parton_jets_600", is_qcd));
hist_b_parton_jets_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_b_parton_jets_800", is_qcd));
hist_b_parton_jets_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_b_parton_jets_1000", is_qcd));

hist_top_parton_jets.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_top_parton_jets", is_qcd));
hist_top_parton_jets_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_top_parton_jets_200", is_qcd));
hist_top_parton_jets_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_top_parton_jets_400", is_qcd));
hist_top_parton_jets_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_top_parton_jets_600", is_qcd));
hist_top_parton_jets_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_top_parton_jets_800", is_qcd));
hist_top_parton_jets_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_top_parton_jets_1000", is_qcd));

  hist_matched_jets.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_matched_jets", is_qcd));
  hist_matched_jets_200.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_matched_jets_200", is_qcd));
  hist_matched_jets_400.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_matched_jets_400", is_qcd));
  hist_matched_jets_600.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_matched_jets_600", is_qcd));
  hist_matched_jets_800.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_matched_jets_800", is_qcd));
  hist_matched_jets_1000.reset(new HOTVRWtopHists(ctx, "HOTVRWtopHists_matched_jets_1000", is_qcd));
  if(debug){std::cout << "End Module Constructor" << '\n';}
}
/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/
bool HOTVRWTopMatchingModule::process(Event & event) {
  if(debug){std::cout << "Begin process...!! event " << nevent << '\n';}
  ++nevent;
  auto genparticles = event.genparticles;

  if(isTTbar){
    if (debug) {std::cout << "Event is ttbar!" << '\n';}
    ttgenprod->process(event);
    const auto & ttbargen = event.get(h_ttbargen);
    if (!ttbargen.IsTopHadronicDecay() || !ttbargen.IsAntiTopHadronicDecay()) { // skip leptonic decays
      return false;
    }
  }
  // ------MATCHING get pseudojet collections for the clustering--------
    matching = new Matching();
    if (debug) {std::cout << "New Matching class" << '\n';}

    matching->set_stable_particles(genparticles); // Input for HOTVR jet clustering
    pseudojets = matching->get_stableParticles();
    if (debug) {std::cout << "Get stable particles" << '\n';}

    matching->set_partons(genparticles); // Input for top parton jet clustering
    matching->set_Wb_partons_from_top(genparticles); // Input for W and b parton jets

    parton_pseudojets = matching->get_partons();
    Wb_parton_pseudojets = matching->get_Wb_partons();

    if (debug) {
      std::cout << "Set and get partons" << '\n';
      std::cout << "Wb partons size "<< Wb_parton_pseudojets.size() << '\n';
    }

    // -----CLUSTERING -----
    clustering = new Clustering(m_clustering);
    clustering->cluster_jets(pseudojets); // cluster the pseudojets, possible modes defined in hotvr.config: "HOTVR, HOTVR_SD, VR"
    hotvr_jets = clustering->get_hotvr_jets();
    _top_hotvr_jets=clustering->get_top_hotvr_jets();
    if (debug) {std::cout << "Cluster HOTVR jets " << '\n';}

    for (size_t i = 0; i < _top_hotvr_jets.size(); i++) {
      TopJet jet = _top_hotvr_jets[i];
      hist_hotvr_jets->fill_topjet(event, jet);
      if(jet.pt()>200 &&jet.pt()<400)  hist_hotvr_jets_200->fill_topjet(event, jet);
      if(jet.pt()>400 &&jet.pt()<600)  hist_hotvr_jets_400->fill_topjet(event, jet);
      if(jet.pt()>600 &&jet.pt()<800)  hist_hotvr_jets_600->fill_topjet(event, jet);
      if(jet.pt()>800 &&jet.pt()<1000)  hist_hotvr_jets_800->fill_topjet(event, jet);
      if(jet.pt()>1000 &&jet.pt()<1200)  hist_hotvr_jets_1000->fill_topjet(event, jet);
    }
   // cluster W parton GENJETS
   clustering->cluster_parton_jets(Wb_parton_pseudojets, 24);
   W_parton_jets = clustering->get_parton_pseudojets();
   _top_parton_jets_W = clustering->get_parton_jets();
   if (debug) {
     std::cout << "Cluster W jets " << '\n';
     std::cout << "Number of W jets: " << _top_parton_jets_W.size() << '\n';
   }
   // cluster b parton GENJETS
   clustering->ResetPartonJets();
   clustering->cluster_parton_jets(Wb_parton_pseudojets, 5);
   b_parton_jets = clustering->get_parton_pseudojets();
   _top_parton_jets_b = clustering->get_parton_jets();
   if (debug) {
     std::cout << "Cluster b jets " << '\n';
     std::cout << b_parton_jets.size() << '\n';
   }

   clustering->ResetPartonJets();
   uint keepID=0;
   if (isTTbar){
     keepID=6;
   }
   clustering->cluster_parton_jets(parton_pseudojets, keepID);
   parton_jets = clustering->get_parton_pseudojets();
   _top_parton_jets = clustering->get_parton_jets();
   if (debug) {std::cout << "Cluster top jets " << '\n';}
   // Plot control hists for the 3 collections of parton jets
    for (size_t i = 0; i < _top_parton_jets_W.size(); i++) {
      TopJet jet = _top_parton_jets_W[i];
      hist_W_parton_jets->fill_topjet(event, jet);
      if(jet.pt()>200 &&jet.pt()<400)  hist_W_parton_jets_200->fill_topjet(event, jet);
      if(jet.pt()>400 &&jet.pt()<600)  hist_W_parton_jets_400->fill_topjet(event, jet);
      if(jet.pt()>600 &&jet.pt()<800)  hist_W_parton_jets_600->fill_topjet(event, jet);
      if(jet.pt()>800 &&jet.pt()<1000)  hist_W_parton_jets_800->fill_topjet(event, jet);
      if(jet.pt()>1000 &&jet.pt()<1200)  hist_W_parton_jets_1000->fill_topjet(event, jet);
    }
    for (size_t i = 0; i < _top_parton_jets_b.size(); i++) {
      TopJet jet = _top_parton_jets_b[i];
      hist_b_parton_jets->fill_topjet(event, jet);
      if(jet.pt()>200 &&jet.pt()<400)  hist_b_parton_jets_200->fill_topjet(event, jet);
      if(jet.pt()>400 &&jet.pt()<600)  hist_b_parton_jets_400->fill_topjet(event, jet);
      if(jet.pt()>600 &&jet.pt()<800)  hist_b_parton_jets_600->fill_topjet(event, jet);
      if(jet.pt()>800 &&jet.pt()<1000)  hist_b_parton_jets_800->fill_topjet(event, jet);
      if(jet.pt()>1000 &&jet.pt()<1200)  hist_b_parton_jets_1000->fill_topjet(event, jet);
    }
    for (size_t i = 0; i < _top_parton_jets.size(); i++) {
      TopJet jet = _top_parton_jets[i];
      hist_top_parton_jets->fill_topjet(event, jet);
      if(jet.pt()>200 &&jet.pt()<400)  hist_top_parton_jets_200->fill_topjet(event, jet);
      if(jet.pt()>400 &&jet.pt()<600)  hist_top_parton_jets_400->fill_topjet(event, jet);
      if(jet.pt()>600 &&jet.pt()<800)  hist_top_parton_jets_600->fill_topjet(event, jet);
      if(jet.pt()>800 &&jet.pt()<1000)  hist_top_parton_jets_800->fill_topjet(event, jet);
      if(jet.pt()>1000 &&jet.pt()<1200)  hist_top_parton_jets_1000->fill_topjet(event, jet);
    }
// // ------MATCHING zum W--------
// //run_matching: loop over the parton jets and match them to the hotvr jets
   TopJet closest_top_jet;

   // for HOTVR jets
   for (size_t i = 0; i < _top_hotvr_jets.size(); i++) {
     // for each jet
     TopJet jet = _top_hotvr_jets[i];
     if (jet.pt()>30) {
     // // get closest W jet
     TopJet closest_W_jet = matching->get_closest_jet(jet, _top_parton_jets_W);
     // // get closest top jet
     closest_top_jet = matching->get_closest_jet(jet, _top_parton_jets);
     // // if matched to W jet -> find closest b jet
     if (debug) {
       std::cout << " ------- HOTVR jet pt = " << jet.pt() << '\n';
       std::cout << " ------- HOTVR jet mass = " << jet.v4().mass() << '\n';

       std::cout << "closest W jet pt = " << closest_W_jet.pt() << '\n';
       std::cout << "closest top jet pt = " << closest_top_jet.pt() << '\n';
     }

     if (closest_W_jet.pt() > 0) {
         TopJet closest_b_jet = matching->get_closest_jet(jet, _top_parton_jets_b);
         if(debug){
           std::cout << "closest b jet pt = " << closest_b_jet.pt() << '\n';
           std::cout << "Delta R between top and b " << deltaR(closest_top_jet, closest_b_jet) << '\n';

       }
    //  //   // if no b jet -> fill hists
         if (deltaR(jet, closest_b_jet)<1) { // distance(jet, b) < 1 ??
           if(debug){std::cout << "matched to b jet" << '\n';}
    //  //     // fill hist for W and top parton jet
           hist_W_b_matched_parton_jets->fill_topjet(event, jet, closest_W_jet, closest_b_jet, closest_top_jet);   // kleinsten Abstand W und HOTVR plotten, pt ratio plotten
         }
    //  //   // if b jet -> fill hists
         else{
           if(debug){
             std::cout << "NOT matched to b jet" << '\n';
             if (deltaR(closest_top_jet, closest_b_jet)<1) {
               std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << '\n';
               std::cout << "Not matched to b, but dR(top,b)<1. Event number: " << nevent << '\n';
               std::cout << "HOTVR jet mass=" << jet.v4().mass() << '\n';
               std::cout << "Abstand W b " << deltaR(closest_W_jet, closest_b_jet) << '\n';
             }
           }
          hist_W_matched_parton_jets->fill_topjet(event, jet, closest_W_jet, closest_b_jet, closest_top_jet);   // kleinsten Abstand W und HOTVR jet plotten, pt ratio plotten
         }
       }// end if matched to W
    }// jet pt > 30
   }



// // delete clustering infos
 clustering->Reset();
// decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the HOTVRWTopMatchingModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(HOTVRWTopMatchingModule)
