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

class HOTVRClusteringModule: public AnalysisModule {
public:

  explicit HOTVRClusteringModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
//initialize hist classes
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_200;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_400;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_600;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_800;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_1000;

  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub3;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub3_200;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub3_400;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub3_600;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub3_800;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub3_1000;

  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub1;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub1_200;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub1_400;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub1_600;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub1_800;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_Nsub1_1000;

  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_fpt;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_fpt_200;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_fpt_400;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_fpt_600;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_fpt_800;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_fpt_1000;

  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_mass;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_mass_200;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_mass_400;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_mass_600;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_mass_800;
  std::unique_ptr<HOTVRJetsHists> hist_hotvr_jets_mass_1000;

// initialize event handle
  Event::Handle<vector<TopJet>> h_HOTVR_jets;
  Event::Handle<vector<TopJet>> h_parton_jets;

// for performance study
  std::unique_ptr<AnalysisModule> cl_topjet_hotvr;
  std::unique_ptr<Selection> sel_ntop_hotvr;
  std::unique_ptr<Selection> sel_toptag_hotvr, sel_toptag_softdrop;


// initialize classes
  string m_clustering;
  Matching* matching;
  Clustering* clustering;

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
// ... containing UHH2(top)jets (in order to save the tau32 values)
  vector<TopJet> _top_hotvr_jets;
  vector<TopJet> _top_vr_jets;
  vector<TopJet> _top_parton_jets;

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

/*
 ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
 ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

HOTVRClusteringModule::HOTVRClusteringModule(Context & ctx){
  cout << "Starting HOTVRClusteringModule!" << endl;
  if(debug){cout << "HOTVRClusteringModule: Debugging mode :) " << '\n';}
// get info from xml
  m_clustering = ctx.get("Clustering");
// check for the dataset version (ttbar or QCD)
  dataset_version = ctx.get("dataset_version");
  isTTbar = dataset_version.find("ttbar") == 0;
  is_mc = ctx.get("dataset_type") == "MC";
  is_qcd = (dataset_version.find("QCD") == 0);

  ctx.undeclare_all_event_output();

  h_HOTVR_jets = ctx.declare_event_output<vector<TopJet>>("HOTVR_jets");
  h_parton_jets = ctx.declare_event_output<vector<TopJet>>("parton_jets");

//declare event handle
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

  hist_hotvr_jets_Nsub3.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub3", is_qcd));
  hist_hotvr_jets_Nsub3_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub3_200", is_qcd));
  hist_hotvr_jets_Nsub3_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub3_400", is_qcd));
  hist_hotvr_jets_Nsub3_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub3_600", is_qcd));
  hist_hotvr_jets_Nsub3_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub3_800", is_qcd));
  hist_hotvr_jets_Nsub3_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub3_1000", is_qcd));

  hist_hotvr_jets_Nsub1.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub1", is_qcd));
  hist_hotvr_jets_Nsub1_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub1_200", is_qcd));
  hist_hotvr_jets_Nsub1_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub1_400", is_qcd));
  hist_hotvr_jets_Nsub1_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub1_600", is_qcd));
  hist_hotvr_jets_Nsub1_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub1_800", is_qcd));
  hist_hotvr_jets_Nsub1_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_Nsub1_1000", is_qcd));

  hist_hotvr_jets_fpt.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_fpt", is_qcd));
  hist_hotvr_jets_fpt_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_fpt_200", is_qcd));
  hist_hotvr_jets_fpt_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_fpt_400", is_qcd));
  hist_hotvr_jets_fpt_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_fpt_600", is_qcd));
  hist_hotvr_jets_fpt_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_fpt_800", is_qcd));
  hist_hotvr_jets_fpt_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_fpt_1000", is_qcd));

  hist_hotvr_jets_mass.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_mass", is_qcd));
  hist_hotvr_jets_mass_200.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_mass_200", is_qcd));
  hist_hotvr_jets_mass_400.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_mass_400", is_qcd));
  hist_hotvr_jets_mass_600.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_mass_600", is_qcd));
  hist_hotvr_jets_mass_800.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_mass_800", is_qcd));
  hist_hotvr_jets_mass_1000.reset(new HOTVRJetsHists(ctx, "HOTVRJetsHists_hotvr_jets_mass_1000", is_qcd));

}
/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/
bool HOTVRClusteringModule::process(Event & event) {
  if(debug){std::cout << "Begin process..." << '\n';}
// here define genparticles and read from event
  auto genparticles = event.genparticles;

  if(isTTbar){
    if (debug) {std::cout << "Event is ttbar!" << '\n';}
    ttgenprod->process(event);
    const auto & ttbargen = event.get(h_ttbargen);
    if (!ttbargen.IsTopHadronicDecay() || !ttbargen.IsAntiTopHadronicDecay()) { // skip leptonic decays
      return false;
    }
  }

// //MATCHING -> class to filter the stable particles from the genparticles and create a pseudojet from them
  if (debug) {std::cout << "Start setting pseudojets from stable particles!" << '\n';}
   matching = new Matching();
   if (debug) {std::cout << "Genparticles size "<< genparticles->size() << '\n';}
   matching->set_stable_particles(genparticles);
   pseudojets = matching->get_stableParticles();
   matching->set_partons(genparticles);
   parton_pseudojets = matching->get_partons();
//
// //CLUSTERING -> clusters jets depending on the algorithm that is chosen
 if (debug) {std::cout << "Start Clustering!" << '\n';}
   clustering = new Clustering(m_clustering);
  clustering->cluster_jets(pseudojets); // cluster the pseudojets, possible modes defined in hotvr.config: "HOTVR, HOTVR_SD, VR"
  hotvr_jets = clustering->get_hotvr_jets();
  //hotvr_jets_constituents = clustering->get_hotvr_jet_constituents();

  // cluster and get the clustered parton jets (AK10 jets)
  clustering->cluster_parton_jets(parton_pseudojets, isTTbar);
  parton_jets = clustering->get_parton_pseudojets();
  //get the topjets from the Clustering
  _top_hotvr_jets=clustering->get_top_hotvr_jets();
  _top_parton_jets=clustering->get_parton_jets();

//fill hists with hotvr jets and corresponding rejected subjets
  hist_hotvr_jets->fill_n_jets(event, _top_hotvr_jets);
  if (debug) {std::cout << "Now fill some control hists!" << '\n';}

  for(uint j=0; j<_top_hotvr_jets.size(); ++j){ // loop over hotvr jets
    //std::cout << "Jet pt "<< _top_hotvr_jets[j].pt() << '\n';
    hist_hotvr_jets->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>200 && _top_hotvr_jets[j].pt()<400)  hist_hotvr_jets_200->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>400 && _top_hotvr_jets[j].pt()<600)  hist_hotvr_jets_400->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>600 && _top_hotvr_jets[j].pt()<800)  hist_hotvr_jets_600->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>800 && _top_hotvr_jets[j].pt()<1000)  hist_hotvr_jets_800->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>1000 && _top_hotvr_jets[j].pt()<1200)  hist_hotvr_jets_1000->fill_topjet(event, _top_hotvr_jets[j]);
// fill hists after Nsub3 cut
if (_top_hotvr_jets[j].subjets().size()>2) {
    hist_hotvr_jets_Nsub3->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>200 && _top_hotvr_jets[j].pt()<400)  hist_hotvr_jets_Nsub3_200->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>400 && _top_hotvr_jets[j].pt()<600)  hist_hotvr_jets_Nsub3_400->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>600 && _top_hotvr_jets[j].pt()<800)  hist_hotvr_jets_Nsub3_600->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>800 && _top_hotvr_jets[j].pt()<1000)  hist_hotvr_jets_Nsub3_800->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>1000 && _top_hotvr_jets[j].pt()<1200)  hist_hotvr_jets_Nsub3_1000->fill_topjet(event, _top_hotvr_jets[j]);
  }
  // fill hists for jets with 1 subjet
  if (_top_hotvr_jets[j].subjets().size()==1) {
    if (debug && _top_hotvr_jets[j].v4().M()>40) {
      std::cout << "WARNING: HOTVRClusteringModule - nsubjets == 1 but mjet > 40 " << '\n';
    }
      hist_hotvr_jets_Nsub1->fill_topjet(event, _top_hotvr_jets[j]);
      if(_top_hotvr_jets[j].pt()>200 && _top_hotvr_jets[j].pt()<400)  hist_hotvr_jets_Nsub1_200->fill_topjet(event, _top_hotvr_jets[j]);
      if(_top_hotvr_jets[j].pt()>400 && _top_hotvr_jets[j].pt()<600)  hist_hotvr_jets_Nsub1_400->fill_topjet(event, _top_hotvr_jets[j]);
      if(_top_hotvr_jets[j].pt()>600 && _top_hotvr_jets[j].pt()<800)  hist_hotvr_jets_Nsub1_600->fill_topjet(event, _top_hotvr_jets[j]);
      if(_top_hotvr_jets[j].pt()>800 && _top_hotvr_jets[j].pt()<1000)  hist_hotvr_jets_Nsub1_800->fill_topjet(event, _top_hotvr_jets[j]);
      if(_top_hotvr_jets[j].pt()>1000 && _top_hotvr_jets[j].pt()<1200)  hist_hotvr_jets_Nsub1_1000->fill_topjet(event, _top_hotvr_jets[j]);
    }
// fill hists after fpt cut
  if (_top_hotvr_jets[j].subjets().size()>2 && _top_hotvr_jets[j].hotvr_fpt1()<0.8) {
    hist_hotvr_jets_fpt->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>200 && _top_hotvr_jets[j].pt()<400)  hist_hotvr_jets_fpt_200->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>400 && _top_hotvr_jets[j].pt()<600)  hist_hotvr_jets_fpt_400->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>600 && _top_hotvr_jets[j].pt()<800)  hist_hotvr_jets_fpt_600->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>800 && _top_hotvr_jets[j].pt()<1000)  hist_hotvr_jets_fpt_800->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>1000 && _top_hotvr_jets[j].pt()<1200)  hist_hotvr_jets_fpt_1000->fill_topjet(event, _top_hotvr_jets[j]);
  }
// fill hist after mass cut
  if (_top_hotvr_jets[j].subjets().size()>2 && 140<_top_hotvr_jets[j].v4().M() && _top_hotvr_jets[j].v4().M()<220) {
    hist_hotvr_jets_mass->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>200 && _top_hotvr_jets[j].pt()<400)  hist_hotvr_jets_mass_200->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>400 && _top_hotvr_jets[j].pt()<600)  hist_hotvr_jets_mass_400->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>600 && _top_hotvr_jets[j].pt()<800)  hist_hotvr_jets_mass_600->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>800 && _top_hotvr_jets[j].pt()<1000)  hist_hotvr_jets_mass_800->fill_topjet(event, _top_hotvr_jets[j]);
    if(_top_hotvr_jets[j].pt()>1000 && _top_hotvr_jets[j].pt()<1200)  hist_hotvr_jets_mass_1000->fill_topjet(event, _top_hotvr_jets[j]);
  }
} // end loop over hotvr jets

//set the event handle (store HOTVR jets and parton jets)
  event.set(h_HOTVR_jets, _top_hotvr_jets);
  event.set(h_parton_jets, _top_parton_jets);
  //
  // // delete clustering infos
   clustering->Reset();
// decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the HOTVRClusteringModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(HOTVRClusteringModule)
