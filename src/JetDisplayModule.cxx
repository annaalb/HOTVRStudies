#include <iostream>
#include <memory>
// include general UHH2 classes
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/common/include/PrintingModules.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TTbarGen.h"
// include everything from the HOTVR Studies
#include "UHH2/HOTVRStudies/include/Matching.h"
#include "UHH2/HOTVRStudies/include/Clustering.h"
#include "UHH2/HOTVRStudies/include/JetDisplayHists.h"
#include "UHH2/HOTVRStudies/include/TopTagger.h"
// include from fastjet
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/HOTVRinfo.hh"
#include "fastjet/contrib/HOTVR.hh"

using namespace std;
using namespace uhh2;
using namespace fastjet;
using namespace contrib;

/** AnalysisModule for HOTVR Studies to create Jet Displays.
*
* This is the central class which calls other AnalysisModules, Hists or Selection classes.
* This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
*/

class JetDisplayModule: public AnalysisModule {
public:

  explicit JetDisplayModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
//initialize hist classes
std::unique_ptr<JetDisplayHists> hist_jet_display_00;
std::unique_ptr<JetDisplayHists> hist_jet_display_01;
std::unique_ptr<JetDisplayHists> hist_jet_display_02;
std::unique_ptr<JetDisplayHists> hist_jet_display_03;
std::unique_ptr<JetDisplayHists> hist_jet_display_04;
std::unique_ptr<JetDisplayHists> hist_jet_display_05;
std::unique_ptr<JetDisplayHists> hist_jet_display_06;
std::unique_ptr<JetDisplayHists> hist_jet_display_07;
std::unique_ptr<JetDisplayHists> hist_jet_display_08;
std::unique_ptr<JetDisplayHists> hist_jet_display_09;

std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
std::unique_ptr<uhh2::AnalysisModule> genprinter;

// initialize event handle
Event::Handle<vector<PseudoJet>> h_parts;

Event::Handle<vector<PseudoJet>> h_jets;
Event::Handle<vector<PseudoJet>> h_rejected_cluster;
Event::Handle<vector<PseudoJet>> h_soft_cluster;
Event::Handle<vector<PseudoJet>> h_rejected_subjets;

Event::Handle<TTbarGen> h_ttbargen;

// initialize classes
  string m_clustering;
  Matching* matching;
  Clustering* clustering;
  TopTagger* toptagger;

// initialize vectors of jets containing pseudojets
  vector<fastjet::PseudoJet> parts;

  vector<fastjet::PseudoJet> hotvr_jets;
  vector<fastjet::PseudoJet> rejected_cluster;
  vector<fastjet::PseudoJet> soft_cluster;
  vector<fastjet::PseudoJet> rejected_subjets;

  string dataset_version;

  bool isTTbar;
  bool is_mc, is_qcd;

  uint nevent=0;
  uint _event_max = 10;
};

/*
 ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
 ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

JetDisplayModule::JetDisplayModule(Context & ctx){
  cout << "Starting JetDisplayModule!" << endl;
// get info from xml
  m_clustering = ctx.get("Clustering");
// check for the dataset version (ttbar or QCD)
  dataset_version = ctx.get("dataset_version");
  isTTbar = dataset_version.find("ttbar") == 0;
  is_mc = ctx.get("dataset_type") == "MC";
  is_qcd = (dataset_version.find("QCD") == 0);

  h_parts = ctx.get_handle<vector<PseudoJet>>("parts");

  h_jets = ctx.get_handle<vector<PseudoJet>>("jets");
  h_rejected_cluster = ctx.get_handle<vector<PseudoJet>>("rejected_cluster");
  h_soft_cluster = ctx.get_handle<vector<PseudoJet>>("soft_cluster");
  h_rejected_subjets = ctx.get_handle<vector<PseudoJet>>("rejected_subjets");

  const std::string ttbar_gen_label("ttbargen");
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, true));
  genprinter.reset(new GenParticlesPrinter(ctx));
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");

  hist_jet_display_00.reset( new JetDisplayHists(ctx, "JetDisplay_event00")); //TODO find elegant solution?
  hist_jet_display_01.reset( new JetDisplayHists(ctx, "JetDisplay_event01"));
  hist_jet_display_02.reset( new JetDisplayHists(ctx, "JetDisplay_event02"));
  hist_jet_display_03.reset( new JetDisplayHists(ctx, "JetDisplay_event03"));
  hist_jet_display_04.reset( new JetDisplayHists(ctx, "JetDisplay_event04"));
  hist_jet_display_05.reset( new JetDisplayHists(ctx, "JetDisplay_event05"));
  hist_jet_display_06.reset( new JetDisplayHists(ctx, "JetDisplay_event06"));
  hist_jet_display_07.reset( new JetDisplayHists(ctx, "JetDisplay_event07"));
  hist_jet_display_08.reset( new JetDisplayHists(ctx, "JetDisplay_event08"));
  hist_jet_display_09.reset( new JetDisplayHists(ctx, "JetDisplay_event09"));

  clustering = new Clustering(m_clustering);
  matching = new Matching();
  toptagger = new TopTagger();

}
/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/
bool JetDisplayModule::process(Event & event) {

  std::cout << "--------start processing event--------" << '\n';
  std::cout << "Nevent = " << nevent << std::endl;

  if (nevent<10) {
  // here define genparticles and read from event
  auto genparticles = event.genparticles;
  //genprinter->process(event);

  if(isTTbar){
    ttgenprod->process(event);
  }

  const auto & ttbargen = event.get(h_ttbargen);

  if (!ttbargen.IsTopHadronicDecay() || !ttbargen.IsAntiTopHadronicDecay()) {
    std::cout << "skip leptonic or semi leptonic ttbar decay " << '\n';
    return false;
  }

//MATCHING -> class to filter the stable particles from the genparticles and create a pseudojet from them
  matching->set_stable_particles(genparticles);
  parts = matching->get_stableParticles();

  std::cout << "parts in the event = " << parts.size() << std::endl;

//CLUSTERING -> clusters jets depending on the algorithm that is chosen
  clustering->cluster_jets(parts, nevent); // cluster the pseudojets, possible modes defined in hotvr.config: "HOTVR, HOTVR_SD, VR"

  // here are our HOTVR jets with subjets
  hotvr_jets = clustering->get_hotvr_jets();
  // get the rejected subjets
  rejected_subjets = clustering->get_rejected_subjets();
  rejected_cluster = clustering->get_rejected_cluster();
  soft_cluster = clustering->get_soft_cluster();

  event.set(h_jets, hotvr_jets);
  event.set(h_parts, parts);
  event.set(h_rejected_cluster, rejected_cluster);
  event.set(h_soft_cluster, soft_cluster);
  event.set(h_rejected_subjets, rejected_subjets);

  bool select_event = true;
  uint hotvr_jets_with_subjets=0;

  // for (size_t k = 0; k < hotvr_jets.size(); k++) {
  // //  if (toptagger->Is_tagged("sd", hotvr_jets[k]) && hotvr_jets[k].pt()>200) {
  //   if (hotvr_jets[k].pt()>200) {
  //     select_event=true;
  //   }
  //   //if (hotvr_jets[k].user_info<HOTVRinfo>().nsubjets() > 1){++hotvr_jets_with_subjets;}
  // }

  // if (hotvr_jets_with_subjets>9) {
  //   select_event=false;
  // }

  if(isTTbar && nevent < _event_max && select_event){
    std::cout << "Select event "<< nevent << '\n';
    if (nevent==0) { hist_jet_display_00->fill(event);}
    else if (nevent==1) { hist_jet_display_01->fill(event);}
    else if (nevent==2) { hist_jet_display_02->fill(event);}
    else if (nevent==3) { hist_jet_display_03->fill(event);}
    else if (nevent==4) { hist_jet_display_04->fill(event);}
    else if (nevent==5) { hist_jet_display_05->fill(event);}
    else if (nevent==6) { hist_jet_display_06->fill(event);}
    else if (nevent==7) { hist_jet_display_07->fill(event);}
    else if (nevent==8) { hist_jet_display_08->fill(event);}
    else if (nevent==9) { hist_jet_display_09->fill(event);}
    ++nevent;
  }

  // some couts to see what the jets look like
  cout << "Number of HOTVR jets = " << hotvr_jets.size() << endl;
  for (uint i=0; i<hotvr_jets.size(); ++i){
    cout << "jet " << i << " : " << endl;
    cout << "pt = " << hotvr_jets.at(i).pt() << " mass = " << hotvr_jets.at(i).m() << endl;
    std::vector<fastjet::PseudoJet> subjets = hotvr_jets.at(i).user_info<HOTVRinfo>().subjets();
    cout << "Nsubjets = " << subjets.size() << endl;
    for (uint j=0; j<subjets.size(); ++j){
      cout << "   subjet " << j << " : " << " pt = " << subjets.at(j).pt() << "  m = " << subjets.at(j).m() << endl;
    }
  }

  // delete clustering infos and clear input particles
  clustering->Reset();
  matching->Reset();

  } // end if nevent < 10

  // decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the JetDisplayModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(JetDisplayModule)
