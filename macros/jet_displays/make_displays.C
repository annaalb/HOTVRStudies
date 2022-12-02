#include <TString.h>
#include <TSystem.h>
#include <algorithm>
#include <vector>
#include <iostream>

#include "style.h"
#include "jet_display_helpers.h"
#include "jet_display_hotvr.h"
#include "jet_display_akt.h"

void do_displays(TString prefix);

//TString gdir = "NoCutsHistosGen_mass/"; // use this for HOTVR display
//TString gdir = "NoCutsHistos3Gen/";       // anti-kt jets

int main()
{
    gROOT->SetBatch(true);

    std::vector<TString> prefixes; 
    //prefixes.push_back(TString("history_rho"));
    //prefixes.push_back(TString("history_fancy"));
    prefixes.push_back(TString("history_ghosts"));
    //prefixes.push_back(TString("fancy_R_pt_protection_for_ghosts"));
    //prefixes.push_back(TString("fancyR_wo_turnon_fix_mass"));
    //prefixes.push_back(TString("fancyR_with_turnon_beta0"));
    //prefixes.push_back(TString("fancyR_with_turnon_grooming_fix_mass"));
    //prefixes.push_back(TString("fancyR_with_turnon"));
    //prefixes.push_back(TString("fancyR_with_turnon_R0_1"));
    //prefixes.push_back(TString("fancyR_with_turnon_no_grooming"));    

    for (int i=0; i<prefixes.size(); ++i){
      do_displays(prefixes[i]);
    }

    exit(0);
}

void do_displays(TString prefix)
{
    gROOT->SetBatch(true);

    // 50 events for high and low pt:
    //TFile* file = new TFile("JetDisplayCycle.MC.TTbarEff_gen_tight_selection_subjets.root", "READ");

    // high resolution, high pt
    //   TFile* file = new TFile("JetDisplayCycle.MC.TTbarEff_gen_tight_selection_highres_highpt.root", "READ");
    //TFile* file = new TFile("JetDisplayCycle.MC.TTbarEff_gen_softdrop.root", "READ");
    
    //TFile* file = new TFile("files/HOTVR_SD_fancy_R/uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23.root", "READ");

    TString fname = TString("files/") + prefix + TString("/uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23.root");
    std::cout << "file name = " << fname << std::endl;
    TFile* file = new TFile(fname, "READ");
    //TFile* file = new TFile("files/debug_fancyR/fancyR_with_turnon_no_grooming/uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23.root", "READ");    
    //TFile* file = new TFile("files/new_uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23.root", "READ");
    // high resolution, low pt
    //TFile* file = new TFile("JetDisplayCycle.MC.TTbarEff_gen_tight_selection_highres2.root", "READ");    

  	// loop over directories in file and plot jets for each directory
  	// (each holds one event)
  	file->cd();

  	std::vector<TString> dirs;

  	TKey *key;
  	TIter nextkey( gDirectory->GetListOfKeys() );
  	while ( (key = (TKey*)nextkey())) {

      TObject *obj = key->ReadObj();
      if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
        TString name = obj->GetName();
        if (name.Contains("JetDisplay")) dirs.push_back(name);
      }
    }
	
    // before plotting: delete ROOT's random number generator
    gRandom = NULL;
	  for (Int_t i=0; i<dirs.size(); ++i){
      UInt_t ievent = GetEventNumber(dirs[i]);
      //if (ievent!=0) continue;
      NewEvent(ievent);
		  Int_t num = ievent;
      if (dirs[i].Contains("_akt_")){
        jet_display_akt(file, dirs[i], prefix, num);
      } else if (dirs[i].Contains("_ca_")) {
        jet_display_akt(file, dirs[i], prefix, num);
      } else {
        TString fname = file->GetName();
        if (fname.Contains("no_grooming")){
          std::cout << "\nPlot only large jets without subjets (no grooming case)." << std::endl;
          jet_display_hotvr(file, dirs[i], prefix, num, true); 
        } else {
          std::cout << "\nPlot large jets with subjets (grooming case)" << std::endl;
          jet_display_hotvr(file, dirs[i], prefix, num, false); 
        }
      }
	  }

}

