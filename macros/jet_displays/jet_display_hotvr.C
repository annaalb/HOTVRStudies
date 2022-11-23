#include "jet_display_hotvr.h"
#include "jet_display_helpers.h"
#include "style.h"
#include <TSystem.h>
#include <iostream>
#include <vector>

//TString gdir = "NoCutsHistosGen_mass/"; // use this for HOTVR display
//TString gdir = "NoCutsHistos3Gen/";       // anti-kt jets

void jet_display_hotvr(TFile* file, TString dir, TString prefix, Int_t num, Bool_t bLargeJetsOnly)
{ 

  std::cout << "making jet display for HOTVR" << std::endl;

  SetupGlobalStyle(); 

  //TFile* file = new TFile("RocCycle.MC.TTbarEff_gen_jetdisplay2.root", "READ");
  //TFile* file = new TFile("genevt/RocCycle.MC.TTbarEff_gen_event86.root", "READ");  
  
  // store all jets with their subjets in this vector
  std::vector< std::vector<TH2F*> > subjets;
  subjets = GetSubjets(file, dir);
  for (int i=0; i<subjets.size(); ++i){
    CreateContours(subjets[i]);
  }

  // store all large jets
  std::vector<TH2F*> jets;
  jets = GetJets(file, dir);
  CreateContours(jets);  

  TString canname = TString::Format("canvas_%d", num);
  TCanvas *c = SetupCanvas2d(canname);

  // some colors
  std::vector< std::vector<int> > cols;
  cols = GetColours(); 

  // gray hues for rejected jets
  int lgray = TColor::GetColor( "#e6e6e6");
  int mgray = TColor::GetColor( "#dedede");
  int dgray = TColor::GetColor( "#808080");

  // here are the jets without subjets - ie the ones that have never fulfilled softdrop
  TH2F* jets_rejected;
  jets_rejected = GetJetsWithoutSubjets(file, dir);
  jets_rejected->GetXaxis()->SetTitle("#phi");
  jets_rejected->GetYaxis()->SetTitle("#eta");
  jets_rejected->SetTitle("");
  Cosmetics2d(jets_rejected);
  FixHistForContour(jets_rejected);
  DrawContour(jets_rejected, lgray);  

  // here are the subjets (and jets) with low pT
  //TH2F* subjets_rejected;
  //subjets_rejected = GetJetsWithLowPt(file, dir);

  // here are the rejected subjets/clusters by the SD criterion
  TH2F* subjets_rejected_sd;
  subjets_rejected_sd = GetSubjetsRejectedBySD(file, dir);
  if (subjets_rejected_sd){
    Cosmetics2d(subjets_rejected_sd);
    FixHistForContour(subjets_rejected_sd);
    DrawContour(subjets_rejected_sd, dgray);   
  } 


  if (bLargeJetsOnly){
    // draw large jets
    for (int i=0; i<jets.size(); ++i){
      DrawContour(jets[i], cols[i][1]);
    }

  } else {
    // draw subjets
    for (int i=0; i<subjets.size(); ++i){
      for (int j=0; j<subjets[i].size(); ++j){  
        FixHistForContour(subjets[i][j]);
        DrawContour(subjets[i][j], cols[i][j]);
      }
    }
  }
  
  // draw all stable particles
  DrawStableParts(file, dir);

  // draw top's decay partons
  //DrawTopDecay(file);

  DrawTopDecay_withb(file, dir);
  DrawAntiTopDecay_withb(file, dir);
  DrawTops(file, dir);

  HOTVR_label_jetdisplay();

  Eventnumber_jetdisplay(num, "HOTVR");

  gPad->RedrawAxis();

  TString plot_dir = TString("plots/") + prefix;
  TString cmd = TString("mkdir -p ") + plot_dir; 
  gSystem->Exec(cmd.Data()); 

  TString outname = plot_dir + "/" + dir + "_hotvr"; 
  PrintInfo(file, outname, dir);

  //c->SaveAs("plots/" + prefix + dir + "_hotvr.pdf");
  c->SaveAs(outname + ".pdf");

}
  
