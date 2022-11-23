#include "jet_display_akt.h"
#include "jet_display_helpers.h"
#include "style.h"
#include <vector>

void jet_display_akt(TFile* file, TString dir, TString prefix, Int_t num)
{ 
  SetupGlobalStyle(); 

  //TFile* file = new TFile("RocCycle.MC.TTbarEff_gen_jetdisplay.root", "READ");
  
  // store all jets in this vector
  std::vector<TH2F*> jets;
  jets = GetJets(file, dir);
  CreateContours(jets);

  TCanvas *c = SetupCanvas2d();

  std::vector< std::vector<int> > cols;
  cols = GetColours(); 

  // draw jets
  for (int i=0; i<jets.size(); ++i){
    jets[i]->GetXaxis()->SetTitle("#phi");
    jets[i]->GetYaxis()->SetTitle("#eta");
    jets[i]->SetTitle("");
    Cosmetics2d(jets[i]);
    FixHistForContour(jets[i]);
    DrawContour(jets[i], cols[i][0]);
    
  }

  // draw all stable particles
  DrawStableParts(file, dir);

  // draw top's decay partons
  DrawTopDecay(file, dir);

  Eventnumber_jetdisplay(num, "CA");

  gPad->RedrawAxis();

  c->SaveAs("plots/softdrop_" + prefix + dir + ".pdf");

}
  
