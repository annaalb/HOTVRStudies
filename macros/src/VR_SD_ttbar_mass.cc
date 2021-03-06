#include "../include/CentralInclude.h"


int main(int argc, char* argv[]){
  SetStyle();
  TString pt = argv[1];
  TString pt_max = argv[2];

  TFile *input1 = new TFile(" /nfs/dust/cms/user/albrecha/uhh2_102X_v2/HOTVRStudiesOutput/root/VR_SD/uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23_VR_SD_24_02_z01_b01.root", "READ");
  TFile *input2 = new TFile(" /nfs/dust/cms/user/albrecha/uhh2_102X_v2/HOTVRStudiesOutput/root/VR_SD/uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23_VR_SD_27_02_z001_b01.root", "READ");

  //only HOTVR
    TH1F* hotvr = (TH1F*)input2->Get("HOTVRStudiesHists/mass_jets_" + pt);
    TH1F* vr = (TH1F*)input2->Get("HOTVRStudiesHists/mass_jets_noMJ_" + pt);
  //beta=1
    TH1F* b1z01 = (TH1F*)input1->Get("HOTVRStudiesHists/mass_jets_SD_" +pt);
    TH1F* b1z001 = (TH1F*)input2->Get("HOTVRStudiesHists/mass_jets_SD_" + pt);

    TCanvas* canvas = new TCanvas("c", "c", 600, 600);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.2);

    hotvr->SetTitle("");

  //only HOTVR
    hotvr->SetLineColor(kBlack);
    hotvr->SetLineWidth(2);
    hotvr->Scale(1/hotvr->Integral());
    hotvr->Draw("H");
    TF1* fit = new TF1("fit", "gaus", 160, 190);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    fit->SetLineStyle(kDotted);
    hotvr->Fit("fit", "R");

    vr->SetLineColor(kRed);
    vr->SetLineWidth(2);
    vr->Scale(1/vr->Integral());
    vr->Draw("H same");
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    fit->SetLineStyle(kDotted);
    vr->Fit("fit", "R");

  //beta 1
      b1z01->SetLineColor(kAzure+7);
      b1z01->SetLineWidth(2);
      b1z01->Scale(1/b1z01->Integral());
      b1z01->Draw("H same");
      fit->SetLineColor(kAzure+7);
      fit->SetLineWidth(2);
      fit->SetLineStyle(kDotted);
      b1z01->Fit("fit", "R");

      b1z001->SetLineColor(kGreen+2);
      b1z001->SetLineWidth(2);
      b1z001->Scale(1/b1z001->Integral());
      b1z001->Draw("H same");
      fit->SetLineColor(kGreen+2);
      fit->SetLineWidth(2);
      fit->SetLineStyle(kDotted);
      b1z001->Fit("fit", "R");


    hotvr->GetYaxis()->SetRangeUser(0,0.14);
    hotvr->GetXaxis()->SetLabelSize(0.04);
    hotvr->GetXaxis()->SetTitleSize(0.06);
    hotvr->GetXaxis()->SetRangeUser(150, 200);
    hotvr->GetYaxis()->SetLabelSize(0.04);
    hotvr->GetYaxis()->SetTitleSize(0.06);

    hotvr->GetXaxis()->SetTitle("m_{jet} [GeV]");
    hotvr->GetYaxis()->SetTitle("fraction of jets");

    auto legend = new TLegend(0.25, 0.6, 0.85, 0.85);
    legend->SetTextSize(.04);
    legend->SetHeader("t#bar{t}, "+ pt +" < p_{T} [GeV] < " + pt_max);
    legend->AddEntry(hotvr,"HOTVR", "l");
    legend->AddEntry(vr,"VR", "l");
    legend->AddEntry(b1z01,"VR+SD, z_{cut} = 0.1, #beta = 1", "l");
    legend->AddEntry(b1z001,"VR+SD, z_{cut} = 0.01, #beta = 1", "l");


    legend->Draw();

    canvas->SaveAs("/nfs/dust/cms/user/albrecha/uhh2_102X_v2/HOTVRStudiesOutput/plots/VR_SD/ttbar_flat_VR_SD_mass_norm_top_"+ pt +".eps");

  return 0;
}
