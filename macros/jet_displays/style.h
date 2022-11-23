#ifndef __style_h_
#define __style_h_

#include <TStyle.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

void SetupGlobalStyle();

TCanvas* SetupCanvas(TString name = "canvas");

TCanvas* SetupCanvas2d(TString name = "canvas");

void Cosmetics(TH1* hist);

void Cosmetics2d(TH2* hist);

double NormaliseToArea(TH1* hist);

TLegend* SetupLegend(double shiftleft=0.);

TLegend* SetupLegend2(double shiftleft=0.);

TLegend* SetupROCLegend();

TLegend* SetupROCLegendSD();

TLegend* SetupEffLegend();

void DrawInfo(TString info);

void HOTVR_label(); 

void HOTVR_label_jetdisplay(double left=0.);

void Eventnumber_jetdisplay(int num, TString algo = "");

void FlatSpectrum(TH2F* h); 

#endif //__style_h_
