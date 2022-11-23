#include "style.h"

void SetupGlobalStyle()
{
  // general appearance and style

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(2);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(2);

  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.05, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(700); 
  gStyle->SetCanvasDefW(800); 

  // Margins:
  gStyle->SetPadTopMargin(    0.05  );
  gStyle->SetPadBottomMargin( 0.16  );
  gStyle->SetPadRightMargin(  0.05  );
  gStyle->SetPadLeftMargin(   0.19  );

  gStyle->UseCurrentStyle();

}

TCanvas* SetupCanvas(TString name)
{
  // optimised plots for including in theses or publications and documents

  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 500;
  CanHeight = 400;

  // set up the canvas
  TCanvas* can = new TCanvas(name,"Control Plots", CanWidth, CanHeight);
      
  can->SetTopMargin(    0.05  ); 
  can->SetBottomMargin( 0.16  );  
  can->SetRightMargin(  0.05  );
  can->SetLeftMargin(   0.16  ); 

  return can;
}

TCanvas* SetupCanvas2d(TString name)
{
  // optimised plots for including in theses or publications and documents

  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 500;
  CanHeight = 500;

  // set up the canvas
  TCanvas* can = new TCanvas(name,"Control Plots", CanWidth, CanHeight);
      
  can->SetTopMargin(    0.05  ); 
  can->SetBottomMargin( 0.14  );  
  can->SetRightMargin(  0.05  );
  can->SetLeftMargin(   0.14  ); 

  return can;
}

void Cosmetics(TH1* hist)
{

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetNdivisions(510);    
    
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.04);
    hist->GetYaxis()->SetLabelOffset(0.011);
    hist->GetYaxis()->SetNdivisions(505);    
  
}

void Cosmetics2d(TH2* hist)
{

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetNdivisions(510);    
    
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTickLength(0.03);
    hist->GetYaxis()->SetLabelOffset(0.011);
    hist->GetYaxis()->SetNdivisions(510);    
  
}

double NormaliseToArea(TH1* hist)
{
  // noramlise the histogram according to its area
  // the resulting histogram has area = 1

  double area = 0.;
  int ibeg = 1; // integrate from bin 1 
  int iend = 1;

  iend = hist->GetNbinsX();
  area = hist->Integral(ibeg,iend);
  if (area>0){
   	hist->Scale(1./area);
  }

  return area;
}

TLegend* SetupLegend(double shiftleft)
{
  Double_t xl1=.21+shiftleft, yl1=0.65, xl2=xl1+0.3, yl2=yl1+.15;
  TLegend* leg = new TLegend(xl1,yl1,xl2,yl2, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.055);
  return leg;
}

TLegend* SetupLegend2(double shiftleft)
{
  Double_t xl1=.21+shiftleft, yl1=0.62, xl2=xl1+0.3, yl2=yl1+.19;
  TLegend* leg = new TLegend(xl1,yl1,xl2,yl2, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->SetEntrySeparation(0.15);
  return leg;
}

TLegend* SetupROCLegend()
{
  //0.58
  Double_t xl1=.58, yl1=0.21, xl2=xl1+0.4, yl2=yl1+.3;
  TLegend* leg = new TLegend(xl1,yl1,xl2,yl2, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  return leg;
}

TLegend* SetupROCLegendSD()
{
  //0.58
  Double_t xl1=.50, yl1=0.21, xl2=xl1+0.4, yl2=yl1+.3;
  TLegend* leg = new TLegend(xl1,yl1,xl2,yl2, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  return leg;
}


TLegend* SetupEffLegend()
{
  Double_t xl1=.21, yl1=0.75, xl2=xl1+0.55, yl2=yl1+.16;
  TLegend* leg = new TLegend(xl1,yl1,xl2,yl2, NULL,"brNDC");
  leg-> SetNColumns(2);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  return leg;
}

void DrawInfo(TString info)
{
  TLatex *text1 = new TLatex(3.5, 24, info);
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.22);
  text1->SetTextFont(42);
  text1->SetTextSize(0.053);
  text1->SetY(0.892);
  text1->Draw();
}

void HOTVR_label()
{
  Int_t c_LGray = TColor::GetColor( "#e6e6e6" );

  TLatex *text1 = new TLatex(3.5, 24, "HOTVR");
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.72);
  text1->SetTextFont(62);
  text1->SetTextSize(0.075);
  text1->SetY(0.892);
  text1->SetTextColor(c_LGray);
  text1->Draw();
}

void HOTVR_label_jetdisplay(double left)
{
  Int_t c_LGray = TColor::GetColor( "#e6e6e6" );

  TLatex *text1 = new TLatex(3.5, 24, "HOTVR");
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.69+left);
  text1->SetTextFont(62);
  text1->SetTextSize(0.075);
  text1->SetY(0.91);
  text1->SetTextColor(c_LGray);
  text1->Draw();
}

void Eventnumber_jetdisplay(int num, TString algo)
{
  if (num < 1) return;

  Int_t c_vdGray = TColor::GetColor( "#404040" );
  TString st = TString::Format("Event %d", num);
  TLatex *text1 = new TLatex(3.5, 24, st.Data());
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.18);
  text1->SetTextFont(62);
  text1->SetTextSize(0.05);
  text1->SetY(0.91);
  text1->SetTextColor(c_vdGray);
  text1->Draw();

  if (algo!=""){
    Int_t c_vdGray = TColor::GetColor( "#404040" );
   TString st = TString::Format("clustered with %s", algo.Data());
   TLatex *text1 = new TLatex(3.5, 24, st.Data());
   text1->SetNDC();
   text1->SetTextAlign(13);
   text1->SetX(0.18);
   text1->SetTextFont(62);
   text1->SetTextSize(0.035);
   text1->SetY(0.86);
   text1->SetTextColor(c_vdGray);
   text1->Draw();    
  }

}

void FlatSpectrum(TH2F* h)
{
  // get sums and maximum
  double sums[1000];
  double maxs[1000];
  double tot = 0;
  
  for (int i=1; i<h->GetNbinsX(); ++i)
  { 
    maxs[i-1] = 0.;
    for (int j=1; j<h->GetNbinsY(); ++j)
    {
      sums[i-1] += h->GetBinContent(i,j);
      tot += h->GetBinContent(i,j);
      if (h->GetBinContent(i,j)>maxs[i-1]) maxs[i-1] = h->GetBinContent(i,j);
    }
  }

  // rescale
  for (int i=1; i<h->GetNbinsX(); ++i)
  {   
    //double fac = 1000.*sums[i-1]/tot;
    double fac = maxs[i-1]/tot;    
    if (!(fac>0)) continue;

    for (int j=1; j<h->GetNbinsY(); ++j)
    {
      h->SetBinContent(i,j, h->GetBinContent(i,j)/fac);
    }
  }

}

