#include "jet_display_helpers.h"

std::vector< std::vector<TH2F*> > GetSubjets(TFile* file, TString gdir)
{

  std::vector< std::vector<TH2F*> > subjets;

  file->cd();
  gDirectory->cd(gdir);

  TKey *key;
  TIter nextkey( gDirectory->GetListOfKeys() );
  while ( (key = (TKey*)nextkey())) {

      TObject *obj = key->ReadObj();
      if ( obj->IsA()->InheritsFrom( TH2::Class() ) ) {
        TH2F* h = (TH2F*) obj;
        TString name = h->GetName();
        //cout << "name = " << name << endl;

        // now get the right histograms depending on their names
        if (name.Contains("JetDisplay")){
          TObjArray* arr = name.Tokenize("_");

          //for (Int_t i=0; i<arr->GetEntries(); ++i){
          //  TString s = ((TObjString*) arr->At(i))->GetString();
          //  cout << "i = " << i << " s = " << s << endl;
          //}
          if (arr->GetEntries()<3) continue;
          TString jstr = ((TObjString*) arr->At(1))->GetString();
          //cout << "jet string = " << jstr << endl;
          if (!jstr.Contains("jet")) continue;
          jstr.ReplaceAll("jet", "");
          Int_t ijet = jstr.Atoi();
          if (ijet>9) break;

          // additional check: reject subjetpT histogram
          TString infostr = ((TObjString*) arr->At(2))->GetString();
          if (infostr == "subjetpT") continue;
          if (infostr == "subjetmass") continue;

          if (subjets.size()<ijet+1){ 
            std::vector<TH2F*> vec;
            subjets.push_back(vec);
          }
          if (CheckContent(h)){ 
            //RemovePartons(h, file);
            subjets[ijet].push_back(h);
          }
        }
      } else {
        obj->Delete();
      }
  }

  return subjets;

}

std::vector<TH2F*> GetJets(TFile* file, TString gdir)
{

  std::vector<TH2F*> jets;

  file->cd();
  gDirectory->cd(gdir);

  TKey *key;
  TIter nextkey( gDirectory->GetListOfKeys() );
  while ( (key = (TKey*)nextkey())) {

      TObject *obj = key->ReadObj();
      if ( obj->IsA()->InheritsFrom( TH2::Class() ) ) {
        TH2F* h = (TH2F*) obj;
        TString name = h->GetName();

        // now get the right histograms depending on their names
        if (name.Contains("JetDisplay")){
          TObjArray* arr = name.Tokenize("_");

          //for (Int_t i=0; i<arr->GetEntries(); ++i){
          //  TString s = ((TObjString*) arr->At(i))->GetString();
          //  cout << "i = " << i << " s = " << s << endl;
          //}
          if (arr->GetEntries()<2) continue;
          TString jstr = ((TObjString*) arr->At(1))->GetString();
          //cout << "jet string = " << jstr << endl;
          if (!jstr.Contains("pf")) continue;
          if (name.Contains("all")) continue;
          jstr.ReplaceAll("pf", "");
          Int_t ijet = jstr.Atoi();

          if (CheckContent(h)){ 
            //RemovePartons(h, file);
            jets.push_back(h);
          }
        }
      } else {
        obj->Delete();
      }
  }

  return jets;

}

TGraph* GetParticles(TFile* file, TString gdir)
{
  // always use the particles from the HOTVR clustering for display
  TString dirname = gdir;
  dirname.ReplaceAll("_akt", "");
  dirname.ReplaceAll("_ca", "");

  TH2F* pf = (TH2F*) file->Get(dirname + "/JetDisplay_pf_all");
  TGraph* parts = new TGraph();
  int k = 0;
  for (int i=1; i<pf->GetNbinsX()+1; ++i){
    for (int j=1; j<pf->GetNbinsY()+1; ++j){
      if (pf->GetBinContent(i,j)>1e-3){
        // smear the particle's position instead of using the bin center
        double x = pf->GetXaxis()->GetBinCenter(i);
        double y = pf->GetYaxis()->GetBinCenter(j);
        //Smear(x, pf->GetXaxis()->GetBinWidth(i));
        //Smear(y, pf->GetYaxis()->GetBinWidth(i));
        parts->SetPoint(k++, x, y);
      }  
    }
  } 
  return parts;
}

void NewEvent(UInt_t ievent)
{
  is_new_event = true;
  iseed = ievent;
}

UInt_t GetEventNumber(TString dir)
{
  UInt_t event = 0;
  TObjArray* arr = dir.Tokenize("_");
  if (arr->GetEntries()<2) return event;
  TString jstr = ((TObjString*) arr->At(arr->GetLast()))->GetString();
  jstr.ReplaceAll("event", "");
  event = jstr.Atoi();
  return event;
}

void Smear(double &x, double dx)
{
  if (!gRandom){
    gRandom = new TRandom3(iseed);
  }
  if (is_new_event){
    delete gRandom;
    gRandom = new TRandom3(iseed);
    is_new_event = false;  
  }
  double f = gRandom->Rndm();

  x += dx*(1.-2.*f);

}

void DrawTops(TFile* file, TString dir)
{
  TGraph* top = GetTopDecay_specific(file, "top", dir);
  top->SetMarkerStyle(29);
  top->SetMarkerSize(1.4);
  top->SetMarkerColor(kBlack);
  top->SetLineWidth(2);
  top->DrawClone("Psame"); 

  TGraph* atop = GetTopDecay_specific(file, "antitop", dir);
  atop->SetMarkerStyle(29);
  atop->SetMarkerSize(1.4);
  atop->SetMarkerColor(kBlack);
  atop->SetLineWidth(2);
  atop->DrawClone("Psame"); 

  //tops->SetMarkerSize(1.1);   
  //tops->DrawClone("Psame");
  //tops->SetMarkerSize(1.0);   
  //tops->DrawClone("Psame");      
  //tops->SetMarkerSize(0.9);   
  //tops->DrawClone("Psame");        
}

TGraph* GetTopDecay(TFile* file, TString gdir)
{
  TH2F* decay = (TH2F*) file->Get(gdir + "/JetDisplay_decay");  
  TGraph* g = new TGraph();
  int k = 0;
  for (int i=1; i<decay->GetNbinsX()+1; ++i){
    for (int j=1; j<decay->GetNbinsY()+1; ++j){
      if (decay->GetBinContent(i,j)>0){
        // smear the particle's position instead of using the bin center
        double x = decay->GetXaxis()->GetBinCenter(i);
        double y = decay->GetYaxis()->GetBinCenter(j);
        //Smear(x, decay->GetXaxis()->GetBinWidth(i));
        //Smear(y, decay->GetYaxis()->GetBinWidth(i));
        g->SetPoint(k++, x, y);
      }  
    }
  } 
  return g;
}
TGraph* GetTopDecay_specific(TFile* file, TString decayname, TString gdir)
{
  TH2F* decay = (TH2F*) file->Get(gdir + "/JetDisplay_" + decayname);  
  TGraph* g = new TGraph();
  int k = 0;
  for (int i=1; i<decay->GetNbinsX()+1; ++i){
    for (int j=1; j<decay->GetNbinsY()+1; ++j){
      if (decay->GetBinContent(i,j)>0){
        // smear the particle's position instead of using the bin center
        double x = decay->GetXaxis()->GetBinCenter(i);
        double y = decay->GetYaxis()->GetBinCenter(j);
        //Smear(x, decay->GetXaxis()->GetBinWidth(i));
        //Smear(y, decay->GetYaxis()->GetBinWidth(i));
        g->SetPoint(k++, x, y);
      }  
    }
  } 
  return g;
}
void DrawTopDecay(TFile* file, TString dir)
{
  TGraph* tdecay = GetTopDecay(file, dir);
  tdecay->SetMarkerStyle(24);
  tdecay->SetMarkerSize(1.2);
  tdecay->SetMarkerColor(kRed+2);
  tdecay->SetLineWidth(2);
  tdecay->DrawClone("Psame"); 
  tdecay->SetMarkerSize(1.1);   
  tdecay->DrawClone("Psame");
  tdecay->SetMarkerSize(1.0);   
  tdecay->DrawClone("Psame");      
  tdecay->SetMarkerSize(0.9);   
  tdecay->DrawClone("Psame");        
}

void DrawTopDecay_withb(TFile* file, TString dir)
{
  TGraph* wdecay = GetTopDecay_specific(file, "W_top", dir);
  TGraph* bdecay = GetTopDecay_specific(file, "b_top", dir);
  wdecay->SetMarkerStyle(24);
  wdecay->SetMarkerSize(1.2);
  wdecay->SetMarkerColor(kRed+2);
  wdecay->SetLineWidth(2);
  wdecay->DrawClone("Psame"); 
  wdecay->SetMarkerSize(1.1);   
  wdecay->DrawClone("Psame");
  wdecay->SetMarkerSize(1.0);   
  wdecay->DrawClone("Psame");      
  wdecay->SetMarkerSize(0.9);   
  wdecay->DrawClone("Psame");

  bdecay->SetMarkerStyle(32);
  bdecay->SetMarkerSize(1.2);
  bdecay->SetMarkerColor(kRed+2);
  bdecay->SetLineWidth(2);
  bdecay->DrawClone("Psame"); 
  bdecay->SetMarkerSize(1.1);   
  bdecay->DrawClone("Psame");
  bdecay->SetMarkerSize(1.0);   
  bdecay->DrawClone("Psame");      
  bdecay->SetMarkerSize(0.9);   
  bdecay->DrawClone("Psame");

}

void DrawAntiTopDecay_withb(TFile* file, TString dir)
{
  TGraph* wdecay = GetTopDecay_specific(file, "W_antitop", dir);
  TGraph* bdecay = GetTopDecay_specific(file, "b_antitop", dir);
  wdecay->SetMarkerStyle(24);
  wdecay->SetMarkerSize(1.2);
  wdecay->SetMarkerColor(kBlue+2);
  wdecay->SetLineWidth(2);
  wdecay->DrawClone("Psame"); 
  wdecay->SetMarkerSize(1.1);   
  wdecay->DrawClone("Psame");
  wdecay->SetMarkerSize(1.0);   
  wdecay->DrawClone("Psame");      
  wdecay->SetMarkerSize(0.9);   
  wdecay->DrawClone("Psame");        

  bdecay->SetMarkerStyle(32);
  bdecay->SetMarkerSize(1.2);
  bdecay->SetMarkerColor(kBlue+2);
  bdecay->SetLineWidth(2);
  bdecay->DrawClone("Psame"); 
  bdecay->SetMarkerSize(1.1);   
  bdecay->DrawClone("Psame");
  bdecay->SetMarkerSize(1.0);   
  bdecay->DrawClone("Psame");      
  bdecay->SetMarkerSize(0.9);   
  bdecay->DrawClone("Psame");        

}

void DrawStableParts(TFile* file, TString dir)
{
  TGraph* parts = GetParticles(file, dir);
  parts->SetMarkerStyle(20);
  parts->SetMarkerSize(0.4);
  parts->SetMarkerColor(TColor::GetColor( "#666666" ));
  parts->Draw("Psame");
} 

TH2F* GetClustersRejected(TFile* file, TString gdir)
{
  TH2F* h = (TH2F*) file->Get(gdir + "/JetDisplay_rejected_jets");
  for (int ib=1; ib<h->GetNbinsX()+1; ++ib){
    for (int jb=1; jb<h->GetNbinsY()+1; ++jb){
      if (h->GetBinContent(ib,jb)>0) h->SetBinContent(ib,jb,1);
    }
  }
  return h;  
}

// jets without subjets
TH2F* GetJetsWithoutSubjets(TFile* file, TString gdir)
{
  //TH2F* h = (TH2F*) file->Get(gdir + "/JetDisplay_radiation");
  TH2F* h = (TH2F*) file->Get(gdir + "/JetDisplay_rejected_jets");
  for (int ib=1; ib<h->GetNbinsX()+1; ++ib){
    for (int jb=1; jb<h->GetNbinsY()+1; ++jb){
      if (h->GetBinContent(ib,jb)>0) h->SetBinContent(ib,jb,1);
    }
  }
  return h;  
}

// jets and subjets rejected by pT cuts
TH2F* GetJetsWithLowPt(TFile* file, TString gdir)
{
  //TH2F* h = (TH2F*) file->Get(gdir + "/JetDisplay_radiation");
  TH2F* h = (TH2F*) file->Get(gdir + "/JetDisplay_rejected_subjets");
  for (int ib=1; ib<h->GetNbinsX()+1; ++ib){
    for (int jb=1; jb<h->GetNbinsY()+1; ++jb){
      if (h->GetBinContent(ib,jb)>0) h->SetBinContent(ib,jb,1);
    }
  }
  return h;  
}

// print out some information about the event
void PrintInfo(TFile* file, TString fname, TString gdir)
{
  //TH2F* h = (TH2F*) file->Get(gdir + "/JetDisplay_radiation");
  // bug in filling: rejected subjets are called "beam" here  
  TH1F* hpt = (TH1F*) file->Get(gdir + "/JetpT");  
  TH1F* hmass = (TH1F*) file->Get(gdir + "/Jetmass");
  TH1F* hphi = (TH1F*) file->Get(gdir + "/Jetphi");
  TH1F* heta = (TH1F*) file->Get(gdir + "/Jeteta");

  std::ofstream o; 
  o.open(fname + ".txt"); 

  TH2F* tdecay = (TH2F*) file->Get(gdir + "/JetDisplay_" + "top");  
  TH2F* atdecay = (TH2F*) file->Get(gdir + "/JetDisplay_" + "antitop");  
  for (int i=1; i<tdecay->GetNbinsX()+1; ++i){
    for (int j=1; j<tdecay->GetNbinsY()+1; ++j){
      if (tdecay->GetBinContent(i,j)>0){
        o << "top pt = " << tdecay->GetBinContent(i,j) << std::endl;
        o << "phi = " << tdecay->GetXaxis()->GetBinCenter(i) << ", eta = " 
          << tdecay->GetYaxis()->GetBinCenter(j) << std::endl;        
      }
      if (atdecay->GetBinContent(i,j)>0){
        o << "antitop pt = " << atdecay->GetBinContent(i,j) << std::endl;
        o << "phi = " << atdecay->GetXaxis()->GetBinCenter(i) << ", eta = " 
          << atdecay->GetYaxis()->GetBinCenter(j) << std::endl;
      }
    }
  }

  //double PI = 3.1415926535; 

  for (int ib=1; ib<hpt->GetNbinsX()+1; ++ib){  
    double pt = hpt->GetBinContent(ib); 
    double phi = hphi->GetBinContent(ib); 
    if (phi>M_PI) phi -= 2*M_PI;
    if (pt>0){
      std::cout << "Jet " << ib << ": " << std::endl;
      std::cout << "pt = " << hpt->GetBinContent(ib);
      std::cout << " m = "  << hmass->GetBinContent(ib);    
      std::cout << " phi = " << phi;
      std::cout << " eta = " << heta->GetBinContent(ib) << std::endl;

      o << "Jet " << ib << ": " << std::endl;
      o << "pt  = " << hpt->GetBinContent(ib);
      o << " m   = " << hmass->GetBinContent(ib);
      o << " phi = " << phi;
      o << " eta = " << heta->GetBinContent(ib) << std::endl;
    }
  }
  o << "\n"; 
  o << std::endl;
  o.close();

}

// subjets removed by softdrop criterion
TH2F* GetSubjetsRejectedBySD(TFile* file, TString gdir)
{
  TH2F* h = (TH2F*) file->Get(gdir + "/JetDisplay_radiation");
  if (!h) return NULL;
  for (int ib=1; ib<h->GetNbinsX()+1; ++ib){
    for (int jb=1; jb<h->GetNbinsY()+1; ++jb){
      if (h->GetBinContent(ib,jb)>0) h->SetBinContent(ib,jb,1);
    }
  }
  return h;  
}

void RemovePartons(TH2* jh, TFile* file, TString gdir)
{
  // remove partons first (no idea why they are filled)  
  TH2F* decay = (TH2F*) file->Get(gdir + "/JetDisplay_decay");  
  for (int i=1; i<decay->GetNbinsX()+1; ++i){
    for (int j=1; j<decay->GetNbinsY()+1; ++j){
      if (decay->GetBinContent(i,j)>0){ 
        // check surrounding bins
        int k = 0;
        if (jh->GetBinContent(i+1,j)) ++k;
        if (jh->GetBinContent(i-1,j)) ++k;
        if (jh->GetBinContent(i,j+1)) ++k;
        if (jh->GetBinContent(i,j-1)) ++k;
        if (k==0) jh->SetBinContent(i,j,0);
      }
    }
  }

}

bool CheckContent(TH2F* h)
{
  double max = h->GetMaximum();
  int entries = h->GetEntries();
  if (max>0 || entries>0) return true;
  return false;
}

void CreateContours(std::vector<TH2F*> jets)
{
    for (int j=0; j<jets.size(); ++j){         
      TH2F* h = jets[j];
      // loop over bins, create one area with constant height
      for (int ib=1; ib<h->GetNbinsX()+1; ++ib){
        for (int jb=1; jb<h->GetNbinsY()+1; ++jb){
          if (h->GetBinContent(ib,jb)>0) h->SetBinContent(ib,jb,1);
        }
      }
    }
}

void CreateContour(TH2F* h)
{
  // loop over bins, create one area with constant height
  for (int ib=1; ib<h->GetNbinsX()+1; ++ib){
    for (int jb=1; jb<h->GetNbinsY()+1; ++jb){
      if (h->GetBinContent(ib,jb)>0) h->SetBinContent(ib,jb,1);
    }
  }
}

void DrawContour(TH2F* hist, int col)
{

  // remove some trailing islands (bugs?)
  //hist->SetBinContent( 1,16,0);    
  //hist->SetBinContent( 8, 9,0);  
  //hist->SetBinContent(12,11,0);
  //hist->SetBinContent(42, 8,0);
  //hist->SetBinContent(50,12,0);

  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  //float cl;
  //int   fcol = kBlue;

  h->SetContourLevel( 0, 0.5 );

  h->SetLineWidth( 1 );
  h->SetFillColor( col );  
  h->SetLineColor( col );
  //h->Draw( "samecont3" );
  h->Draw( "samecont0" );
  
}

std::vector< std::vector<int> > GetColours()
{
  // orange hues
  std::vector<int> col1;
  col1.push_back(TColor::GetColor( "#ff944d" ));
  col1.push_back(TColor::GetColor( "#ffb380" ));
  col1.push_back(TColor::GetColor( "#ffd1b3" ));
  col1.push_back(TColor::GetColor( "#ffe0cc" ));
  col1.push_back(TColor::GetColor( "#fff0e6" ));
  for (int i=5; i<10; ++i) col1.push_back(TColor::GetColor( "#fff0e6" ));
 
  // green hues
  std::vector<int> col2;
  col2.push_back(TColor::GetColor( "#79d279" ));
  col2.push_back(TColor::GetColor( "#9fdf9f" ));
  col2.push_back(TColor::GetColor( "#c6ecc6" ));
  col2.push_back(TColor::GetColor( "#d9f2d9" ));
  col2.push_back(TColor::GetColor( "#ecf9ec" ));
  for (int i=5; i<10; ++i) col2.push_back(TColor::GetColor( "#ecf9ec" ));

  // blue hues
  std::vector<int> col3;
  col3.push_back(TColor::GetColor( "#4d94ff" ));
  col3.push_back(TColor::GetColor( "#80b3ff" ));
  col3.push_back(TColor::GetColor( "#b3d1ff" ));
  col3.push_back(TColor::GetColor( "#cce0ff" ));
  col3.push_back(TColor::GetColor( "#e6f0ff" ));
  for (int i=5; i<10; ++i) col3.push_back(TColor::GetColor( "#e6f0ff" ));

  // pink hues
  std::vector<int> col4;
  col4.push_back(TColor::GetColor( "#e600e6" ));
  col4.push_back(TColor::GetColor( "#ff1aff" ));
  col4.push_back(TColor::GetColor( "#ff99ff" ));
  col4.push_back(TColor::GetColor( "#ffccff" ));
  col4.push_back(TColor::GetColor( "#ffe6ff" ));
  for (int i=5; i<10; ++i) col4.push_back(TColor::GetColor( "#ffe6ff" ));

  std::vector< std::vector<int> > cols;
  cols.push_back(col1); // orange
  cols.push_back(col3); // blue
  cols.push_back(col2); // green
  //cols.push_back(col4); // pink
  for (int i=3; i<10; ++i) cols.push_back(col2);

  return cols;

}


void FixHistForContour(TH2F* h)
{

    // fix plot if N-sigma outside of scan range
    // scan along x-axis first
    double df = 0.49;

    for (int i=1; i<h->GetNbinsX()+1; ++i)
    {     
       Double_t v = h->GetBinContent(i,1);
       if (v >= 1.){
         h->SetBinContent(i, 1, df);
       }
       v = h->GetBinContent(i,h->GetNbinsY());
       if (v >= 1.){
         h->SetBinContent(i, h->GetNbinsY(), df);
       }
    }
    // now the y-axis
    for (int i=1; i<h->GetNbinsY()+1; ++i)
    {     
       Double_t v = h->GetBinContent(1,i);
       if (v >= 1.){
         h->SetBinContent(1, i, df);
       }
       v = h->GetBinContent(h->GetNbinsX(), i);
       if (v >= 1.){
         h->SetBinContent(h->GetNbinsX(), i, df);
         //cout << "fixing bin content for v = " << v << " in bin " << i << endl;
       }
    }    

}



