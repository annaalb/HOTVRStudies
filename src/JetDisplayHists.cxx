#include "UHH2/HOTVRStudies/include/JetDisplayHists.h"
#include <math.h>

using namespace std;
using namespace uhh2;
using namespace fastjet;
using namespace contrib;

//TODO add phishift
JetDisplayHists::JetDisplayHists(Context & ctx, const string & dirname):
Hists(ctx, dirname)
{
  // book all histograms here

  //TODO if needed use vector of hists here
  //jet display
  int bins=50;
  for(int o=0;o<10;o++) { // loop over 10 possible jets
    TString hname3 = TString::Format("JetDisplay_jet%i_subjetmass",o);
    book<TH2F>(hname3,hname3,10,0,9,300,0,300);
    hname3 = TString::Format("JetDisplay_jet%i_subjetpT",o);
    book<TH2F>(hname3,hname3,10,0,9,300,0,3000);

    for(int p=0;p<10;p++) { // loop over 10 possible subjets
      TString hname2 = TString::Format("JetDisplay_jet%i_subjet%i", o,p);
      book<TH2F>(hname2,"Jet event display",bins,-PI,PI,bins,-PI,PI);
    } // end loop over subjets
  } // end loop over 10 jets

  beam = book<TH2F>("JetDisplay_rejected_jets","HOTVR jets without subjets",bins,-PI,PI,bins,-PI,PI);
  radiation = book<TH2F>("JetDisplay_radiation","Soft and wide angle radiation, rejected by SD veto",bins,-PI,PI,bins,-PI,PI);
  hist_rejected_subjets = book<TH2F>("JetDisplay_rejected_subjets","Jets and subjets discarded because of pT cuts.",bins,-PI,PI,bins,-PI,PI);
  //book<TH2F>("JetDisplay","Jet event display",bins,-PI,PI,bins,-PI,PI);
  //book<TH2F>("JetDisplay1","Jet event display",bins,-PI,PI,bins,-PI,PI);
  fatjet = book<TH2F>("JetDisplay_fatjet","Jet event display",bins,-PI,PI,bins,-PI,PI);
  //book<TH2F>("JetDisplay_subjets","Jet event display",bins,-PI,PI,bins,-PI,PI);
  pf0 = book<TH2F>("JetDisplay_pf0","Jet event display",bins,-PI,PI,bins,-PI,PI);
  pf1 = book<TH2F>("JetDisplay_pf1","Jet event display",bins,-PI,PI,bins,-PI,PI);
  pf2 = book<TH2F>("JetDisplay_pf2","Jet event display",bins,-PI,PI,bins,-PI,PI);
  pf_all = book<TH2F>("JetDisplay_pf_all","Jet event display",bins,-PI,PI,bins,-PI,PI);
  //book<TH2F>("JetDisplay_gluon","Jet event display",bins,-PI,PI,bins,-PI,PI);
  //book<TH2F>("JetDisplay_light","Jet event display",bins,-PI,PI,bins,-PI,PI);
  decay = book<TH2F>("JetDisplay_decay","Jet event display",bins,-PI,PI,bins,-PI,PI);

  jetmass = book<TH1F>("Jetmass","m_{jet}",10,0,10);
  jetpt = book<TH1F>("JetpT","jet p_{T}",10,0,10);
  jetphi = book<TH1F>("Jetphi","jet #phi",10,0,10);
  jeteta = book<TH1F>("Jeteta","jet #eta",10,0,10);

  top = book<TH2F>("JetDisplay_top","Jet event display",bins,-PI,PI,bins,-PI,PI);
  b_top = book<TH2F>("JetDisplay_b_top","Jet event display",bins,-PI,PI,bins,-PI,PI);
  W_top = book<TH2F>("JetDisplay_W_top","Jet event display",bins,-PI,PI,bins,-PI,PI);
  antitop = book<TH2F>("JetDisplay_antitop","Jet event display",bins,-PI,PI,bins,-PI,PI);
  b_antitop = book<TH2F>("JetDisplay_b_antitop","Jet event display",bins,-PI,PI,bins,-PI,PI);
  W_antitop = book<TH2F>("JetDisplay_W_antitop","Jet event display",bins,-PI,PI,bins,-PI,PI);

  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
  h_parts=ctx.get_handle<std::vector<PseudoJet>>("parts");

  h_jets=ctx.get_handle<std::vector<PseudoJet>>("jets");
  h_rejected_cluster=ctx.get_handle<std::vector<PseudoJet>>("rejected_cluster");
  h_soft_cluster=ctx.get_handle<std::vector<PseudoJet>>("soft_cluster");
  h_rejected_subjets=ctx.get_handle<std::vector<PseudoJet>>("rejected_subjets");
}

void JetDisplayHists::fill(const Event & event){

  //---------get all objects from event handle----------
  cout << "Filling info for jet displays" << endl;
  const auto & ttbargen = event.get(h_ttbargen);
  const auto & parts = event.get(h_parts);

  // these are the jets without subjets, i.e. which have not fulfilled the soft drop criterion
  const auto & rejected_cluster = event.get(h_rejected_cluster);
  // soft cluster: these are the jets which got discarded by the soft drop algorithm (the softer one)
  const auto & soft_cluster = event.get(h_soft_cluster);
  // these are the subjets rejected by the pt_sub cut, not filled at the moment
  const auto & rejected_subjets = event.get(h_rejected_subjets);
  const auto & jets = event.get(h_jets);


  //-----get tops and decay products------------
  std::vector<GenParticle> Had_Tops;
  std::vector<GenParticle> Had_AntiTops;

  std::vector<GenParticle> decay_products;
  std::vector<GenParticle> decay_products_Antitop;

  std::vector<GenParticle> W_Top;
  std::vector<GenParticle> W_AntiTop;
  GenParticle b_Top;
  GenParticle b_AntiTop;

  if(ttbargen.IsTopHadronicDecay()){
    Had_Tops.push_back(ttbargen.Top());
    decay_products.push_back(ttbargen.bTop());
    decay_products.push_back(ttbargen.Wdecay1());
    decay_products.push_back(ttbargen.Wdecay2());
    b_Top=ttbargen.bTop();
    W_Top.push_back(ttbargen.Wdecay1());
    W_Top.push_back(ttbargen.Wdecay2());
  //  std::cout << " top "<<ttbargen.Top().pt() << '\n';

  }
  if(ttbargen.IsAntiTopHadronicDecay()){
    //Had_Tops.push_back(ttbargen.Antitop());
    decay_products.push_back(ttbargen.bAntitop());
    decay_products.push_back(ttbargen.WMinusdecay1());
    decay_products.push_back(ttbargen.WMinusdecay2());
    Had_AntiTops.push_back(ttbargen.Antitop());
    b_AntiTop=ttbargen.bAntitop();
    W_AntiTop.push_back(ttbargen.WMinusdecay1());
    W_AntiTop.push_back(ttbargen.WMinusdecay2());
  //  std::cout << "anti top "<< ttbargen.Antitop().pt() << '\n';

  }
  else if(!(ttbargen.IsTopHadronicDecay() || ttbargen.IsAntiTopHadronicDecay())){
    throw runtime_error("no hadronic Top found");
  }


  //-------fill hists-----------
  for(uint y=0;y<rejected_cluster.size();y++){
    for(uint k=0; k<rejected_cluster[y].constituents().size();k++){
      double phi = rejected_cluster[y].constituents().at(k).phi_std();
      double eta = rejected_cluster[y].constituents().at(k).eta();
      double pt = rejected_cluster[y].constituents().at(k).pt();
      if (pt<1.) pt = 1;
      beam->Fill(phi, eta, pt );
    }
  }

  for(uint y=0;y<soft_cluster.size();y++){
    for(uint k=0;k<soft_cluster[y].constituents().size();k++){
      double phi = soft_cluster[y].constituents().at(k).phi_std();
      double eta = soft_cluster[y].constituents().at(k).eta();
      double pt  = soft_cluster[y].constituents().at(k).pt();
      if (pt<1.) pt = 1;
      radiation->Fill(phi, eta, pt);
    }
  }

  for(uint y=0;y<rejected_subjets.size();y++){
    for(uint k=0;k<rejected_subjets[y].constituents().size();k++){
      double phi = rejected_subjets[y].constituents().at(k).phi_std();
      double eta = rejected_subjets[y].constituents().at(k).eta();
      double pt  = rejected_subjets[y].constituents().at(k).pt();
      if (pt<1.) pt = 1;
      hist_rejected_subjets->Fill(phi, eta, pt);
    }
  }

  uint Njets = jets.size();
  if (Njets>10){
    Njets = 10;
  }
  for(uint o=0;o<Njets;o++) { // loop over jets

    if(jets[o].has_user_info<HOTVRinfo>()) {
      std::vector<fastjet::PseudoJet> SortedSubJets=jets[o].user_info<HOTVRinfo>().subjets();
      fatjet->Fill(jets[o].phi_std(),jets[o].eta(),jets[o].pt());

      jetmass->SetBinContent(o+1,jets[o].m());
      jetpt->SetBinContent(o+1,jets[o].pt());
      jetphi->SetBinContent(o+1,jets[o].phi());
      jeteta->SetBinContent(o+1,jets[o].eta());

      if (SortedSubJets.size()>1) { // select only jets with more than one subjet
        uint Nsubjets = SortedSubJets.size();
        if (Nsubjets>10){
          Nsubjets = 10;
        }
        for(uint p=0;p<Nsubjets;p++) { // loop over subjets
          TString hname3 = TString::Format("JetDisplay_jet%i_subjetmass",o);
          hist(hname3)->Fill(p,SortedSubJets.at(p).m());
          hname3 = TString::Format("JetDisplay_jet%i_subjetpT",o);
          hist(hname3)->Fill(p,SortedSubJets.at(p).pt());
          TString hname2 = TString::Format("JetDisplay_jet%i_subjet%i", o,p);

          for(uint h=0;h<SortedSubJets.at(p).constituents().size();h++){
            double phi = SortedSubJets.at(p).constituents().at(h).phi_std();
            double eta = SortedSubJets.at(p).constituents().at(h).eta();
            double pt  = SortedSubJets.at(p).constituents().at(h).pt();
            if (pt<1.) pt = 1;
            ((TH2D*)hist(hname2))->Fill(phi, eta, pt);
          }
        } // end loop over subjets
      }
    }
    for(uint h=0;h<jets[o].constituents().size();h++) { // loop over jet constituents
      double phi = jets[o].constituents().at(h).phi_std();
      double eta = jets[o].constituents().at(h).eta();
      double pt  = jets[o].constituents().at(h).pt();
      if (pt<1.) pt = 1;
      if(o==0) pf0->Fill(phi, eta, pt);
      if(o==1) pf1->Fill(phi, eta, pt);
      if(o==2) pf2->Fill(phi, eta, pt);

    } // end loop over jet constituents
  } // end loop over jets

  // ------fill hists with genparticles (tops and decay products)-------
  for(uint j=0; j< Had_Tops.size();j++){
    top->Fill(Had_Tops[j].phi(),Had_Tops[j].eta(),Had_Tops[j].pt());
    for(uint o=0;o<decay_products.size();o++){
      decay->Fill(decay_products[o].phi(),decay_products[o].eta(),decay_products[o].pt());
    }
    for(uint o=0;o<W_Top.size();o++){
      W_top->Fill(W_Top[o].phi(),W_Top[o].eta(),W_Top[o].pt());
    }
    b_top->Fill(b_Top.phi(),b_Top.eta(),b_Top.pt());

  }
  // same for anti top
  for(uint j=0; j< Had_AntiTops.size();j++){
    antitop->Fill(Had_AntiTops[j].phi(),Had_AntiTops[j].eta(),Had_AntiTops[j].pt());
    for(uint o=0;o<decay_products.size();o++){
      decay->Fill(decay_products[o].phi(),decay_products[o].eta(),decay_products[o].pt());
    }
    for(uint o=0;o<W_Top.size();o++){
      W_antitop->Fill(W_AntiTop[o].phi(),W_AntiTop[o].eta(),W_AntiTop[o].pt());
    }
    b_antitop->Fill(b_AntiTop.phi(),b_AntiTop.eta(),b_AntiTop.pt());
  }

  for(uint o=0;o<parts.size();o++){
    double phi = parts[o].phi_std();
    double eta = parts[o].eta();
    double pt  = parts[o].pt();
    //if (pt<1.) pt = 1;
    pf_all->Fill(phi, eta, pt);
  }

} // end fill function


    JetDisplayHists::~JetDisplayHists(){}
