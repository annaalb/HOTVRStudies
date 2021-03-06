#include "UHH2/HOTVRStudies/include/HOTVRJetsHists.h"

using namespace std;
using namespace uhh2;

HOTVRJetsHists::HOTVRJetsHists(Context & ctx, const string & dirname, bool is_qcd):
Hists(ctx, dirname),
b_is_qcd(is_qcd)
{
// book all histograms here
// processed events
processed_events_ttbar = book<TH1F>("processed_events_ttbar", "processed events ttbar", 1,0,1);
processed_events_qcd = book<TH1F>("processed_events_qcd", "processed events qcd", 1,0,1);
// general hists
hist_pt = book<TH1F>("p_{T}", "p_{T} [GeV]", 200, 0, 2000);
hist_pt_sum_subjets = book<TH1F>("sumsubjetsp_{T}", "sumsubjetsp_{T} [GeV]", 200, 0, 2000);

hist_mass = book<TH1F>("mass", "mass [GeV]", 100, 0, 300);
hist_eta = book<TH1F>("eta", "#eta", 100, -6, 6);
hist_phi = book<TH1F>("phi", "#phi", 40, -4, 4);
hist_energy = book<TH1F>("energy", "energy [GeV]", 100, 0, 2200);
// HOTVR specific parameters
hist_max_distance = book<TH1F>("max_distance", "max_distance [GeV]", 100, 0, 2);
hist_mmin = book<TH1F>("mmin", "mmin [GeV]", 100, 0, 200);
hist_fpt1 = book<TH1F>("fpt", "fpt", 20, 0, 1);
hist_nsubjets = book<TH1F>("nsubjets", "nsubjets", 20, -0.5, 19.5);
// N-subjettiness
hist_tau1 = book<TH1F>("tau1", "tau1", 100, 0, 1);
hist_tau2 = book<TH1F>("tau2", "tau2", 100, 0, 1);
hist_tau3 = book<TH1F>("tau3", "tau3", 100, 0, 1);
hist_tau21 = book<TH1F>("tau21", "tau21", 100, 0, 1);
hist_tau32 = book<TH1F>("tau32", "tau32", 100, 0, 1);

hist_matching_radius = book<TH1F>("matching_radius", "matching_radius", 100, 0, 2);
hist_max_distance_minus_matching_radius = book<TH1F>("max_distance_minus_matching_radius", "max_distance minus matching_radius", 100, 0, 2);

hist_pt_subjet1 = book<TH1F>("subjet1_p_{T}", "subjet1 p_{T} [GeV]", 200, 0, 1000);
hist_pt_subjet2 = book<TH1F>("subjet2_p_{T}", "subjet2 p_{T} [GeV]", 200, 0, 1000);
hist_pt_subjet3 = book<TH1F>("subjet3_p_{T}", "subjet3 p_{T} [GeV]", 200, 0, 1000);
hist_pt_subjet4 = book<TH1F>("subjet4_p_{T}", "subjet4 p_{T} [GeV]", 200, 0, 1000);
hist_pt_subjet5 = book<TH1F>("subjet5_p_{T}", "subjet5 p_{T} [GeV]", 200, 0, 1000);

hist_njets = book<TH1F>("njets", "njets", 20, -0.5, 19.5);
}

void HOTVRJetsHists::fill(const Event & event){}

void HOTVRJetsHists::fill_topjet(const Event & event, const TopJet & jet){
  if (b_is_qcd) {processed_events_qcd->Fill(0);}
  else{processed_events_ttbar->Fill(0);}

    hist_pt->Fill(jet.pt());
  //  hist_mass->Fill(jet.v4().M());
    LorentzVector subjet_sum;
    for (const auto s : jet.subjets()) {
    subjet_sum += s.v4();
    }
    hist_mass->Fill(subjet_sum.M());
    hist_pt_sum_subjets->Fill(subjet_sum.pt());

    hist_eta->Fill(jet.eta());
    hist_phi->Fill(jet.phi());
    hist_energy->Fill(jet.v4().E());
    hist_max_distance->Fill(jet.max_distance());
    hist_mmin->Fill(jet.hotvr_mmin());
    hist_fpt1->Fill(jet.hotvr_fpt1());
    hist_nsubjets->Fill(jet.subjets().size());

    hist_tau1->Fill(jet.tau1_groomed());
    hist_tau2->Fill(jet.tau2_groomed());
    hist_tau3->Fill(jet.tau3_groomed());
    hist_tau21->Fill(jet.tau2_groomed()/jet.tau1_groomed());
    hist_tau32->Fill(jet.tau3_groomed()/jet.tau2_groomed());

    double rho = 600;
    double pt = jet.pt();
    double matching_radius = rho/pt;
    if(matching_radius<0.1){matching_radius=0.1;};
    if(matching_radius>1.5){matching_radius=1.5;};

    hist_matching_radius->Fill(matching_radius);

    double x = abs(jet.max_distance()-matching_radius);
    hist_max_distance_minus_matching_radius->Fill(x);

    // fill pt of subjets
    std::vector<Jet> subjets = jet.subjets();
    if (subjets.size()>0) {hist_pt_subjet1->Fill(subjets[1].pt());}
    if (subjets.size()>1) {hist_pt_subjet2->Fill(subjets[2].pt());}
    if (subjets.size()>2) {hist_pt_subjet3->Fill(subjets[3].pt());}
    if (subjets.size()>3) {hist_pt_subjet4->Fill(subjets[4].pt());}
    if (subjets.size()>4) {hist_pt_subjet5->Fill(subjets[5].pt());}

}

void HOTVRJetsHists::fill_n_jets(const Event & event, const vector<TopJet> & jets){
  hist_njets->Fill(jets.size());
}



HOTVRJetsHists::~HOTVRJetsHists(){}
