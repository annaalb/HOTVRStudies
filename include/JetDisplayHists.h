#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/TTbarGen.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/HOTVRinfo.hh"
#include "fastjet/contrib/HOTVR.hh"

#include <iostream>

#include "TH1F.h"
#include "TH2F.h"

/**  Class to fill general hists for the Jet Displays.
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class JetDisplayHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    JetDisplayHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;

    virtual ~JetDisplayHists();

protected:
    TH2F *beam, *radiation, *hist_rejected_subjets, *fatjet, *pf0, *pf1, *pf2, *pf_all, *top, *decay;
    TH1F *jetmass, *jetpt, *jetphi, *jeteta;
    TH2F *b_top, *W_top, *b_antitop, *W_antitop, *antitop;

    uhh2::Event::Handle<TTbarGen>h_ttbargen;
    uhh2::Event::Handle<std::vector<fastjet::PseudoJet>>h_parts;

    uhh2::Event::Handle<std::vector<fastjet::PseudoJet>>h_jets;
    uhh2::Event::Handle<std::vector<fastjet::PseudoJet>>h_rejected_cluster;
    uhh2::Event::Handle<std::vector<fastjet::PseudoJet>>h_soft_cluster;
    uhh2::Event::Handle<std::vector<fastjet::PseudoJet>>h_rejected_subjets;

    // uhh2::Event::Handle<std::vector<std::vector<fastjet::PseudoJet>>>h_jets;
    // uhh2::Event::Handle<std::vector<std::vector<fastjet::PseudoJet>>>h_rejected_cluster;
    // uhh2::Event::Handle<std::vector<std::vector<fastjet::PseudoJet>>>h_soft_cluster;
    // uhh2::Event::Handle<std::vector<std::vector<fastjet::PseudoJet>>>h_rejected_subjets;

    // uhh2::Event::Handle<std::vector<std::vector<fastjet::PseudoJet>>>h_jet0_subjets;
    // uhh2::Event::Handle<std::vector<std::vector<fastjet::PseudoJet>>>h_jet1_subjets;
    // uhh2::Event::Handle<std::vector<std::vector<fastjet::PseudoJet>>>h_jet2_subjets;
    // uhh2::Event::Handle<std::vector<std::vector<fastjet::PseudoJet>>>h_jet3_subjets;

};
