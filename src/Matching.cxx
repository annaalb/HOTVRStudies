#include "UHH2/HOTVRStudies/include/Matching.h"

using namespace std;
using namespace fastjet;
using namespace contrib;

//constructor
Matching::Matching()
{

}
//converts a genparticle into a fastjet::pseudojet
PseudoJet Matching::convert_particle(GenParticle* genparticle)
{
  PseudoJet gen_particle(genparticle->v4().Px(),genparticle->v4().Py(),genparticle->v4().Pz(),genparticle->v4().E());
  return gen_particle;
}
/*
██████   ██████   ██████  ██
██   ██ ██    ██ ██    ██ ██
██████  ██    ██ ██    ██ ██
██   ██ ██    ██ ██    ██ ██
██████   ██████   ██████  ███████
*/
// keep all neutrinos
bool Matching::IsNeutrino(GenParticle* p)
{
  int id = abs(p-> pdgId());
  if (id==12 || id==14 || id==16) return true;
  else return false;
}
// keep all tops
bool Matching::IsTop(GenParticle* p)
{
  int id = p->pdgId();
  if (abs(id)==6) return true;
  else return false;
}
// keep all W bosons
bool Matching::IsW(GenParticle* p)
{
  int id = p->pdgId();
  if (abs(id)==24) return true;
  else return false;
}
// keep all b quarks
bool Matching::Isb(GenParticle* p)
{
  int id = p->pdgId();
  if (abs(id)==5) return true;
  else return false;
}
// keep all quarks, gluons, W, Z bosons
bool Matching::IsParton(GenParticle* p)
{
	int fid = abs(p->pdgId());
	bool isp = false;
	if (fid<7)   isp = true;
	if (fid==21) isp = true;
	if (fid==23) isp = true;
	if (fid==24) isp = true;
 	return isp;
}
//keep all particles on hadron level
bool Matching::IsStableParticle(GenParticle* p)
{
  int st = p->status();
  if (st==1) return true;
  else return false;
}
// check if this top quark decays into W+b
bool Matching::BeforeTopDecay(GenParticle* p, vector<GenParticle>* genparticles)
{
	// get daughters
	GenParticle* d1 = GetDaughter(p, genparticles, 1);
	GenParticle* d2 = GetDaughter(p, genparticles, 2);

	// daughters found?
	if (!d1) return false;
	if (!d2) return false;

	bool flag = false;
	// check ID of daughters
	int id1 = abs(d1->pdgId());
	int id2 = abs(d2->pdgId());

	if (id1==5  && id2==24) flag = true;
	if (id1==24 && id2==5)  flag = true;

	return flag;
}
//keep all partons that do not decay further or hadronize -> final state partons
bool Matching::FinalStateParton(GenParticle* p, vector<GenParticle>* genparticles)
{
  GenParticle* d1 = GetDaughter(p, genparticles, 1);
	GenParticle* d2 = GetDaughter(p, genparticles, 2);
	bool result = false;

  if(std::abs(p->pdgId())<100){//is a parton
    if (!d1 && !d2) result = true; //has it no or hadronic daughters? -> final state parton

	   if (d1){
	      int fid = abs(d1->pdgId());
	      if (fid>100) result = true;
	      }
	   if (d2){
	      int fid = abs(d2->pdgId());
	      if (fid>100) result = true;
	     }
  }

  // protect against beam remnant
	LorentzVector v = p->v4();
	if (abs(v.Eta()) > 8) result = false;
	if (isinf(v.Eta())) result = false;
	if (isnan(v.Eta())) result = false;

  return result;
}
//check if the top decays hadronically
bool Matching::IsHadronic(GenParticle* p, vector<GenParticle>* genparticles)
{
  GenParticle* d1 = GetDaughter(p, genparticles, 1);
  GenParticle* d2 = GetDaughter(p, genparticles, 2);
  GenParticle* d;
  if(abs(d1->pdgId())==24) d=d1;
  else d=d2;
  while(abs(d->pdgId())==24){ //check daughters of W -> are they hadrons?
    d1= GetDaughter(d, genparticles, 1);
    d2 = GetDaughter(d, genparticles, 2);
    if(abs(d1->pdgId())==24) d=d1;
    else d=d2;
  }
  if(abs(d->pdgId())<6) return true;
  else return false;
}
/*
 ██████  ███████ ████████     ██████   █████  ██    ██  ██████  ██   ██ ████████ ███████ ██████
██       ██         ██        ██   ██ ██   ██ ██    ██ ██       ██   ██    ██    ██      ██   ██
██   ███ █████      ██        ██   ██ ███████ ██    ██ ██   ███ ███████    ██    █████   ██████
██    ██ ██         ██        ██   ██ ██   ██ ██    ██ ██    ██ ██   ██    ██    ██      ██   ██
 ██████  ███████    ██        ██████  ██   ██  ██████   ██████  ██   ██    ██    ███████ ██   ██
*/

//returns the daughter 1 or 2 of a genparticle
GenParticle* Matching::GetDaughter(GenParticle* p, vector<GenParticle>* genparticles, int n)
{
  uint index = 0;
  GenParticle* d = NULL;
  if(n<=1){
    index = p->daughter1();
  } else {
    index = p->daughter2();
  }
  if (index>0 && index<genparticles->size()){
		d = &genparticles->at(index-1);  }
    return d;
}
/*
██    ██ ██████  ██████   █████  ████████ ███████ ███████ ██   ██ ██ ██████  ██      ██ ███████ ████████
██    ██ ██   ██ ██   ██ ██   ██    ██    ██      ██      ██  ██  ██ ██   ██ ██      ██ ██         ██
██    ██ ██████  ██   ██ ███████    ██    █████   ███████ █████   ██ ██████  ██      ██ ███████    ██
██    ██ ██      ██   ██ ██   ██    ██    ██           ██ ██  ██  ██ ██      ██      ██      ██    ██
 ██████  ██      ██████  ██   ██    ██    ███████ ███████ ██   ██ ██ ██      ███████ ██ ███████    ██
*/

// go through all decays from the top and flag those to be skipped
// particles to be skipped: daughters from the top
void Matching::UpdateSkipList(GenParticle* top, vector<GenParticle>* genparticles, vector<bool>& skip)
{
	GenParticle* d1 = GetDaughter(top, genparticles, 1);
	GenParticle* d2 = GetDaughter(top, genparticles, 2);

	vector<GenParticle*> daughters;
	if (d1){
		skip[d1->index()-1] = true;
		daughters.push_back(d1);
	}
	if (d2){
		skip[d2->index()-1] = true;
		daughters.push_back(d2);
	}

	int nd = daughters.size();
	for (int i=0; i<nd; ++i){ // flag all further daughters from the daughters of the top
		GenParticle* p = daughters[i];
		GenParticle* d1 = GetDaughter(p, genparticles, 1);
		GenParticle* d2 = GetDaughter(p, genparticles, 2);
		if (d1){
		  if (!IsParton(d1)) continue;
			skip[d1->index()-1] = true;
			daughters.push_back(d1);
		}
		if (d2){
		  if (!IsParton(d2)) continue;
			skip[d2->index()-1] = true;
			daughters.push_back(d2);
		}
		nd = daughters.size();
	}
}
/*
███████ ███████ ████████ ████████ ███████ ██████
██      ██         ██       ██    ██      ██   ██
███████ █████      ██       ██    █████   ██████
     ██ ██         ██       ██    ██      ██   ██
███████ ███████    ██       ██    ███████ ██   ██
*/

// stores a list of stable particles (converted to pseudojets) for the clustering (particle level)
void Matching::set_stable_particles(vector<GenParticle>* genparticles)
{
	  for (unsigned i=0; i<genparticles->size(); ++i)
	    {
	      GenParticle* part = &(genparticles->at(i));
	      if (IsStableParticle(part) && !IsNeutrino(part)){
		      _stableParticles.push_back(convert_particle(part));
	      }
	    }
}
//stores a list of final state parton jets for the clustering (parton level)
void Matching::set_partons(vector<GenParticle>* genparticles)
{
	//check for final state partons from top decay
	    vector<GenParticle*> partons;
	    // here is a list of partons that should be skipped (daughters from top decay)
	    vector<bool> skip;
	    int nparts = genparticles->size();
	    skip.reserve(nparts);
	    //set all skips to false
	    for (int i=0; i<nparts; ++i) skip[i]=false;
	    //loop over all particles
	    for(int i=0; i<nparts; ++i){
	      // check if this parton should be skipped
	      if (skip[i]) continue; //the particle should be skipped
	  	  GenParticle* p = &(genparticles->at(i));
	  	  // is it a parton?
	      if (!IsParton(p)) continue; //if it is a particle->continue with next particle
	      // check if it's a top quark
	      if (IsTop(p)){
	  	     // is it the one that decays?
	  	      if (BeforeTopDecay(p,genparticles)){
	  	            // keep it
	  	            partons.push_back(p);
	  	            // flag all daughter partons to be skipped
	  	            UpdateSkipList(p, genparticles, skip);
	                continue;
	  	          }
	              else {
	  	            continue;
	  	          }
	      }
	      if (FinalStateParton(p, genparticles)){ //add all final state partons to the partons list
	  	     partons.push_back(p);
	      }
	    }

	//save the partons that should be clustered
	     for(unsigned i=0;i<partons.size();i++){
	       GenParticle* p = partons[i];
	       PseudoJet pseudojet=convert_particle(p);
	       if(abs(partons[i]->pdgId())==6 && IsHadronic(partons[i],genparticles)) pseudojet.set_user_index(6); //hadronically decaying top
	       else  pseudojet.set_user_index(0);
	       _partons_to_cluster.push_back(pseudojet);
	     }
}
//stores a list of final state W partons for the clustering (parton level)
void Matching::set_W_partons(vector<GenParticle>* genparticles)
{
	//check for final state partons from top decay
	    vector<GenParticle*> partons;
	    // here is a list of partons that should be skipped (daughters from top decay)
	    vector<bool> skip;
	    int nparts = genparticles->size();
	    skip.reserve(nparts);
	    //set all skips to false
	    for (int i=0; i<nparts; ++i) skip[i]=false;
	    //loop over all particles
	    for(int i=0; i<nparts; ++i){
	      // check if this parton should be skipped
	      if (skip[i]) continue; //the particle should be skipped
	  	  GenParticle* p = &(genparticles->at(i));
	  	  // is it a parton?
	      if (!IsParton(p)) continue; //if it is a particle->continue with next particle
	      // check if it's a top quark
	      if (IsTop(p)){
	  	     // is it the one that decays?
	  	      if (BeforeTopDecay(p,genparticles)){
							GenParticle* W_cand;
							GenParticle* b_cand;
									// get the W
									W_cand = GetDaughter(p, genparticles, 1);
                  b_cand = GetDaughter(p, genparticles, 2);
									if (!IsW(W_cand)) {
										W_cand = GetDaughter(p, genparticles, 2);
										b_cand = GetDaughter(p, genparticles, 1);
									}
									// std::cout << "----------Matching----------" << '\n';
									// std::cout << "top pt = " << p->pt() << '\n';
									// std::cout << "W pt = " << W_cand->pt() << '\n';
									// std::cout << "b pt = " << b_cand->pt() << '\n';

									//remove all final state partons from the b decay
									UpdateSkipList(b_cand, genparticles, skip);
                //  if (BeforeWDecay(W_cand, genparticles)) { // is it the W that decays? // TODO implement this method
                    // keep it
                   partons.push_back(W_cand);
                   // flag all daughter partons to be skipped
                   UpdateSkipList(W_cand, genparticles, skip);
                   continue;
								// } // end W that decays
                  // else {
									// 	// TODO take the daughter of the W (again check if this is the one that decays)
                  //   continue;
                  // }
	  	          } // end top that decays
	              else {
	  	            continue;
	  	          }
	      }
	      if (FinalStateParton(p, genparticles)){ //add all final state partons to the partons list
	  	     partons.push_back(p);
	      }
	    }

	//save the partons that should be clustered
	     for(unsigned i=0;i<partons.size();i++){
	       GenParticle* p = partons[i];
	       PseudoJet pseudojet=convert_particle(p);
	       if(IsW(p) && IsHadronic(p,genparticles)) pseudojet.set_user_index(6); //hadronically decaying W
	       else  pseudojet.set_user_index(0);
	       _W_partons_to_cluster.push_back(pseudojet);
	     }
}
//stores a list of final state b partons for the clustering (parton level)
void Matching::set_b_partons(vector<GenParticle>* genparticles)
{
	//check for final state partons from top decay
	    vector<GenParticle*> partons;
	    // here is a list of partons that should be skipped (daughters from top decay)
	    vector<bool> skip;
	    int nparts = genparticles->size();
	    skip.reserve(nparts);
	    //set all skips to false
	    for (int i=0; i<nparts; ++i) skip[i]=false;
	    //loop over all particles
	    for(int i=0; i<nparts; ++i){
	      // check if this parton should be skipped
	      if (skip[i]) continue; //the particle should be skipped
	  	  GenParticle* p = &(genparticles->at(i));
	  	  // is it a parton?
	      if (!IsParton(p)) continue; //if it is a particle->continue with next particle
	      // check if it's a top quark
	      if (IsTop(p)){
	  	     // is it the one that decays?
	  	      if (BeforeTopDecay(p,genparticles)){
							GenParticle* W_cand;
							GenParticle* b_cand;
									// get the W
									W_cand = GetDaughter(p, genparticles, 1);
                  b_cand = GetDaughter(p, genparticles, 2);
									if (!IsW(W_cand)) {
										W_cand = GetDaughter(p, genparticles, 2);
										b_cand = GetDaughter(p, genparticles, 1);
									}
									// std::cout << "----------Matching----------" << '\n';
									// std::cout << "top pt = " << p->pt() << '\n';
									// std::cout << "W pt = " << W_cand->pt() << '\n';
									// std::cout << "b pt = " << b_cand->pt() << '\n';

									//remove all final state partons from the b decay
									UpdateSkipList(W_cand, genparticles, skip);
                    // keep it
                   partons.push_back(b_cand);
                   // flag all daughter partons to be skipped
                   UpdateSkipList(b_cand, genparticles, skip);
                   continue;
	  	          } // end top that decays
	              else {
	  	            continue;
	  	          }
	      }
	      if (FinalStateParton(p, genparticles)){ //add all final state partons to the partons list
	  	     partons.push_back(p);
	      }
	    }
	//save the partons that should be clustered
	     for(unsigned i=0;i<partons.size();i++){
	       GenParticle* p = partons[i];
	       PseudoJet pseudojet=convert_particle(p);
				 if(Isb(p)) pseudojet.set_user_index(6); //b particle
		 		else  pseudojet.set_user_index(0);
				_b_partons_to_cluster.push_back(pseudojet);
	     }
}
/*
██████  ██    ██ ███    ██     ███    ███  █████  ████████  ██████ ██   ██ ██ ███    ██  ██████
██   ██ ██    ██ ████   ██     ████  ████ ██   ██    ██    ██      ██   ██ ██ ████   ██ ██
██████  ██    ██ ██ ██  ██     ██ ████ ██ ███████    ██    ██      ███████ ██ ██ ██  ██ ██   ███
██   ██ ██    ██ ██  ██ ██     ██  ██  ██ ██   ██    ██    ██      ██   ██ ██ ██  ██ ██ ██    ██
██   ██  ██████  ██   ████     ██      ██ ██   ██    ██     ██████ ██   ██ ██ ██   ████  ██████
*/
// match the parton level jets containing a top to the closest particle level jet
void Matching::run_matching(vector<TopJet> particle_jets, vector<TopJet> parton_jets)
{
  // clear all used vectors
_matched_jets.clear();
_matched_parton_jets.clear();
_matched_pairs.clear();

	std::vector<TopJet> jets = particle_jets;

//loop over parton jets
	for(uint j=0; j<parton_jets.size(); j++){
    if(jets.size()==0) continue; // skip empty jets

		TopJet matched_jet;
	  double minDeltaR=1000;
	  double delta_R;
//loop over particle jets, find closest
	  for(uint i=0; i<jets.size(); i++){
	    delta_R=deltaR(jets[i],parton_jets[j]);
	    if(delta_R<minDeltaR){
	      minDeltaR=delta_R;
	      matched_jet=jets[i];
				jets.erase(jets.begin()+i);
	    }
	  }
		// set the matching radius
		// double rho = 600;
		// double pt = matched_jet.pt();
		// double matching_radius = rho/pt;

    double matching_radius= 2; // take max distance of the HOTVR jet as the matching radius

	//	double matching_radius= matched_jet.max_distance(); // take max distance of the HOTVR jet as the matching radius
		// if(matching_radius<0.1){matching_radius=0.1;};
		// if(matching_radius>1.5){matching_radius=1.5;};

	//store the matched jets if it is close enough
		if (IsMatched(matched_jet,matching_radius,parton_jets[j])){
      _matched_jets.push_back(matched_jet);
			_matched_parton_jets.push_back(parton_jets[j]);
			_matched_pairs.push_back(make_pair(matched_jet, parton_jets[j]));
    }
  }
}

// match the parton level jets containing a top to the closest particle level jet
void Matching::run_matching_W_top(vector<TopJet> particle_jets, vector<TopJet> parton_jets_top, vector<TopJet> parton_jets_W)
{
  // clear all used vectors
_matched_jets.clear();
_matched_parton_jets.clear();
_matched_pairs.clear();

	std::vector<TopJet> jets = particle_jets;

//loop over parton jets
	for(uint j=0; j<parton_jets_W.size(); j++){
    if(jets.size()==0) continue; // skip empty jets

		TopJet matched_jet;
	  double minDeltaR=1000;
	  double delta_R;
//loop over particle jets, find closest
	  for(uint i=0; i<jets.size(); i++){
	    delta_R=deltaR(jets[i],parton_jets_W[j]);
	    if(delta_R<minDeltaR){
	      minDeltaR=delta_R;
	      matched_jet=jets[i];
				jets.erase(jets.begin()+i);
	    }
	  }

		double matching_radius= matched_jet.max_distance(); // take max distance of the HOTVR jet as the matching radius
		// if(matching_radius<0.1){matching_radius=0.1;};
		// if(matching_radius>1.5){matching_radius=1.5;};

	//store the matched jets if it is close enough
		if (IsMatched(matched_jet,matching_radius,parton_jets_W[j])){
      _matched_jets.push_back(matched_jet);
			_matched_parton_jets.push_back(parton_jets_W[j]);
			_matched_pairs.push_back(make_pair(matched_jet, parton_jets_W[j]));
    }
  }
}


// run matching for pseudojets (match to parton level AK4 jets)
void Matching::run_matching(std::vector<fastjet::PseudoJet> jets,std::vector<fastjet::PseudoJet> denominator_jets){
_matched_pseudojets.clear();
fastjet::PseudoJet matched_jet;
double minDeltaR=1000;
  double deltaR;
	for(uint j=0;j<denominator_jets.size();j++){
		for(uint i=0;i<jets.size();i++){ // get the closest jet
    	deltaR=jets[i].delta_R(denominator_jets[j]);
    	if(deltaR<minDeltaR) {
      	minDeltaR=deltaR;
      	matched_jet=jets[i];
    	}
 		}
		// now set the matching radius
		double rho = 600;
		double pt = matched_jet.pt();
		double matching_radius = rho/pt;
		if(matching_radius<0.1){matching_radius=0.1;};
		if(matching_radius>1.5){matching_radius=1.5;};
		// check if the jet is closer than the matching radius
		if (IsMatched(matched_jet,matching_radius,denominator_jets[j]))
			{
			_matched_pseudojets.push_back(matched_jet);
		}
	}
}
// for pseudojets with storing the jet constituents
void Matching::run_matching(std::vector<fastjet::PseudoJet> jets, vector<vector<PseudoJet>> constituents, std::vector<fastjet::PseudoJet> denominator_jets){
_matched_pseudojets.clear();
fastjet::PseudoJet matched_jet;
std::vector<fastjet::PseudoJet> jet_constituents;
double minDeltaR=1000;
  double deltaR;
	for(uint j=0;j<denominator_jets.size();j++){
		for(uint i=0;i<jets.size();i++){ // loop over jets
    	deltaR=jets[i].delta_R(denominator_jets[j]);
    	if(deltaR<minDeltaR) {
      	minDeltaR=deltaR;
      	matched_jet=jets[i];
				jet_constituents=constituents[i];
    	}
 		} // end loop jets
	//	if(!matched_jet.has_user_info<HOTVRinfo>()) matched_jet.set_user_info(new HOTVRinfo(matched_jet,matched_jet.constituents()));
		// set the matching radius
		double rho = 600;
		double pt = matched_jet.pt();
		double matching_radius = rho/pt;
		if(matching_radius<0.1){matching_radius=0.1;};
		if(matching_radius>1.5){matching_radius=1.5;};
		// check if the jet is close enough
		if (IsMatched(matched_jet,matching_radius,denominator_jets[j]))
			{
			_matched_pseudojets.push_back(matched_jet);
			_matched_pseudojets_and_constituents.push_back(make_pair(matched_jet, jet_constituents));
		}
	}// end loop denominator_jets
}

TopJet Matching::get_closest_jet(double radius, TopJet jet, std::vector<TopJet> parton_jets)
{
   TopJet closest_jet;
  double minDeltaR=radius;
  double delta_R;
//loop over particle jets, find closest
  for(uint i=0; i<parton_jets.size(); i++){
    delta_R=deltaR(jet,parton_jets[i]);
    if(delta_R<minDeltaR){
      minDeltaR=delta_R;
      closest_jet=parton_jets[i];
    }
  }
   return closest_jet;
}


//is the matched jet close enough?
bool Matching::IsMatched(TopJet jet, double matching_radius, TopJet denominator_jet)
{
  double delta_R;
  delta_R=deltaR(jet, denominator_jet);
  if(delta_R<matching_radius) return true;
  else return false;
}

//is the matched jet close enough?
bool Matching::IsMatched(fastjet::PseudoJet jet, double matching_radius, fastjet::PseudoJet denominator_jet)
{
  double delta_R;
  delta_R=jet.delta_R(denominator_jet);
  if(delta_R<matching_radius) return true;
  else return false;
}
// reset all the info, clear vectors
void Matching::Reset(){
  _stableParticles.clear();
  _partons_to_cluster.clear();
  _matched_jets.clear();
  _matched_parton_jets.clear();
  _matched_pairs.clear();
  _matched_pseudojets.clear();
  _matched_pseudojets_and_constituents.clear();
}
