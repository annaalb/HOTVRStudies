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

// identify neutrinos
bool Matching::IsNeutrino(GenParticle* p)
{
  int id = abs(p-> pdgId());
  if (id==12 || id==14 || id==16) return true;
  else return false;
}

// identify tops
bool Matching::IsTop(GenParticle* p)
{
  int id = p->pdgId();
  if (abs(id)==6) return true;
  else return false;
}

// identify W bosons
bool Matching::IsW(GenParticle* p)
{
  int id = p->pdgId();
  if (abs(id)==24) return true;
  else return false;
}

// identify Z bosons
bool Matching::IsZ(GenParticle* p)
{
  int id = p->pdgId();
  if (abs(id)==23) return true;
  else return false;
}

// identify Higgs bosons
bool Matching::IsH(GenParticle* p)
{
  int id = p->pdgId();
  if (abs(id)==25) return true;
  else return false;
}

// Identify b quarks
bool Matching::Isb(GenParticle* p)
{
  int id = p->pdgId();
  if (abs(id)==5) return true;
  else return false;
}

// checks for quarks, gluons, W, Z, H bosons
bool Matching::IsParton(GenParticle* p)
{
	int fid = abs(p->pdgId());
	bool isp = false;
	if (fid<7)   isp = true;
	if (fid==21) isp = true;
	if (fid==23) isp = true;
	if (fid==24) isp = true;
	if (fid==25) isp = true;
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

// check if this W or Z decays into qq
bool Matching::BeforeWZDecay(GenParticle* p, vector<GenParticle>* genparticles)
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

  // condition: both daughters are quarks or leptons (not photons, Z or W)
	if (id1<19  && id2<19) flag = true; // decay to quarks

	return flag;
}

// check if this H decays into qq or gg (let's not consider more complicated 4q decays)
bool Matching::BeforeHDecay(GenParticle* p, vector<GenParticle>* genparticles)
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

  // condition: both daughters are quarks or gluons
	if (id1<19  && id2<19) flag = true; // decay to quarks (mostly bb, but cc also allowed) or leptons (mostly tautau)
	if (id1==21  && id2==21) flag = true; // decay to gluons

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
bool Matching::IsHadronicTop(GenParticle* p, vector<GenParticle>* genparticles)
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

//check if the WZH decay hadronically
bool Matching::IsHadronicWZH(GenParticle* p, vector<GenParticle>* genparticles)
{
  GenParticle* d1 = GetDaughter(p, genparticles, 1);
  GenParticle* d2 = GetDaughter(p, genparticles, 2);

  // check ID of daughters
	int id1 = abs(d1->pdgId());
	int id2 = abs(d2->pdgId());

  bool flag = false;
  // condition: both daughters are quarks or gluons
	if (id1<6  && id2<6) flag = true; // decay to quarks (mostly bb, but cc also allowed)
	if (id1==21  && id2==21) flag = true; // decay to gluons

  return flag;
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

// go through all decays from the input particle (top here, but can be W, Z, H)
// and flag those to be skipped
// particles to be skipped: all daughters from the decay of the input particle
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

// check what event we analyse: ttbar, WW, ZZ, HH or QCD
// returns the PDGID of the main particle. Returns 0 if no heavy particle is found (=QCD)
int Matching::PrimaryEventPDGID(std::vector<GenParticle> *genparticles)
{
  bool found_top = false;
  bool found_W = false;
  bool found_Z = false;
  bool found_H = false;

  int nparts = genparticles->size();

  //loop over all particles
  for(int i=0; i<nparts; ++i){

    GenParticle* p = &(genparticles->at(i));

    // is it a parton?
    if (!IsParton(p)) continue;

    // check what particle it is
    if (IsTop(p)) found_top = true;
    else if (IsZ(p)) found_Z = true;
    else if (IsW(p)) found_W = true;
    else if (IsH(p)) found_H = true;
    else continue;
  }

  // note: order is important, since the Higgs can decay to WW or ZZ
  int id = 0;
  if (found_top) id = 6;
  else if (found_H) id = 25;
  else if (found_Z) id = 23;
  else if (found_W) id = 24;
  else id = 0;

  return id;

}

// check if the event is fully hadronic, works on every type of event
bool Matching::IsGravitonEvent(std::vector<GenParticle> *genparticles)
{
  bool is_grav = false;

  //loop over all particles
  int nparts = genparticles->size();
  for(int i=0; i<nparts; ++i){
    GenParticle* p = &(genparticles->at(i));
    if (p->pdgId() > 5000000) is_grav = true;
  }
  return is_grav;
}

// check if the event is fully hadronic, works on every type of event
bool Matching::IsFullyHadronicEvent(std::vector<GenParticle> *genparticles)
{
  bool event_is_had = true;

  //std::cout << "Check if it's fully hadronic" << std::endl;

  //loop over all particles
  int nparts = genparticles->size();
  for(int i=0; i<nparts; ++i){

    GenParticle* p = &(genparticles->at(i));

    // is it a parton?
    if (!IsParton(p)) continue;

    // check what particle it is
    // we are only interested in W, Z, H and top
    int pdgid=0;
    if (IsTop(p)) pdgid = 6;
    else if (IsZ(p)) pdgid = 23;
    else if (IsW(p)) pdgid = 24;
    else if (IsH(p)) pdgid = 25;
    else {
      continue;
    }

    // if it's W, Z, H and top, check if it's the particle before the decay
    if (pdgid>0){
       // is it the one that decays?
       bool beforedecay=false;
       if (pdgid==6) beforedecay = BeforeTopDecay(p, genparticles);
       else if (pdgid==23) beforedecay = BeforeWZDecay(p, genparticles);
       else if (pdgid==24) beforedecay = BeforeWZDecay(p, genparticles);
       else if (pdgid==25) beforedecay = BeforeHDecay(p, genparticles);

       //std::cout << "here is the W/Z/H/top we are looking at: " << std::endl;
       //p->Print(genparticles);
       //std::cout << "Is it before its decay? " << beforedecay << std::endl;

       if (beforedecay){
         bool had_decay = false;
         // check if it decays hadronically
         if (pdgid==6) had_decay = IsHadronicTop(p,genparticles);
         else if (pdgid==23) had_decay = IsHadronicWZH(p, genparticles);
         else if (pdgid==24) had_decay = IsHadronicWZH(p, genparticles);
         else if (pdgid==25) had_decay = IsHadronicWZH(p, genparticles);

         //std::cout << "does it decay hadronically? " << had_decay << std::endl;

         // if we found a particle which does not decay hadronically, flag the event as "not-had"
         if (!had_decay){
           event_is_had = false;
         }

        } else {
          continue;
        }
      } // pdg ID of W, Z, H, top
    } // loop over all particles

  return event_is_had;

}


// stores a list of final state parton jets for the clustering (parton level)
// the target is the PDG ID of the particle that we are looking for
// for example, if the final state should be top quarks + other partons, then
// the target = 6. This can be done for any target, like 24 (=W), 23 (=Z), 25 (=H)
// default value of target: 6 (see header)
// if all final state partons should be kept, the target should be 0 (used for QCD)
// if the target is set to 5 (=b), all final state partons are kept as well (same as for 0)
void Matching::set_partons(vector<GenParticle>* genparticles, int target)
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
        //p->Print(genparticles);
	  	  // is it a parton?
	      if (!IsParton(p)) continue; //if it is a particle->continue with next particle
	      // check if it's the target particle
        bool istarget=false;
        if (abs(target)==6) istarget = IsTop(p);
        else if (abs(target)==23) istarget = IsZ(p);
        else if (abs(target)==24) istarget = IsW(p);
        else if (abs(target)==25) istarget = IsH(p);
        else if (abs(target)==5) istarget = false; // special case for b: it does not decay, take the final state parton b
        else if (abs(target)==0) istarget = false;
        else {
          std::cout << "Matching::set_partons: Error, unknown target ID" << std::endl;
        }

        // if it's the target particle, keep it or go on (depending if it's the last copy of it)
	      if (istarget){
           // is it the one that decays?
           bool beforedecay=false;
           if (abs(target)==6) beforedecay = BeforeTopDecay(p, genparticles);
           else if (abs(target)==23) beforedecay = BeforeWZDecay(p, genparticles);
           else if (abs(target)==24) beforedecay = BeforeWZDecay(p, genparticles);
           else if (abs(target)==25) beforedecay = BeforeHDecay(p, genparticles);

	  	     if (beforedecay){
	  	        // it's the target before its decay: keep it
              //std::cout << "--> This is the particle before the decay. Keep!" << std::endl;
	  	        partons.push_back(p);
	  	        // flag all daughter partons to be skipped
	  	        UpdateSkipList(p, genparticles, skip);
	            continue;
	  	     } else {
	  	        continue;
	  	     }
	      }

        // still here? then it was not the target particle before its decay
	      if (FinalStateParton(p, genparticles)){ //add all final state partons to the partons list
	  	     partons.push_back(p);
	      }
	    } // end loop over all particles

      //std::cout << "\n\nHere is the list of selected particles: " << std::endl;

	    //save the partons that should be clustered
	    for(unsigned i=0;i<partons.size();i++){
	       GenParticle* p = partons[i];
         //p->Print(genparticles);
	       PseudoJet pseudojet=convert_particle(p);
	       if (abs(partons[i]->pdgId())==6 && IsHadronicTop(partons[i],genparticles)) pseudojet.set_user_index(6); //hadronically decaying top
	       else if (abs(partons[i]->pdgId())==5) pseudojet.set_user_index(5); // b quark before its fragmentation
         else if (abs(partons[i]->pdgId())==23 && IsHadronicWZH(partons[i],genparticles)) pseudojet.set_user_index(23); //hadronically decaying Z
         else if (abs(partons[i]->pdgId())==24 && IsHadronicWZH(partons[i],genparticles)) pseudojet.set_user_index(24); //hadronically decaying W
         else if (abs(partons[i]->pdgId())==25 && IsHadronicWZH(partons[i],genparticles)) pseudojet.set_user_index(25); //hadronically decaying H
         else pseudojet.set_user_index(0);
	       _partons_to_cluster.push_back(pseudojet);
	     }
}

//stores a list of final state W partons for the clustering (parton level)
void Matching::set_Wb_partons_from_top(vector<GenParticle>* genparticles)
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

                  // keep the W and the b
                  partons.push_back(W_cand);
                  partons.push_back(b_cand);
                  // remove all final state partons from the W and b decays (avoid double counting, since we store the two)
                  UpdateSkipList(W_cand, genparticles, skip);
									UpdateSkipList(b_cand, genparticles, skip);
                  // also remove the W and b, such that we don't store them by accident
                  skip[W_cand->index()-1] = true;
                  skip[b_cand->index()-1] = true;
                  continue;
	  	      } // end top that decays
	          else {
	  	        continue;
	  	      }
	      } // end if(top)

	      if (FinalStateParton(p, genparticles)){ //add all final state partons to the partons list
	  	     partons.push_back(p);
	      }
	    }

	    // save the partons that should be clustered
	    for(unsigned i=0;i<partons.size();i++){
	       GenParticle* p = partons[i];
	       PseudoJet pseudojet=convert_particle(p);
	       if(IsW(p) && IsHadronicWZH(p,genparticles)) pseudojet.set_user_index(24); // hadronically decaying W
	       if(Isb(p)) pseudojet.set_user_index(5); // b quark
	       else pseudojet.set_user_index(0);
	       _Wb_partons_to_cluster.push_back(pseudojet);
	     }
}

/*
██████  ██    ██ ███    ██     ███    ███  █████  ████████  ██████ ██   ██ ██ ███    ██  ██████
██   ██ ██    ██ ████   ██     ████  ████ ██   ██    ██    ██      ██   ██ ██ ████   ██ ██
██████  ██    ██ ██ ██  ██     ██ ████ ██ ███████    ██    ██      ███████ ██ ██ ██  ██ ██   ███
██   ██ ██    ██ ██  ██ ██     ██  ██  ██ ██   ██    ██    ██      ██   ██ ██ ██  ██ ██ ██    ██
██   ██  ██████  ██   ████     ██      ██ ██   ██    ██     ██████ ██   ██ ██ ██   ████  ██████
*/

// returns the index of the closest jet in the list jets, calculated with respect to the
// momentum vector of the reference jet refjet.
// the metric for deciding the "closeness" is DeltaR in eta and phi
uint Matching::i_closest_particlejet(TopJet refjet, vector<TopJet> jets)
{
  uint iclosest = 999;
  double minDeltaR=1000;
  double delta_R;

  // loop over jets to find the one that is closest to the refjet
  for(uint i=0; i<jets.size(); i++){
    delta_R = deltaR(jets[i],refjet);
    if(delta_R<minDeltaR){
      minDeltaR = delta_R;
      iclosest = i;
    }
  }
  return iclosest;
}

TopJet Matching::get_closest_jet(TopJet refjet, std::vector<TopJet> jets)
{
  TopJet closest_jet;
  double minDeltaR=1000;
  double delta_R;
  // loop over jets, find closest to the refjet
  for(uint i=0; i<jets.size(); i++){
    delta_R=deltaR(jets[i], refjet);
    if(delta_R<minDeltaR){
      minDeltaR=delta_R;
      closest_jet=jets[i];
    }
  }
  return closest_jet;
}

//is the matched jet close enough?
bool Matching::IsMatched(TopJet jet, double matching_radius, TopJet refjet)
{
  double delta_R;
  delta_R=deltaR(jet, refjet);
  if(delta_R<matching_radius) return true;
  else return false;
}

//is the matched jet close enough?
bool Matching::IsMatched(fastjet::PseudoJet jet, double matching_radius, fastjet::PseudoJet refjet)
{
  double delta_R;
  delta_R=jet.delta_R(refjet);
  if(delta_R<matching_radius) return true;
  else return false;
}


// match the parton level jets to the closest particle level jets
void Matching::run_matching(vector<TopJet> particle_jets, vector<TopJet> parton_jets)
{
  // clear all used vectors
  _matched_jets.clear();
  _matched_parton_jets.clear();
  _matched_pairs.clear();

	std::vector<TopJet> jets = particle_jets;

  // loop over parton jets
	for(uint j=0; j<parton_jets.size(); j++){

    // get the closest particle level jet
    uint i = i_closest_particlejet(parton_jets[j], particle_jets);
    if (i>900) continue; // no matching jet found

		TopJet matched_jet = particle_jets[i];

		// set the matching radius
		// double rho = 600;
		// double pt = matched_jet.pt();
		// double matching_radius = rho/pt;
    //	double matching_radius= matched_jet.max_distance(); // take max distance of the HOTVR jet as the matching radius
		// if(matching_radius<0.1){matching_radius=0.1;};
		// if(matching_radius>1.5){matching_radius=1.5;};

    // check if the closest jet is "close enough"
    double matching_radius = 1.5; // allow for a large matching radius

	  // store the matched jet if it is close enough
		if (IsMatched(matched_jet,matching_radius,parton_jets[j])){
      _matched_jets.push_back(matched_jet);
			_matched_parton_jets.push_back(parton_jets[j]);
			_matched_pairs.push_back(make_pair(matched_jet, parton_jets[j]));
      // remove the matched jet, such that it can not be matched again
      jets.erase(jets.begin()+i);
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
void Matching::run_matching(std::vector<fastjet::PseudoJet> jets, vector<vector<PseudoJet>> constituents, std::vector<fastjet::PseudoJet> denominator_jets)
{
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
