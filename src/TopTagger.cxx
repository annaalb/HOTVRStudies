#include "UHH2/HOTVRStudies/include/TopTagger.h"

using namespace std;
using namespace fastjet;
using namespace contrib;

TopTagger::TopTagger()
{
  //?
}

bool TopTagger::Is_tagged(string tagger, TopJet &jet)
{
  if(tagger=="HOTVR"||tagger=="hotvr") return Is_HOTVR_tagged(jet);
  else if(tagger=="SD"||tagger=="Soft Drop" ||tagger=="softdrop" || tagger=="sd") return Is_SD_tagged(jet);
  else if(tagger=="W") return Is_W_tagged(jet);
  else if(tagger=="Z") return Is_Z_tagged(jet);
  else if(tagger=="H") return Is_H_tagged(jet);
  else return false;
}

bool TopTagger::Is_HOTVR_tagged(TopJet &jet)
{
  double mjet=jet.v4().M();
  int nsubjets=jet.subjets().size();
  double ptfraction=jet.hotvr_fpt1();
  double mmin=jet.hotvr_mmin();
  bool tagged=false;
//std::cout << "Mjet "<<mjet<< " nsubjets "<<nsubjets<< " mmin "<<mmin<<" ptf "<<ptfraction << '\n';
  if(mjet>140 && mjet<220 && nsubjets>2 && mmin>50 && ptfraction<0.8) tagged=true;
  return tagged;
}

bool TopTagger::Is_SD_tagged(TopJet &jet)
{ // TODO : choose tagging parameters for HOTVR SD
  double mjet=jet.v4().M();
  int nsubjets=jet.subjets().size();
  double ptfraction=jet.hotvr_fpt1();
  double mmin=jet.hotvr_mmin();
  bool tagged=false;
//std::cout << "Mjet "<<mjet<< " nsubjets "<<nsubjets<< " mmin "<<mmin<<" ptf "<<ptfraction << '\n';
  if(mjet>140 && mjet<220 && nsubjets>2 && mmin>50 && ptfraction<0.8) tagged=true;
  return tagged;
}
//--------same for pseudojets-------------
bool TopTagger::Is_tagged(string tagger, PseudoJet &jet)
{
  if(tagger=="HOTVR"||tagger=="hotvr") return Is_HOTVR_tagged(jet);
  else if(tagger=="SD"||tagger=="Soft Drop" ||tagger=="softdrop" || tagger=="sd") return Is_SD_tagged(jet);
  else return false;
}

bool TopTagger::Is_HOTVR_tagged(PseudoJet &jet)
{
  double mjet=jet.m();
  int nsubjets=jet.user_info<HOTVRinfo>().nsubjets();
  double ptfraction=jet.user_info<HOTVRinfo>().ptfraction(0);
  double mmin=jet.user_info<HOTVRinfo>().mmin();
  bool tagged=false;
  if(mjet>140 && mjet<220 && nsubjets>2 && mmin>50 && ptfraction<0.8) tagged=true;
  return tagged;
}

bool TopTagger::Is_SD_tagged(PseudoJet &jet)
{ // TODO : choose tagging parameters for HOTVR SD
  double mjet=jet.m();
  int nsubjets=jet.user_info<HOTVRinfo>().nsubjets();
  double ptfraction=jet.user_info<HOTVRinfo>().ptfraction(0);
  double mmin=jet.user_info<HOTVRinfo>().mmin();
  bool tagged=false;
  if(mjet>140 && mjet<220 && nsubjets>2 && mmin>50 && ptfraction<0.8) tagged=true;
  return tagged;
}

bool TopTagger::Is_W_tagged(TopJet &jet)
{
  double mjet=jet.v4().M();
  int nsubjets=jet.subjets().size();
  bool tagged=false;
  if(mjet>65 && mjet<95) tagged=true;
  return tagged;
}

bool TopTagger::Is_Z_tagged(TopJet &jet)
{
  double mjet=jet.v4().M();
  int nsubjets=jet.subjets().size();
  bool tagged=false;
  if(mjet>75 && mjet<105) tagged=true;
  return tagged;
}

bool TopTagger::Is_H_tagged(TopJet &jet)
{
  double mjet=jet.v4().M();
  int nsubjets=jet.subjets().size();
  bool tagged=false;
  if(mjet>110 && mjet<140) tagged=true;
  return tagged;
}
