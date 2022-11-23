#ifndef __jet_display_hotvr_h_
#define __jet_display_hotvr_h_

#include <TFile.h>
#include <TString.h>

void jet_display_hotvr(TFile* file, TString dir, TString prefix = "", Int_t num=-1, Bool_t bLargeJetsOnly = false);

#endif // __jet_display_hotvr_h_