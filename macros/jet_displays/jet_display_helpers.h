#ifndef __jet_display_helpers_h_
#define __jet_display_helpers_h_

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <TKey.h>
#include <TString.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TColor.h>

static UInt_t iseed;
static Bool_t is_new_event;

std::vector< std::vector<TH2F*> > GetSubjets(TFile* file, TString gdir);

std::vector<TH2F*> GetJets(TFile* file, TString gdir);

TGraph* GetParticles(TFile* file, TString gdir);

void NewEvent(UInt_t ievent);

UInt_t GetEventNumber(TString dir);

void Smear(double &x, double dx);

void DrawTops(TFile* file, TString gdir);

TGraph* GetTopDecay(TFile* file, TString gdir);

TGraph* GetTopDecay_specific(TFile* file, TString decayname, TString gdir);

void DrawTopDecay(TFile* file, TString gdir);

void DrawTopDecay_withb(TFile* file, TString gdir);

void DrawAntiTopDecay_withb(TFile* file, TString gdir);

void DrawStableParts(TFile* file, TString gdir);

TH2F* GetClustersRejected(TFile* file, TString gdir);

TH2F* GetJetsWithoutSubjets(TFile* file, TString gdir);

TH2F* GetJetsWithLowPt(TFile* file, TString gdir);

void PrintInfo(TFile* file, TString fname, TString gdir);

TH2F* GetSubjetsRejectedBySD(TFile* file, TString gdir);

void RemovePartons(TH2* jh, TFile* file, TString gdir);

bool CheckContent(TH2F* h);

void CreateContours(std::vector<TH2F*> jets);

void CreateContour(TH2F* h);

void DrawContour(TH2F* hist, int col);

std::vector< std::vector<int> > GetColours();

void FixHistForContour(TH2F* h);

#endif //__jet_display_helpers_h_


