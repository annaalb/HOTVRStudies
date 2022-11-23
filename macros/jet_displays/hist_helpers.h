#ifndef __hist_helpers_h_
#define __hist_helpers_h_

#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>
#include <map>
#include <TH2F.h>
#include <TGraph.h>
#include <TColor.h>
#include <TString.h>
#include "style.h"
#include "jet_display_helpers.h"

struct part{
	uint ind; 
	float pt;
	float eta;
	float phi;
	int color;
};

struct leave{
	uint ind;
	int parent1;
	int parent2;
	int child;
	double dij;
};

extern TString _jet_cols_fname; 
extern uint _event;

std::vector<part> GetInputParticles(TString fname);

std::vector<leave> GetHistoryTree(TString fname);

void PlotAxes();
void PlotPart(part p);

void PlotInputParts(std::vector<part> ips);

void DoClusteringStep(int step, std::vector<part>& ips, std::vector<leave> tree, std::map<int, TH2F*>& chists);

TH2F* ClusterTwoParticles(std::vector<part>& ips, leave tl);

void ClusterParticleAndCluster(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists);

void ClusterTwoClusters(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists);

void CreateJet(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists);

Int_t GetNextColor();

Int_t GetJetColor(TH2F* jethist);

int CompareHists(TH2F* h1, TH2F* h2);

inline void SetJetColsFile(TString in){_jet_cols_fname = in;}

inline void SetJetColsEvent(uint in){_event = in;}

void DrawAll(std::vector<part> ips, std::map<int, TH2F*> chists);

void PrintStep(int step, std::vector<part> ips, std::vector<leave> tree);

#endif //__hist_helpers_h_














