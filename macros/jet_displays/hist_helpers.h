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

void PlotPart(part p, int color=0);

void PlotInputParts(std::vector<part> ips, int color=0);

bool DoClusteringStep(int step, std::vector<part>& ips, std::vector<leave> tree, std::map<int, TH2F*>& chists);

bool ClusterTwoParticles(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists);

bool ClusterParticleAndCluster(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists);

bool ClusterTwoClusters(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists);

bool CreateJet(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists);

void UpdateParticleColors(int ch, std::vector<part>& ips, std::vector<leave> tree, std::map<int, TH2F*> chists, int col=0);

Int_t GetNextColor();

Int_t GetJetColor(TH2F* jethist);

int CompareHists(TH2F* h1, TH2F* h2);

inline void SetJetColsFile(TString in){_jet_cols_fname = in;}

inline void SetJetColsEvent(uint in){_event = in;}

void DrawAll(std::vector<part> ips, std::map<int, TH2F*> chists);

void PrintStep(int step, std::vector<part> ips, std::vector<leave> tree);

#endif //__hist_helpers_h_














