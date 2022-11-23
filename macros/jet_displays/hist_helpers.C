#include "hist_helpers.h"

std::vector<part> GetInputParticles(TString fname)
{

	using namespace std;
	vector<part> iparts; 
	ifstream file;
	file.open(fname.Data());
	string line;
	if ( file.is_open() ) {
	  while ( file ) {
		getline (file, line);
		//std::cout << line << '\n';		
		if (strncmp(&(line.front()),"#",1)==0){ 
			// ignore, this is a comment
		} else if (strncmp(&(line.front()),"\n",1)==0){
			// ignore, new line
		} else {
			part ap; 
			sscanf( line.data(), "%d %f %f  %f", &ap.ind, &ap.pt, &ap.eta, &ap.phi );
			if (ap.phi>M_PI) ap.phi -= 2*M_PI;
			//ap.color = kGray+2;
			ap.color = kBlack;
			//std::cout << "ind = " << ap.ind << ", pt = " << ap.pt << ", eta = " << ap.eta << ", phi = " << ap.phi << std::endl;
			if (iparts.size()>1){
				if (ap.ind == iparts.at(iparts.size()-1).ind){
				  // do nothing, we have this index already	
				} else {
				  iparts.push_back(ap);
				} 
			} else { // take the first one
			  iparts.push_back(ap); 
			}
		}
	  }
	
	} else {
		std::cout << "GetInputParticles: Couldn't open file " << fname << "\n";
	}
	return iparts;
}

std::vector<leave> GetHistoryTree(TString fname)
{
	using namespace std;

	vector<leave> tree; 

	ifstream file;
	file.open(fname.Data());
	string line;
	if ( file.is_open() ) {
	  while ( file ) {
		getline (file, line);
		//std::cout << line << '\n';		
		if (strncmp(&(line.front()),"#",1)==0){ 
			// ignore, this is a comment
		} else if (strncmp(&(line.front()),"\n",1)==0){
			// ignore, new line
		} else {
			leave ol; 
			sscanf( line.data(), "%d %d %d %d %lf", &ol.ind, &ol.parent1, &ol.parent2, &ol.child, &ol.dij );
			//std::cout << "ind = " << ol.ind << ", p1 = " << ol.parent1 << ", p2 = " << ol.parent2 << ", child = " << ol.child << ", dij = " << ol.dij << std::endl;
			if (tree.size()>1){
				if (ol.ind == tree.at(tree.size()-1).ind){
				  // do nothing, we have this index already	
				} else {
				  tree.push_back(ol);
				} 
			} else { // take the first one
			  tree.push_back(ol); 
			}
		}
	  }
	
	} else {
		std::cout << "GetInputParticles: Couldn't open file " << fname << "\n";
	}

	return tree;
}

void PlotAxes()
{
	static TH2F* h = new TH2F("h", "", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
	h->GetXaxis()->SetTitle("#phi");
	h->GetYaxis()->SetTitle("#eta");	
	Cosmetics2d(h);
	h->Draw();
}

void PlotPart(part p)
{
	TGraph* g = new TGraph(1); 
	g->SetPoint(0, p.phi, p.eta); 
	g->SetMarkerStyle(20);
	g->SetMarkerColor(p.color); 		
	double size = 0.4;
	double x = TMath::Log(p.pt);
	if (x>0) size += 0.2*x; // minimum size of 0.4, then increase with log(pt)
	g->SetMarkerSize(size);
	g->DrawClone("Psame");
	g->Delete();
}

void PlotInputParts(std::vector<part> ips)
{
	// convert ips to graphs with different marker sizes (depending on pT)
	// and plot these on the current pad

	for (int i=0; i<ips.size(); ++i){
		//std::cout << "Plotting part " << i << " with color " << ips[i].color << std::endl;
		PlotPart(ips[i]);
	}

} 

void UpdateColors(std::vector<part>& ips, int oldc, int newc)
{
	// update all particles with old color 'oldc'
	// to new color 'newc'
	for (int i=0; i<ips.size(); ++i){
		if (ips[i].color == oldc){
			ips[i].color = newc;
		}
	}

}

void DoClusteringStep(int step, std::vector<part>& ips, std::vector<leave> tree, std::map<int, TH2F*>& chists)
{
	// plot one step: check entry in history. start with Np (funny 0-counting in c++)
    // for each history entry: look up in the list of input particles which particles are combined 
    // create a histogram which combines the particles. 
    // store the histogram in a map with its identifier 
    // if there appears an entry in the history which has an ID >= Np (equal sign important), then check 
    // the map with histograms. Either combine a particle with an existing histogram, 
    // or combine two histograms (need two different functions for this)

    // todo: check if one parent becomes a negative number. then dib was smallest and the cluster became a jet. 

	using namespace std;

	// index in the clustering history, which we look at
	uint Np = ips.size();
	int ch = step - 1 + Np; // reduce step by 1, such that we can start with step 1

	// get the current place in the history
	leave tl = tree[ch]; 

	cout << "ch = " << ch << " parent1 = " << tl.parent1 << " parent2 = " << tl.parent2 << endl;

	// get the clusters we want to combine. three possibilities: 
	// (1) both are input particles, 
	// (2) one is an input particle, the other one a combined cluster
	// (3) both are combined clusters
	enum Combination {
      PP,
      PC,
      CC,
      JET
    };
    Combination comb; 

    int pa1 = abs(tl.parent1);
    int pa2 = abs(tl.parent2);

    if ((pa1 < Np) && (pa2 < Np)){
      // a single particle can become a jet if diB was smallest
      if (tl.parent2<0){
      	comb = JET;
      } else {
    	comb = PP; 
      }
    }
    if ( ((pa1 < Np) && (pa2 >= Np)) || ((pa1 >= Np) && (pa2 < Np)) ){
      // two possibilities: either combine a cluster with a particle, 
      // or a cluster becomes a jet
      if (tl.parent2<0){
      	comb = JET;
      } else {
        comb = PC; 
      }
    }  
    if ((pa1 >= Np) && (pa2 >= Np)){
      comb = CC;	
    } 


    switch (comb){
      case PP: {
      	cout << "clustering two particles" << endl;
      	// this is the easy case: we need to create a new histogram 
      	// which shows the clustering of these two
      	TH2F* combohist = ClusterTwoParticles(ips, tl);
      	// store combination in the map, with the child as unique index
      	chists[tl.ind] = combohist; 
      	break;
      }

      case PC: {
      	cout << "clustering a particle and a cluster" << endl;
      	// here, we will combine the particle with the histogram of the cluster. 
      	// always take the color of the histogram 
      	ClusterParticleAndCluster(ips, tl, chists);
      	break;
      }

      case CC: {
      	cout << "clustering two clusters" << endl;
      	// now, let's combine two histograms. the color should be the one 
      	// with the larger number (this cluster has seen larger dij values)
      	ClusterTwoClusters(ips, tl, chists);
      	break;
      }

  	  case JET: {
  	  	cout << "This cluster becomes a jet" << endl;
  	  	// create a jet by giving this a new color, chosen from the jet color profile
  	  	CreateJet(ips, tl, chists);
  	  	break;
      }

    }

}

TH2F* ClusterTwoParticles(std::vector<part>& ips, leave tl)
{
	// first: get the two particles we cluster and color them
    part p1 = ips[tl.parent1];
    part p2 = ips[tl.parent2];
    Int_t color = GetNextColor();
    p1.color = color;
    p2.color = color;

	// make sure to store the changes (only the color) in the list of inputs
    ips[tl.parent1].color = color;
    ips[tl.parent2].color = color;

    //cout << "p1.phi = " << p1.phi << " p1.eta = " << p1.eta << endl;
    //cout << "p2.phi = " << p1.phi << " p2.eta = " << p2.eta << endl;  

    TString hname = TString::Format("PPcombo_%d_%d", tl.parent1, tl.parent2);

    TH2F* h = new TH2F(hname, "", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
	Cosmetics2d(h);
	h->Fill(p1.phi, p1.eta);
	h->Fill(p2.phi, p2.eta);
	CreateContour(h);
	FixHistForContour(h);
	h->SetFillColor(color);
	h->SetLineColor(color);

	return h;

}

void ClusterParticleAndCluster(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists)
{
	// let's find out which one is the particle, and which one the cluster
	TH2F* cluster; 
	int cluster_id = 0;
	part p;
	int pid = 0; 
	uint Np = ips.size();
	if ( (tl.parent1 < Np) && (tl.parent2 >= Np) ){
		p = ips[tl.parent1];
		pid = tl.parent1;
		cluster = chists[tl.parent2];
		cluster_id = tl.parent2;
		//cout << "Case 1: Index of particle = " << p.ind << " index of cluster = " << tl.parent2 << " name of hist = " << cluster->GetName() << endl;
	} else if ((tl.parent1 >= Np) && (tl.parent2 < Np)) {
		p = ips[tl.parent2];
		pid = tl.parent2;
		cluster = chists[tl.parent1];
		cluster_id = tl.parent1; 
		//cout << "Case 2: Index of particle = " << p.ind << " index of cluster = " << tl.parent2 << " name of hist = " << cluster->GetName() << endl;
	} else {
		std::cout << "ERROR in ClusterParticleAndCluster! parent 1: " << tl.parent1 << " parent 2:" << tl.parent2 << std::endl;
	}

	// now add the particle to the cluster
	p.color = cluster->GetFillColor();
	cluster->Fill(p.phi, p.eta);
	CreateContour(cluster);
	FixHistForContour(cluster);
	// rename the hist
    TString hname = TString::Format("PCcombo_%d_%d", tl.parent1, tl.parent2);
    cluster->SetName(hname);

    // now replace the histogram with the new one (unique index needs to change to the index of the leave)
    auto node = chists.extract(cluster_id);
    node.key() = tl.ind;
    chists.insert(std::move(node));


	// store the particle color in the list of input particles
	ips[pid].color = p.color; 
    
    return;
}


void ClusterTwoClusters(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists)
{
	// let's get the clusters
	TH2F* cluster1; 
	TH2F* cluster2; 

	cluster1 = chists[tl.parent1];
	cluster2 = chists[tl.parent2];

	// decide on a color
	Int_t col1 = cluster1->GetFillColor();
	Int_t col2 = cluster2->GetFillColor();	
	Int_t colnew;
	Int_t colold; 
	if (col1>col2){ 
		colnew = col1;
		colold = col2;
	} else {
		colnew = col2;
		colold = col1;
	}

	// merge the two into one
	cluster1->Add(cluster2); 
	cluster1->SetFillColor(colnew);
	CreateContour(cluster1);
	FixHistForContour(cluster1);
	// rename the hist
    TString hname = TString::Format("CCcombo_%d_%d", tl.parent1, tl.parent2);
    cluster1->SetName(hname);

    // now replace the histogram of cluster1 with the new one (unique index needs to change to the index of the leave)
    // and remove the histogram 2
    auto node = chists.extract(tl.parent1);
    node.key() = tl.ind;
    chists.insert(std::move(node));
    chists.extract(tl.parent2);

    // update all particles to have the new color
    UpdateColors(ips, colold, colnew);

    return;
}

void CreateJet(std::vector<part>& ips, leave tl, std::map<int, TH2F*>& chists)
{

	// decide on a color
	Int_t col = 0;
	// give the jet a name
	TString jname = TString::Format("JET_%d", tl.ind);
	std::cout << "Creating jet " << jname << "." << std::endl;

	// first check if it's a special case, where a single input particle became a jet
	uint Np = ips.size();
    int pa1 = tl.parent1;

    if (pa1 < Np){
    	// create a histogram for this single-particle jet
	    ips[tl.parent1].color = col; // store the new color of this particle
    	part p = ips[tl.parent1];

    	// create a histogram for this jet
    	TH2F* h = new TH2F(jname, "", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
		Cosmetics2d(h);
		h->Fill(p.phi, p.eta);
		CreateContour(h);
		FixHistForContour(h);
		col = GetJetColor(h);
		h->SetFillColor(col);
		h->SetLineColor(col);
		chists[tl.ind] = h;



    } else {
    	// replace the histogram of cluster with the new one (unique index needs to change to the index of the leave)
    	auto node = chists.extract(tl.parent1);
    	node.key() = tl.ind;
    	chists.insert(std::move(node));
    	col = GetJetColor(chists.at(tl.ind));
    	int colold = chists[tl.ind]->GetFillColor();
    	chists[tl.ind]->SetFillColor(col);
    	chists[tl.ind]->SetLineColor(col);    
    	chists[tl.ind]->SetName(jname);

    	// update all particles to have the new jet color
    	UpdateColors(ips, colold, col);
    }

    return;
}

Int_t GetNextColor()
{
	static int num=0;
	num++;
	Int_t color = gStyle->GetColorPalette(num); 	
	return color; 
}

Int_t GetJetColor(TH2F* jethist)
{
	// first, set some default color
	static int num=0;
	Int_t color = kGreen+num; 

	// now, let's try to get jet histograms from the file provided
	//std::cout << "file = " << _jet_cols_fname << std::endl;
	std::vector< std::vector<int> > allcols = GetColours();
	TString dir = TString::Format("JetDisplay_event%02d", _event);
	TFile* file = new TFile(_jet_cols_fname, "READ");
	std::vector<TH2F*> jets;
	if (file && file->IsOpen()){
		jets = GetJets(file, dir);
	} else {
		std::cout << "GetJetColor: Couldn't open file with jet display, returning standard jet color. " << std::endl;
		num++; // increase num by one, because we have used this color
		return color;
	}

	// now compare each of these jets with the jet we want to color
	int Nmatch = -1; 
	int imatch = -1;
	for (int i=0; i<jets.size(); ++i){
		Nmatch = CompareHists(jethist, jets[i]);
		if (Nmatch>2){
			imatch = i; 
			break;
		}
	}

	// success if NMatch > 0
	if (Nmatch>0){
		if (imatch<=2) color = allcols[imatch][0];
		else color = allcols[2][0];
		std::cout << "GetJetColor: Found a jet matching to this cluster. Returning color: " << color << std::endl;
		return color;
	}

	// still here? then try jets rejected by SD criterion
	TH2F* subjets_rejected_sd = GetSubjetsRejectedBySD(file, dir);
	Nmatch = CompareHists(jethist, subjets_rejected_sd);
	if (Nmatch>0){
		color = TColor::GetColor( "#808080"); // dark gray
		std::cout << "GetJetColor: Found a SD-rejected cluster matching to this cluster. Returning color #808080: " << color << std::endl;		
		return color; 
	}

	// still here? now try jets rejected because they don't have a subjet
    TH2F* jets_rejected = GetJetsWithoutSubjets(file, dir);
	Nmatch = CompareHists(jethist, jets_rejected);
	if (Nmatch>0){
		color = TColor::GetColor( "#e6e6e6"); // light gray
		std::cout << "GetJetColor: Found a jet without subjets matching to this cluster. Returning color #e6e6e6: " << color << std::endl;		
		return color; 
	}

	// lastly: all jets with pt smaller than ptsub and ptjet
	TH2F* jets_low_pt = GetJetsWithLowPt(file, dir);
	Nmatch = CompareHists(jethist, jets_low_pt);
	if (Nmatch>0){
		color = TColor::GetColor( "#c4edc4"); // very light green
		std::cout << "GetJetColor: Found a low-pt (sub)jet matching to this cluster. Returning color #c4edc4: " << color << std::endl;		
		return color; 
	}	

	// still here? then there is no match to any of the jets in the display
	// maybe the histogram bin sizes don't match, or something else went wrong.
	// in this case, return a standard color
	num++; // increase num by one, because we have used this color
	std::cout << "GetJetColor: Did not find any jet in the jet display histograms. Returning default color: " << color << std::endl;		
	return color;

}

int CompareHists(TH2F* h1, TH2F* h2)
{
	int Nmatch = 0; 
	for (int ix=1; ix<h1->GetNbinsX()+1; ++ix){
		for (int iy=1; iy<h1->GetNbinsY()+1; ++iy){
			if (h1->GetBinContent(ix,iy)>0){
				double x = h1->GetXaxis()->GetBinCenter(ix);
				double y = h1->GetYaxis()->GetBinCenter(iy);
				//std::cout << "found content in bin: " << ix << ", " << iy << " with x = " << x << ", y = " << y << std::endl;
				int jx = h2->GetXaxis()->FindBin(x);
				int jy = h2->GetYaxis()->FindBin(y);				
				if (h2->GetBinContent(jx,jy) > 0) Nmatch++;
			}
		}
	}

	return Nmatch; 
}

void DrawAll(std::vector<part> ips, std::map<int, TH2F*> chists)
{	
	// draws the input particles and clusters on the current pad
	
	// first: draw all input particles
	PlotInputParts(ips); 

	// now draw the clusters by looping over all stored histograms
	for (const auto& [key, cluster] : chists){
      //std::cout << "Drawing cluster with index " << '[' << key << "] = " << cluster->GetName() << "; Entries = " << cluster->GetEntries() 
      //          << " color = " << cluster->GetFillColor() <<  std::endl;
      DrawContour(cluster, cluster->GetFillColor()); 
    }

}

void PrintStep(int step, std::vector<part> ips, std::vector<leave> tree)
{
  Int_t c_Gray = TColor::GetColor( "#8c8c8c" );
  uint Np = ips.size();
  int ch = step - 1 + Np; 

  TString info = TString::Format("Clustering step: %04d, d_{ij} = %f", ch, tree[ch].dij);
  TLatex *text1 = new TLatex(3.5, 24, info);
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.64);
  text1->SetTextFont(62);
  text1->SetTextSize(0.02);
  text1->SetY(0.985);
  //text1->SetTextColor(c_Gray);
  text1->SetTextColor(kBlack);
  text1->Draw();

}













