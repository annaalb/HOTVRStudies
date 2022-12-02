#include <TString.h>
#include <TSystem.h>
#include <algorithm>

#include "style.h"
#include "jet_display_helpers.h"
#include "hist_helpers.h"


void do_histories(TString prefix);

//TString gdir = "NoCutsHistosGen_mass/"; // use this for HOTVR display
//TString gdir = "NoCutsHistos3Gen/";       // anti-kt jets

TString _jet_cols_fname;
uint _event;

int main()
{
    gROOT->SetBatch(true);
    //gStyle->SetPalette(kLightTemperature);
    gStyle->SetPalette(kBlueRedYellow);
    TColor::InvertPalette();

    std::vector<TString> prefixes; 
    //prefixes.push_back(TString("history_rho"));
    //prefixes.push_back(TString("history_ET"));
    prefixes.push_back(TString("history_ghosts"));

    for (int i=0; i<prefixes.size(); ++i){
      do_histories(prefixes[i]);
    }

    exit(0);
}

void do_histories(TString prefix)
{
    gROOT->SetBatch(true);

    _event = 0;

    TString fname_tree = TString("files/") + prefix + TString("/history_") + TString::Format("%d",_event) + TString(".txt");
    TString fname_input = TString("files/") + prefix + TString("/pseudojets_input_") + TString::Format("%d",_event) + TString(".txt");

    // here is the file for the display (optional)
    TString hpre = "mytests2";
    TString fname_display = TString("files/") + prefix + TString("/uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23.root");
    SetJetColsFile(fname_display); 
    SetJetColsEvent(_event);

    std::cout << "getting info from clustering tree: " << fname_tree << std::endl;
    std::cout << "input particles: " << fname_input << std::endl;
    std::cout << "get the jet colors from matching jets in this file: " << fname_display << std::endl;

    TString plot_dir = TString("plots/") + prefix;
    TString cmd = TString("mkdir -p ") + plot_dir; 
    gSystem->Exec(cmd.Data()); 
    cmd.Append("/steps");
    gSystem->Exec(cmd.Data()); 
    TString outname = plot_dir + "/steps" + "/history" + "_hotvr"; 

    std::vector<part> inputs = GetInputParticles(fname_input);
    //for (int i=0; i<inputs.size(); ++i){
    //  std::cout << "ind = " << inputs[i].ind << " pt = " << inputs[i].pt << std::endl;
    //}

    std::vector<leave> tree = GetHistoryTree(fname_tree); 
    //for (int i=0; i<tree.size(); ++i){
    //  std::cout << "ind = " << tree[i].ind << " parent1 = " << tree[i].parent1 << std::endl;
    //}

    uint Np = inputs.size();
    std::cout << "Size of inputs = " << Np << std::endl;

    SetupGlobalStyle(); 
    TString canname = "canvas"; // TString::Format("canvas_%d", num);
    TCanvas *c = SetupCanvas2d(canname);
    PlotAxes(); 

    // first: plot all input particles (maybe with different sizes of dots)
    // and safe as 0th step
    PlotInputParts(inputs); 
    PrintStep(0, inputs, tree);
    gPad->RedrawAxis();
    TString outname_step = outname + "_00000";      
    c->SaveAs(outname_step + ".gif");
    c->Clear();

    // now go through history and plot it, safe every step as single graphic
    std::map<int, TH2F*> chists;
    uint Nsteps = tree.size() - Np + 1;
    //Nsteps = 814; 
    for (uint i=1; i<Nsteps; ++i){
      bool draw_step; 
      draw_step = DoClusteringStep(i, inputs, tree, chists);       
      
      std::cout << "size of map with hists = " << chists.size() << "   draw_step = " << draw_step << std::endl;
      //for (const auto& [key, value] : chists){
      //  std::cout << '[' << key << "] = " << value->GetName() << "; Entries = " << value->GetEntries() << " color = " << value->GetFillColor() << std::endl;
      //}
      //for (uint j=0; j<chists.size(); ++j){
      //  cout << "UID = " << chists[j] << ", name of hist = " << chists[j]->GetName() << endl;
      //}

      //draw_step = false;
      //if (i<1516) draw_step = false;
      //if (i>3177) draw_step = true;
      if (draw_step){
        PlotAxes();             
        DrawAll(inputs, chists);
        PrintStep(i, inputs, tree);
        gPad->RedrawAxis();
        TString outname_step = TString::Format(outname + "_%05d", i);      
        c->SaveAs(outname_step + ".gif");  
        c->Clear();
      }

      //if (i==280) exit(0);

      //c->SaveAs(outname_step + ".");  
      /*
      if (i==1 || i%10==0){
        TString pdfname = TString::Format(plot_dir + "/history_hotvr" + "_%d", i);
        c->SaveAs(pdfname + ".pdf");
      }
      */
    }

    // last picture: plot all clusters and on top the input particles in gray
    PlotAxes();             
    DrawAll(inputs, chists);
    PlotInputParts(inputs, kBlack);
    gPad->RedrawAxis();
    outname_step = outname + TString("_final");      
    c->SaveAs(outname_step + ".gif");  
    c->Clear();


}

