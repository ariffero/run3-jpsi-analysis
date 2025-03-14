//
// macro used ot study the reweighting of the pT distribution of tracks
//

// C++ headers
#include <iostream>
#include <string>
#include <vector>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TStyle.h"

// Custom headers: used to read the trees
#include "../library/dimuVarsCommon.h"

//global values: kine cuts
float lowPt = 0;
float upPt = 0.25;
float lowMass = 2.9;
float upMass = 3.3;
float lowRap = -4;
float upRap = -2.5;

bool applyKine = true;

// -------------------------------------------------------------
// function to fill the vector passed with the values in the tree
void fillVector(string filePath, string treePath, string dataNature, vector<double>& vec, string varName ="ptTracks"){

  // open the file that stores the tree 
  TFile *file = new TFile(filePath.c_str());
  // call the function that allow to access the TTree accordingly to the tipe of the data in the file
  if(dataNature=="data")  createDataTree(file,treePath.c_str());
  if(dataNature=="reco")  createRecoTree(file,treePath.c_str());

  for(Long64_t ev=0; ev<globalTree->GetEntries(); ev++){
    globalTree->GetEvent(ev);

    if(applyKine){
      if(fPt < lowPt) continue;
      if(fPt > upPt) continue;

      if(fM < lowMass) continue;
      if(fM > upMass) continue;

      if(fRap < lowRap) continue;
      if(fRap > upRap) continue;
    }

    if(varName=="ptTracks")   vec.push_back(fPtn);
    if(varName=="ptTracks")   vec.push_back(fPtp);

    if(varName=="genPtTracks")   vec.push_back(fGenPtn);
    if(varName=="genPtTracks")   vec.push_back(fGenPtp);

  }
}

// -------------------------------------------------------------------------------
// function to get the y corresponding the the passed x of the passed histogram 
double getYfromX(TH1 *h, double x){
  double y = h->GetBinContent(h->FindBin(x));
  return y;
}

// function to fill a TH1 with weighting on another variable
void fillTH1OneWeight(TH1 *h, vector<double> v,vector <double> vec_weight, TH1 *hWeights, string xTitle = "#it{p}_{T} (GeV/#it{c})"){
  double w = 1;
  for(unsigned int i=0; i<v.size(); i++){
      w = getYfromX(hWeights,vec_weight[i]);
      h->Fill(v[i],w);
  }

  h->GetXaxis()->SetTitle(xTitle.c_str());
  h->GetYaxis()->SetTitle("#counts");
}

// -------------------------------------------
// function to normalize a TH1
void normalize(TH1 *h){
  h->Scale(1/h->Integral());
}

// -------------------------------------------
// entry point: main function
void studyPtTracksRew(){
  
  // strings to read the trees
  string filePath = "";
  string treePath = "";

  // for the data
  filePath = "../Data/merged-data.root";
  treePath = "DF_2336518079279520/O2dimu";

  // vectors for the pT of the tracks
  vector<double> vec_ptTracksGen;
  vector<double> vec_ptTracksReco;
  vector<double> vec_ptTracksData;

  // fill the vector with the pT of tracks in data
  fillVector(filePath, treePath,"data",vec_ptTracksData,"ptTracks");

  // for reco MC
  filePath = "../MonteCarlo/reco_tree.root";
  treePath = "DF_2336518085565631/dimu"; 

  fillVector(filePath, treePath,"reco",vec_ptTracksReco,"ptTracks");
  fillVector(filePath, treePath,"reco",vec_ptTracksGen,"genPtTracks");

  cout<<vec_ptTracksGen.size()<<" "<<vec_ptTracksReco.size()<<" "<<vec_ptTracksData.size()<<endl;

  // declare the histos with the pT of tracks
  double lowPtTrk = 0.7;
  double upPtTrk = 2.5;
  TH1D *hPtTrksData = new TH1D("hPtTrksData","hPtTrksData",100,lowPtTrk,upPtTrk);
  TH1D *hPtTrksReco = new TH1D("hPtTrksReco","hPtTrksReco",100,lowPtTrk,upPtTrk);
  TH1D *hPtTrksGen = new TH1D("hPtTrksGen","hPtTrksGen",100,lowPtTrk,upPtTrk);
  // rew histo
  TH1D *hPtTrksRecoWeighted = (TH1D*)hPtTrksReco->Clone("hPtTrksRecoWeighted");
  hPtTrksRecoWeighted->SetTitle("hPtTrksRecoWeighted");
  // histo with the weights
  TH1D *hWeights = new TH1D("hWeights","hWeights",100,lowPtTrk,upPtTrk);

  // fill the histos
  for(unsigned int i=0; i<vec_ptTracksData.size(); i++){
    hPtTrksData->Fill(vec_ptTracksData[i]);
  }
  for(unsigned int i=0; i<vec_ptTracksReco.size(); i++){
    hPtTrksReco->Fill(vec_ptTracksReco[i]);
    hPtTrksGen->Fill(vec_ptTracksGen[i]);
  }

  // compute and save the weights
  hWeights->Divide(hPtTrksData,hPtTrksReco);
  hWeights->SaveAs("hWeightsPtTracks.root");

  // apply the weights
  fillTH1OneWeight(hPtTrksRecoWeighted,vec_ptTracksReco,vec_ptTracksGen,hWeights);

  // do some cosmetics and then check the results

  hPtTrksData->SetLineColor(kRed);
  hPtTrksReco->SetLineColor(kBlue);
  hPtTrksRecoWeighted->SetLineColor(kBlack);

  hPtTrksData->SetLineWidth(2);
  hPtTrksReco->SetLineWidth(2);
  hPtTrksRecoWeighted->SetLineWidth(2);

  normalize(hPtTrksData);
  normalize(hPtTrksReco);
  normalize(hPtTrksRecoWeighted);

  hPtTrksReco->Divide(hPtTrksData);
  hPtTrksRecoWeighted->Divide(hPtTrksData);

  TCanvas *c = new TCanvas();

  hPtTrksReco->Draw("histo same");
  hPtTrksData->Draw("histo same");
  hPtTrksRecoWeighted->Draw("histo same");

  TLegend *leg = new TLegend();
  leg->AddEntry(hPtTrksData,"data");
  leg->AddEntry(hPtTrksReco,"reco");
  leg->AddEntry(hPtTrksRecoWeighted,"reco rew");
  leg->Draw();

}