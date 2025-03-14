#include <Riostream.h>
#include <filesystem>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TSystem.h"
#include "TDatabasePDG.h"
#include <string>
#include "TDirectory.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "../library/dimuVarsCommon.h"
#include "TNamed.h"

//global values: kine cuts
float lowPt = 0;
float upPt = 0.2;
float lowMass = 2.9;
float upMass = 3.3;
float lowRap = -4;
float upRap = -2.5;

// custon rew function
double pcFunc(double *x, double *par){
  if(x[0]<0) return par[0]*(1 + par[1]*TMath::Sin(x[0]));
  else       return par[0]*(1 + par[2]*TMath::Sin(x[0]));
}

// gloabla pointer to the weights
TH1D *hWeights = NULL;

// global pointer to the rew function
TF1 *gRewFuc = NULL; 

// function for the weight
double rewFunction(double phi){
  return gRewFuc->Eval(phi);
}

// function that returns the weight from the histogram
double WeightFunction(double* x, double*) {
  int bin = hWeights->FindBin(x[0]);
  return hWeights->GetBinContent(bin);
}

// do some cosmetcis to TH1
void cosmetics(TH1D *h, string dataType){
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.5);
  
  if(dataType == "recoMC"){
    h->SetMarkerColor(kGreen+1);
    h->SetLineColor(kGreen+1);
  }
  
  h->GetYaxis()->SetTitle("#counts");
}


// entry point
void produceRewHistos(string dataType = "recoMC", bool applyKine = true, bool useWeightsFromFile = true){

  // read the rew function obtained from the fit to the phi distrib. of data
  if(useWeightsFromFile){
    TFile *wFile = new TFile("hWeightsPtTracks.root");
    hWeights = (TH1D*)wFile->Get("hWeights");

    // Create a TF1 that wraps the WeightFunction function.
    gRewFuc = new TF1("gRewFuc", WeightFunction, hWeights->GetXaxis()->GetXmin(),
                            hWeights->GetXaxis()->GetXmax(), 0);
  }
  
  // exit if there is no fuction
  if(gRewFuc==NULL){
    cout<<"There is not a rew function. Bye!"<<endl;
    return;
  }

  gStyle->SetOptStat(1);
  string filePath = "";
  string treePath = "";
  string outputFile = "";

  if(dataType=="recoMC"){
    filePath = "../MonteCarlo/reco_tree.root";
    treePath = "DF_2336518085565631/dimu"; 
    outputFile = "controlRewRecoMC.root";
  }

  // opening the file that stores the tree 
  TFile *file = new TFile(filePath.c_str());
  // create the tree
  if(dataType=="recoMC") createRecoTree(file,treePath.c_str());

  //DEFINING THE HISTOGRAM OF THE CONTROL PLOTS
  char htitle[100]; 
  char hname[100];

  // 1) track distributions: pt, eta, phi
  TH1D *trackPtDistrib[2]  = {};
  TH1D *trackEtaDistrib[2] = {}; 
  TH1D *trackPhiDistrib[2] = {};
  // track pt below and above the bump in pair phi
  TH1D *trackPtLowPhiDistrib[2]  = {};
  TH1D *trackPtHighPhiDistrib[2]  = {};

  string id[2] = {"positive","negative"};
  for(int i=1; i<=2; i++){
    //pt
    sprintf(hname,"trackPtDistrib%d",i);
    sprintf(htitle,"Distribution of the transverse momentum of the %s tracks - %s",id[i-1].c_str(), dataType.c_str());
    trackPtDistrib[i-1] = new TH1D(hname,htitle,100,0.7,2.5);
    trackPtDistrib[i-1]->GetXaxis()->SetTitle("#it{p}_{T, trk} (GeV/#it{c})");
    //eta
    sprintf(hname,"trackEtaDistrib%d",i);
    sprintf(htitle,"Distribution of the pseudorapidity of the %s tracks - %s",id[i-1].c_str(), dataType.c_str());
    trackEtaDistrib[i-1] = new TH1D(hname,htitle,100,-4,-2.5);
    trackEtaDistrib[i-1]->GetXaxis()->SetTitle("#eta_{trk}");
    //phi   
    sprintf(hname,"trackPhiDistrib%d",i);
    sprintf(htitle,"Distribution of the #varphi of the %s tracks - %s",id[i-1].c_str(), dataType.c_str());
    trackPhiDistrib[i-1] = new TH1D(hname,htitle,100,-TMath::Pi(),TMath::Pi());
    trackPhiDistrib[i-1]->GetXaxis()->SetTitle("#varphi_{trk}");

    //pt
    sprintf(hname,"trackPtLowPhiDistrib%d",i);
    sprintf(htitle,"Distribution of the transverse momentum of the %s tracks, #phi < 0.1 - %s",id[i-1].c_str(), dataType.c_str());
    trackPtLowPhiDistrib[i-1] = new TH1D(hname,htitle,100,0.7,2.5);
    trackPtLowPhiDistrib[i-1]->GetXaxis()->SetTitle("#it{p}_{T, trk} (GeV/#it{c})");

    sprintf(hname,"trackPtHighPhiDistrib%d",i);
    sprintf(htitle,"Distribution of the transverse momentum of the %s tracks, #phi > 0.1 - %s",id[i-1].c_str(), dataType.c_str());
    trackPtHighPhiDistrib[i-1] = new TH1D(hname,htitle,100,0.7,2.5);
    trackPtHighPhiDistrib[i-1]->GetXaxis()->SetTitle("#it{p}_{T, trk} (GeV/#it{c})");
  }

  // 2) pair distributions: energy, mass, pt, phi, rapidity
  TH1D *hMass = new TH1D("hMass", Form("Invariant mass distribution - %s", dataType.c_str()), 100, lowMass, upMass);
  hMass->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/#it{c}^{2})");
  TH1D *hPt = new TH1D("hPt", Form("Transverse momentum distribution - %s", dataType.c_str()), 100, lowPt, upPt);
  hPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  TH1D *hPhi = new TH1D("hPhi", Form("#varphi distribution - %s", dataType.c_str()), 100, -TMath::Pi(), TMath::Pi());
  hPhi->GetXaxis()->SetTitle("#varphi");
  TH1D *hRapidity = new TH1D("hRapidity", Form("Rapidity distribution - %s", dataType.c_str()), 100, lowRap, upRap);
  hRapidity->GetXaxis()->SetTitle("y");

  // 3) azimuth angle: average, charge
  TH1D *hPhiAverage = new TH1D("hPhiAverage", Form("#phi #it{average} distribution - %s", dataType.c_str()), 100, -TMath::Pi(), TMath::Pi());
  hPhiAverage->GetXaxis()->SetTitle("#phi #it{average}");
  TH1D *hPhiCharge = new TH1D("hPhiCharge", Form("#phi #it{charge} distribution - %s", dataType.c_str()), 100, -TMath::Pi(), TMath::Pi());
  hPhiCharge->GetXaxis()->SetTitle("#phi #it{charge}");

  // fill the histos
  for(Long64_t ev=0; ev<globalTree->GetEntries(); ev++){
    globalTree->GetEvent(ev);

    // aply kinematics cuts if requested
    if(applyKine){

      if(fPt < lowPt) continue;
      if(fPt > upPt) continue;

      if(fM < lowMass) continue;
      if(fM > upMass) continue;

      if(fRap < lowRap) continue;
      if(fRap > upRap) continue;
    }
    
    // trk distribs
    trackPtDistrib[0]->Fill(fPtp, rewFunction(fGenPhi));
    trackEtaDistrib[0]->Fill(fEtap, rewFunction(fGenPhi));
    if(fPhip<TMath::Pi()) trackPhiDistrib[0]->Fill(fPhip, rewFunction(fGenPhi));
    if(fPhip>TMath::Pi()) trackPhiDistrib[0]->Fill(fPhip-2*TMath::Pi(), rewFunction(fGenPhi));

    trackPtDistrib[1]->Fill(fPtn, rewFunction(fGenPhi));
    trackEtaDistrib[1]->Fill(fEtan, rewFunction(fGenPhi));
    if(fPhin<TMath::Pi()) trackPhiDistrib[1]->Fill(fPhin, rewFunction(fGenPhi));
    if(fPhin>TMath::Pi()) trackPhiDistrib[1]->Fill(fPhin-2*TMath::Pi(), rewFunction(fGenPhi));

    if(fPhi < 0.1) trackPtLowPhiDistrib[0]->Fill(fPtp, rewFunction(fGenPhi));
    if(fPhi < 0.1) trackPtLowPhiDistrib[1]->Fill(fPtp, rewFunction(fGenPhi));
    if(fPhi > 0.1) trackPtHighPhiDistrib[0]->Fill(fPtn, rewFunction(fGenPhi));
    if(fPhi > 0.1) trackPtHighPhiDistrib[1]->Fill(fPtn, rewFunction(fGenPhi));

    // pair distribs
    hMass->Fill(fM, rewFunction(fGenPhi));
    hPt->Fill(fPt, rewFunction(fGenPhi));
    hPhi->Fill(fPhi, rewFunction(fGenPhi));
    hRapidity->Fill(fRap, rewFunction(fGenPhi));

    hPhiAverage->Fill(fPhiAv, rewFunction(fGenPhi));
    hPhiCharge->Fill(fPhiCh, rewFunction(fGenPhi));
  }

  // aplly cosmetics
  cosmetics(trackPtDistrib[0], dataType);
  cosmetics(trackEtaDistrib[0], dataType);
  cosmetics(trackPhiDistrib[0], dataType);
  cosmetics(trackPtDistrib[1], dataType);
  cosmetics(trackEtaDistrib[1], dataType);
  cosmetics(trackPhiDistrib[1], dataType);
  cosmetics(trackPtLowPhiDistrib[0], dataType);
  cosmetics(trackPtLowPhiDistrib[1], dataType);
  cosmetics(trackPtHighPhiDistrib[0], dataType);
  cosmetics(trackPtHighPhiDistrib[1], dataType);
  cosmetics(hMass, dataType);
  cosmetics(hPt, dataType);
  cosmetics(hPhi, dataType);
  cosmetics(hRapidity, dataType);
  cosmetics(hPhiAverage, dataType);
  cosmetics(hPhiCharge, dataType);

    //save the results in a root file
    TFile *fRes = new TFile(outputFile.c_str(),"recreate");

    trackPtDistrib[0]->Write();
    trackEtaDistrib[0]->Write();
    trackPhiDistrib[0]->Write();
    trackPtDistrib[1]->Write();
    trackEtaDistrib[1]->Write();
    trackPhiDistrib[1]->Write();
    trackPtLowPhiDistrib[0]->Write();
    trackPtLowPhiDistrib[1]->Write();
    trackPtHighPhiDistrib[0]->Write();
    trackPtHighPhiDistrib[1]->Write();
    hMass->Write();
    hPt->Write();
    hPhi->Write();
    hRapidity->Write();
    hPhiAverage->Write();
    hPhiCharge->Write();
  
    fRes->Close();
}