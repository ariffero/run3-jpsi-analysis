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

//global values: kine cuts
float lowPt = 0;
float upPt = 0.2;
float lowMass = 2.2;
float upMass = 4.0;
float lowRap = -4;
float upRap = -2.5;


// fit function
double func(double *x, double *par){
  return par[0]*(1 + par[1]*TMath::Sin(x[0]) + par[2]*TMath::Cos(x[0]));
}

// fit function
double pcFunc(double *x, double *par){
  if(x[0]<0) return par[0]*(1 + par[1]*TMath::Sin(x[0]));
  else       return par[0]*(1 + par[2]*TMath::Sin(x[0]));
}



// entry point
void fitPhi(int reg = 3, bool applyKine = true){

  string filePath = "";
  string treePath = "";

  filePath = "../Data/merged-data.root";
  treePath = "DF_2336518079279520/O2dimu";

  TH1D *hPhi = new TH1D("hPhi","#varphi distribution", 100, -TMath::Pi(), TMath::Pi());
  hPhi->GetXaxis()->SetTitle("#varphi");

  // opening the file that stores the tree 
  TFile *file = new TFile(filePath.c_str());
  // create the tree
  createDataTree(file,treePath.c_str());

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

    hPhi->Fill(fPhi);
  }

  TF1 *fitFunc = NULL;
  if(reg==1) fitFunc = new TF1("fitFunc",func,-TMath::Pi(),0,2);
  if(reg==2) fitFunc = new TF1("fitFunc",func,0,TMath::Pi(),2);
  if(reg==3) fitFunc = new TF1("fitFunc",pcFunc,-TMath::Pi(),TMath::Pi(),3);
  fitFunc->SetParNames("norm","amp1","amp2");
  fitFunc->SetParameter(1,0.2);
  fitFunc->SetParameter(2,0.4);

  hPhi->Fit(fitFunc,"RMS+");

  gStyle->SetOptFit(1);
  TCanvas *c = new TCanvas();
  hPhi->Draw("ep");

  TFile *saveFile = new TFile("rewFunc.root","recreate");
  // normalize the function, computing the normalization by hand
  double norm = 2*TMath::Pi() + 2*(fitFunc->GetParameter(2) - fitFunc->GetParameter(1));
  fitFunc->SetParameter(0,1/norm);
  // check that the function is correctly normalized
  cout<<"integral: "<<fitFunc->Integral(-TMath::Pi(),TMath::Pi())<<endl;
  fitFunc->Write();
  saveFile->Close();
}