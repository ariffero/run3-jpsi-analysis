//
// it computes the resolution for angular variables, such as phi, phiAv, etc...
// plus some other distributions to understand why the resolution is as it is
//

// C++ headers
#include <iostream>
#include <string>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

// my headers
#include "../library/dimuVarsCommon.h"

// global varaibles
const float massMuon = 0.105658; // gev

// -------------
// entry point
void angularResolution(string var = "phi"){

  // string used to get the tree
  gStyle->SetOptStat(1);
  string filePath = "";
  string treePath = "";

  filePath = "../MonteCarlo/reco_tree.root";
  treePath = "DF_2336518085565631/dimu"; 

  // histo for the resolution
  TH1D *hResolution = new TH1D("hResolution",Form("%s resolution",var.c_str()),1000,-TMath::Pi(),TMath::Pi());
  hResolution->GetYaxis()->SetTitle("#counts");

  // histogram for the correlation gen-reco
  TH2D *hCorrelation = new TH2D("hCorrelation","hCorrelation",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
  hCorrelation->GetXaxis()->SetTitle("reco");
  hCorrelation->GetYaxis()->SetTitle("gen");

  // histo for the correlation of pos/neg
  TH2D *hPosNegCorr = new TH2D("hPosNegCorr","hPosNegCorr",5000,-TMath::Pi(),TMath::Pi(),5000,-TMath::Pi(),TMath::Pi());
  hPosNegCorr->GetXaxis()->SetTitle("#varphi_{trk pos, gen} - #varphi_{trk pos, reco}");
  hPosNegCorr->GetYaxis()->SetTitle("#varphi_{trk neg, gen} - #varphi_{trk neg, reco}");

  // histo for resultion vs pT
  TH2D *hResVsPt = new TH2D("hResVsPt","hResVsPt",500,-TMath::Pi(),TMath::Pi(),500,0,5);
  hResVsPt->GetXaxis()->SetTitle("#phi_{gen} - #phi_{reco}");
  hResVsPt->GetYaxis()->SetTitle("reco #it{p}_{T} (GeV/#it{c})");
  // vs gen pT
  TH2D *hResVsPt2 = (TH2D*)hResVsPt->Clone("hResVsPt");
  hResVsPt2->GetYaxis()->SetTitle("gen #it{p}_{T} (GeV/#it{c})");

  double sum = 1;
  double dif = 1;
  double res = 1;

  // opening the file that stores the tree 
  TFile *file = new TFile(filePath.c_str());
  createRecoTree(file,treePath.c_str());

  // define TLorentzVectors
  TLorentzVector vGen;
  TLorentzVector vGenDiff;
  TLorentzVector vGenDiffOpp;
  TLorentzVector vGenP;
  TLorentzVector vGenN;
  TLorentzVector vRec;
  TLorentzVector vRecDiff;
  TLorentzVector vRecDiffOpp;
  TLorentzVector vRecP;
  TLorentzVector vRecN;

  // fill the histos
  for(Long64_t ev=0; ev<globalTree->GetEntries(); ev++){
    globalTree->GetEvent(ev);

    // fill the TLorentzVector of reco events
    vRecP.SetPtEtaPhiM(fPtp,fEtap,fPhip,massMuon);
    vRecN.SetPtEtaPhiM(fPtn,fEtan,fPhin,massMuon);
    vRec = vRecP + vRecN;
    vRecDiff = vRecP - vRecN;
    // fill the TLorentzVector of gen events
    vGenP.SetPtEtaPhiM(fGenPtp,fGenEtap,fGenPhip,massMuon);
    vGenN.SetPtEtaPhiM(fGenPtn,fGenEtan,fGenPhin,massMuon);
    vGen = vGenP + vGenN;
    vGenDiff = vGenP - vGenN;

    // compute the resolution for different variables
    // phi of the j/psi
    if(var=="phiJpsi"){
      sum = 1;
      dif = vGen.DeltaPhi(vRec);
      hResolution->GetXaxis()->SetTitle("#varphi_{gen} - #varphi_{reco}");
      hResolution->SetTitle("Resolution of #varphi of J/#psi");

      hCorrelation->Fill(vRec.Phi(),vGen.Phi());
    }
    // phi of the state formed by the difference of 4-vectors
    else if(var=="phiDiff"){
      sum = 1;
      dif = vGenDiff.DeltaPhi(vRecDiff);
      hResolution->GetXaxis()->SetTitle("#varphi_{gen} - #varphi_{reco}");
      hResolution->SetTitle("Resolution of #varphi of (#mu^{+}-#mu^{-}) state");

      hCorrelation->Fill(vRecDiff.Phi(),vGenDiff.Phi());
    }
    // phi of positive tracks
    else if(var=="phiTrkPos"){
      sum = 1;
      dif = vGenP.DeltaPhi(vRecP);
      hResolution->GetXaxis()->SetTitle("#varphi_{trk pos, gen} - #varphi_{trk pos, reco}");
      hResolution->SetTitle("Resolution of positive tracks #varphi");

      // fill the histo with the pos neg correlation
      hPosNegCorr->Fill(vGenP.DeltaPhi(vRecP), vGenN.DeltaPhi(vRecN));
    }
    // phi of negative tracks
    else if(var=="phiTrkNeg"){
      sum = 1;
      dif = vGenN.DeltaPhi(vRecN);
      hResolution->GetXaxis()->SetTitle("#varphi_{trk neg, gen} - #varphi_{trk neg, reco}");
      hResolution->SetTitle("Resolution of negative tracks #varphi");

      // fill the histo with the pos neg correlation
      hPosNegCorr->Fill(vGenP.DeltaPhi(vRecP), vGenN.DeltaPhi(vRecN));
    }
    // phi average
    else if(var=="phiAv"){
      sum = 1;
      double genPhiAv;
      double phiAv;
      if(gRandom->Rndm()>0.5){
        genPhiAv = vGen.DeltaPhi(vGenDiff);
        phiAv = vRec.DeltaPhi(vRecDiff);
      }
      else{
        genPhiAv = vGen.DeltaPhi(vGenDiffOpp);
        phiAv = vRec.DeltaPhi(vRecDiffOpp);
      }
      dif = genPhiAv - phiAv;
      hResolution->GetXaxis()->SetTitle("#phi_{average, gen} - #phi_{average reco}");
      hResolution->SetTitle("Resolution of #phi_{average}");
    }
    // phi charge
    else if(var=="phiCh"){
      sum = 1;
      double genPhiCh = vGen.DeltaPhi(vGenDiff);
      dif = genPhiCh - fPhiCh;
      hResolution->GetXaxis()->SetTitle("#phi_{charge, gen} - #phi_{charge reco}");
      hResolution->SetTitle("Resolution of #phi_{charge}");
    }

    // fill the histos
    res = (dif/sum);
    hResolution->Fill(res);
    hResVsPt->Fill(res,fPt);
    hResVsPt2->Fill(res,fGenPt);
  }  

  // draw the plot with the resolution
  TCanvas *c = new TCanvas("c","c",1920,1080);
  hResolution->Draw();

  // draw the histo with the correlation between the difference of reco and gen phi_trk
  // in positive and negative tracks (if not empty)
  if(hPosNegCorr->GetEntries()!=0){
    TCanvas *cpn = new TCanvas("cpn","cpn",1920,1080);
    hPosNegCorr->GetXaxis()->SetRangeUser(-0.1,0.1);
    hPosNegCorr->GetYaxis()->SetRangeUser(-0.1,0.1);
    hPosNegCorr->Draw("colz");
  }

  // draw the gen-reco correlation plot if it not empty
  if(hCorrelation->GetEntries()!=0){
    TCanvas *c2 = new TCanvas("c2","c2",1920,1080);
    hCorrelation->Draw("colz");
  }

  // draw the resolution of phi vs pT if requested
  if(var=="phiJpsi"){
    TCanvas *cPt = new TCanvas("cPt","cPt");
    hResVsPt->Draw("colz");  
  }

}