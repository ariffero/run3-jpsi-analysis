//
// makes a fit to the invariant mass distribution of muon pairs (at fwd rapidity)
// fit function is J/psi + psi(2s) (optional)
// it is possible to choose the number of bins in phi and the neutron emission class
//

// -----------------------------------------------------------------
// all headers are defined here
// -----------------------------------------------------------------

// c++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"

// RooFit headers
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooAbsReal.h"
#include "RooGenericPdf.h"
#include "RooCrystalBall.h"
#include "RooFormulaVar.h"
#include "RooHist.h"

using namespace RooFit;

// My headers: needed to save the results in the tree
#include "../library/savedVarInMassFits.h"

// -----------------------------------------------------------------
// global variables, to be set to drive the fits
// -----------------------------------------------------------------
// --> kinematics: will be set to the values in a config
double minPt = 0;
double maxPt = 0;
double minMass = 0;
double maxMass = 0;
double minRapidity = 0;
double maxRapidity = 0;

// --> Crystall ball
double nParL = 0;
double nParR = 0;
double sigmaL = 0;
double sigmaR = 0;
double alphaL = 0;
double alphaR = 0;

// --> exponential
double lambda = 0.0;

// --> switches 
bool isNFixed = false;
bool isSigmaFixed = false;
bool isAlphaFixed = false;
bool includePsi2s = true;
bool isLambdaFixed = false;
bool isChi2Fit = false;
bool excludeJPsi = false;

// position of the legend info
float gXpos = 0.17;
float gYpos = 0.85;

// -----------------------------------------------------------------
// global variables useful in the macro
// -----------------------------------------------------------------
//global variable for number of phi bins
int gPhiBins = 1;
//global variable for neutron classes
string gNeutronClass = "";
//name of the file with the results
string gSaveFileName = "";
//name of the tree with the results
string gSaveTreeName = "";
//name of the (1st) folder that contains the results
string gIdentifier = "";

// -----------------------------------------------------------------
// all functions are defined here
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// do one data fit
void doOneDataFit(TTree *dataTree, TFile *saveFile, float binID[3])
{
  // number of times the function has been called
  static int nCalls = 0;

  // define variables for the tree
  RooRealVar pt("fPt","pt",minPt,maxPt);
  RooRealVar mass("fM","m_{#mu#mu} (GeV/c^{2})",minMass,maxMass);
  RooRealVar rapidity("fRap","rapidity",minRapidity,maxRapidity);
	RooRealVar phiAverage("fPhiAv","phiAverage",binID[1],binID[2]);

  // neutron class from neutron ID
  // 0n0n == 1
  // Xn0n == 2
  // 0nXn == 3
  // XnXn == 4
  RooRealVar *neutronID = NULL;

  if(gNeutronClass=="noSelection"){
    neutronID = new RooRealVar("fNclass","znClass",1,4);
  }
  else if(gNeutronClass=="0n0n"){
    neutronID = new RooRealVar("fNclass","znClass",1,1);
  }
  else if(gNeutronClass=="Xn0n"){
    neutronID = new RooRealVar("fNclass","znClass",2,3);
  }
  else if(gNeutronClass=="XnXn"){
    neutronID = new RooRealVar("fNclass","znClass",4,4);
  }


  // import the tree
  RooDataSet inData("inData", "inData", RooArgSet(*neutronID, pt, mass, phiAverage, rapidity), Import(*dataTree));
  // check the structure of the data
  inData.Print();

  // build pdfs
  // --> j/psi
  RooRealVar m0("m0", "m0", 3.0967,2.9,3.2);
  RooRealVar sL("sigmaL", "sigmaL",sigmaL,0.01,0.2);
  RooRealVar sR("sigmaR", "sigmaR",sigmaR,0.01,0.2);
  if (isSigmaFixed) {
    sL.setConstant(true);
    sR.setConstant(true);    
  }
  RooRealVar nL("nL", "nL", nParL,1,20);    
  RooRealVar nR("nR", "nR", nParR,1,20);    
  if (isNFixed) {
    nL.setConstant(true);
    nR.setConstant(true);    
  }
  RooRealVar aL("alphaL", "alphaL", alphaL,0.1,5);
  RooRealVar aR("alphaR", "alphaR", alphaR,0.1,5);  
  if (isAlphaFixed) {
    aL.setConstant(true);
    aR.setConstant(true);    
  }
  RooCrystalBall jpsi("jpsi", "jpsi", mass, m0, sL, sR, aL, nL, aR, nR);

  // --> psi2s
  RooFormulaVar m02("m02","@0+3.686097-3.096900",RooArgList(m0));
  // 1.09 = sqrt(Mpsi2s/Mjpsi)
  RooFormulaVar sL2("sigmaL2","@0*(1.09)",RooArgList(sL));
  RooFormulaVar sR2("sigmaR2","@0*(1.09)",RooArgList(sR));
  RooFormulaVar nL2("nL2","@0",RooArgList(nL));
  RooFormulaVar nR2("nR2","@0",RooArgList(nR));
  RooFormulaVar aL2("aL2","@0",RooArgList(aL));
  RooFormulaVar aR2("aR2","@0",RooArgList(aR));
  RooCrystalBall psi2s("psi2s", "psi2s", mass, m02, sL2, sR2, aL2, nL2, aR2, nR2);

  // --> non-resonant background
  RooRealVar b("#lambda","exponent",-lambda,-3,-0.2);
  RooExponential nrBg("nrBg","nrBg",mass,b);

  // model and fit
  // --> combine the PDFs
  int nEvents = inData.numEntries();
  RooRealVar nJpsi("N_{J/#psi}","Number of J/psi events",0.8*nEvents,0.05*nEvents,nEvents);
  RooRealVar nPsi2s("N_{#psi'}","Number of psi(2s) events",0.05*nEvents,0,nEvents);
  RooRealVar nBg("N_{bg}","Number of BG events",0.1*nEvents,0,nEvents);

  RooAddPdf *fitData = NULL;
  if (includePsi2s) {
    fitData = new RooAddPdf("fitData","J/psi, psi(2s), and background PDF", RooArgList(jpsi,psi2s,nrBg), RooArgList(nJpsi,nPsi2s,nBg));
  } else {
    fitData = new RooAddPdf("fitData","J/psi and background PDF", RooArgList(jpsi,nrBg), RooArgList(nJpsi,nBg));
  }
  if(excludeJPsi){
    fitData = new RooAddPdf("fitData","Background PDF", RooArgList(nrBg), RooArgList(nBg));
  }

  // --> perform the fit
  RooFitResult *r = NULL;
  r = fitData->fitTo(inData,Save(),Extended(kTRUE),Save());

  // check the status of the fit
  int fitStatus = r->status();
  cout<<"Fit status == "<<fitStatus<<endl;

  // Get the output
  // Create a frame with only the histo and the fit curve
  int massBins = 60;
  RooPlot *frameFit = mass.frame(Title("Frame with data and full fit curve"), Bins(massBins));
  // plot the data
  inData.plotOn(frameFit, Binning(massBins), DataError(RooAbsData::SumW2));
  // Plot the full model
  fitData->plotOn(frameFit, LineColor(kBlack));

  // Now compute the residuals and pulls
  RooHist* residualHist = frameFit->residHist();
  residualHist->SetTitle("Residuals");
  RooHist* pullHist     = frameFit->pullHist();
  pullHist->SetTitle("Pulls");

  // use this fram to compute the chi2
  int nParams = r->floatParsFinal().getSize();
  int ndf = massBins - nParams;
  double chi2_red = frameFit->chiSquare(nParams);

  // Get the number of points in the pull histogram
  int nPoints = pullHist->GetN();
  double chi2Pull = 0.0;
  // Get pointer to the array of y-values (the pulls)
  double* pulls = pullHist->GetY();
  // Loop over all points and sum up the squares of the pulls
  for (int i = 0; i < nPoints; ++i) {
    chi2Pull += pulls[i] * pulls[i];
  }
  std::cout << "Computed chi2 from pulls: " << chi2Pull << std::endl;

  // --> Draw residuals and pulls
  TCanvas *cStat = new TCanvas("cStat", "cStat", 800, 600);
  cStat->Divide(2,0);
  cStat->cd(1);
    residualHist->Draw("AP");
    cStat->Update();
  cStat->cd(2);
    pullHist->Draw("AP");
    cStat->Update();
  

  // --> Now create a fram with all the components
  RooPlot *frame = mass.frame(Title(Form("m_{#mu#mu} mass distr. - %s - #phi %d", gNeutronClass.c_str(), (int)binID[0])),Bins(massBins));
  // draw the data
  inData.plotOn(frame,Binning(massBins),DataError(RooAbsData::SumW2));
  // draw the fit
  fitData->plotOn(frame, LineColor(kBlack));
  // draw the components
  if(!excludeJPsi){
    fitData->plotOn(frame,Name("jpsi"), Components(jpsi), LineColor(kRed+1));
    if (includePsi2s) fitData->plotOn(frame,Name("psi2s"), Components(psi2s), LineColor(kGreen+2));
  }  
  fitData->plotOn(frame,Name("nrBg"), Components(nrBg), LineColor(kBlue+1));
  fitData->paramOn(frame,Layout(0.58, 0.88, 0.8));

  // --> do the plot
  TCanvas *c = new TCanvas(Form("m_{#mu#mu} mass distr. - %s - #phi %d", gNeutronClass.c_str(), (int)binID[0]),
                           Form("m_{#mu#mu} mass distr. - %s - #phi %d", gNeutronClass.c_str(), (int)binID[0]),1920,1080);

  // Create two pads: pad1 for the fit and pad2
  TPad *pad1 = new TPad("pad1", "Fit Pad", 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "Pulls Pad", 0.0, 0.05, 1.0, 0.28);

  // Draw the pads on the canvas
  c->cd();
  pad1->Draw();
  pad2->Draw();

  // Draw fit result on pad1
  pad1->cd();
  pad1->SetBottomMargin(1);
  frame->GetXaxis()->SetTitle();
  frame->Draw();

  // Add chi2 and some info of the fit
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04); // Set text size
  latex.DrawLatex(gXpos, gYpos, Form("#chi^{2}/ndf = %.2f (%.2f/%d)", chi2_red,chi2_red*ndf,ndf));
  latex.DrawLatex(gXpos, gYpos-0.05, Form("# entries = %.d", inData.numEntries()));
  latex.DrawLatex(gXpos, gYpos-0.10, Form("#it{p}_{T} in (%.2f, %.2f) GeV/#it{c}",minPt,maxPt));
  latex.DrawLatex(gXpos, gYpos-0.15, Form("m_{#mu#mu} in (%.2f, %.2f) GeV/#it{c}^{2}",minMass,maxMass));
  latex.DrawLatex(gXpos, gYpos-0.20, Form("Y in (%.2f, %.2f)",minRapidity, maxRapidity));
  latex.DrawLatex(gXpos, gYpos-0.25, Form("neutron class = %s",gNeutronClass.c_str()));

  float degMin = binID[1]*180./TMath::Pi();
  float degMax = binID[2]*180./TMath::Pi();
  latex.DrawLatex(gXpos, gYpos-0.30, Form("#phi ~ in = (%d#circ, %d#circ)", (int)degMin, (int)degMax));

  pad1->Update();

  // Draw the pulls on pad2
  pad2->cd();
  pad2->SetBottomMargin(0.28);
  pullHist->SetTitle("");
  pullHist->GetYaxis()->SetTitleSize(0.1);
  pullHist->GetYaxis()->SetLabelSize(0.1);
  pullHist->GetYaxis()->SetTitleOffset(0.35);
  pullHist->GetYaxis()->SetLabelOffset(0.008);
  pullHist->GetXaxis()->SetTitleSize(0.12);
  pullHist->GetXaxis()->SetLabelSize(0.12);
  pullHist->GetXaxis()->SetTitleOffset(1.1);
  pullHist->GetXaxis()->SetLabelOffset(0.01);

  pullHist->Draw();
  pad2->Update();

  // --> correlation matrix
  TCanvas *cCM = new TCanvas(Form("CM pt in (%.2f,%.2f) - #phi %d",minPt,maxPt,(int)binID[0]),
			     Form("CM pt in (%.2f,%.2f) - #phi %d",minPt,maxPt,(int)binID[0]),
			     200,100,900, 600);
  TH2* h_CorrM = r->correlationHist();
  h_CorrM->SetMarkerSize(1.2);
  h_CorrM->Draw("zcol,text");

  // --> print some values
  cout << " chi2/ndf = " << chi2_red << endl;
  cout << " Entries  = " << inData.numEntries() << endl;
  cout << " nJpsi    = " << nJpsi.getVal() << " ± " << nJpsi.getError() << endl;

  gStyle->SetOptFit(1);

  // --> save a PDF file with the result
  if(nCalls==0) c->Print(Form("%s/mass-fits-%d-%s.pdf[",gIdentifier.c_str(),gPhiBins,gNeutronClass.c_str()));
  c->Print(Form("%s/mass-fits-%d-%s.pdf",gIdentifier.c_str(),gPhiBins,gNeutronClass.c_str()));
  if(nCalls==gPhiBins-1) c->Print(Form("%s/mass-fits-%d-%s.pdf]",gIdentifier.c_str(),gPhiBins,gNeutronClass.c_str()));

  saveFile = new TFile(Form("%s/%s/%s",gIdentifier.c_str(),gNeutronClass.c_str(),gSaveFileName.c_str()),"update");

  // save the results of the fit in a tree
  createSaveFitTree(gSaveTreeName);

  // fill the variables in the tree
  numJPsi = nJpsi.getVal();
  errNumJPsi = nJpsi.getError();
  numBkg = nBg.getVal();
  errNumBkg = nBg.getError();
  entries = inData.numEntries();
  errEntries = TMath::Sqrt(entries);

  meanPhiAverage = (binID[1] + binID[2])/2;
  saveFitTree->Fill();
  saveFile->Write();
  saveFile->Close();

  nCalls++;

} // end of doOneDataFit

// -----------------------------------------------------------------
// entry point: set up and call fitting function
void fitJPsiInPhiBins(string nClass = "noSelection", int nPhiBins = 3, const char *config = "jPsi", bool isMC = false, bool notShow = false)
{
  // choose to show or not the canvas
  if(notShow) gROOT->SetBatch(kTRUE);

  // set the global values to the values in input
  // identifier = name of the config, it identifies the kine for a certain process
  gIdentifier = config;
  // nutron class
  gNeutronClass = nClass;
  // number of bins in phi
  gPhiBins = nPhiBins;

  // get the tree with the data
  TFile *dataFile = NULL;
  TTree *dataTree = NULL;

  if(!isMC){
    dataFile = new TFile("../Data/merged-data.root");
    dataTree = (TTree *) dataFile->Get("DF_2336518079279520/O2dimu");
  }
  else if(isMC){
    dataFile = new TFile("../MC-test/reco.root");
    dataTree = (TTree *) dataFile->Get("DF_2336518085565599/dimu");
    gNeutronClass = "noSelection";
  }
  
  // check if the data are there
  if(dataTree==NULL){
    cout<<"Data file missing. Bye!"<<endl;
    return;
  }
  // read the config file
  const char* configFile = "config.cfg";
  TEnv *env = new TEnv(configFile);
  
  // check if the config is there
  if(!filesystem::exists(configFile)){
    cout<<"Config file missing. Bye!"<<endl;
    return;
  }

  // construct the keys using the input configuration name.
  TString keyPtLow   = TString::Format("%s.minPt", config);
  TString keyPtUp    = TString::Format("%s.maxPt", config);
  TString keyMassLow = TString::Format("%s.minMass", config);
  TString keyMassUp  = TString::Format("%s.maxMass", config);
  TString keyRapLow  = TString::Format("%s.minRapidity", config);
  TString keyRaPUp   = TString::Format("%s.maxRapidity", config);

  // retireve the values: set the values to the ones in the config
  minPt       = env->GetValue(keyPtLow.Data(), 0.0);
  maxPt       = env->GetValue(keyPtUp.Data(), 0.0);
  minMass     = env->GetValue(keyMassLow.Data(), 0.0);
  maxMass     = env->GetValue(keyMassUp.Data(), 0.0);
  minRapidity = env->GetValue(keyRapLow.Data(), 0.0);
  maxRapidity = env->GetValue(keyRaPUp.Data(), 0.0);
  

  // chose parameters for the fit
  // --> crystal ball
  sigmaL = 0.08;
  sigmaR = 0.07;
  alphaL = 1.2;
  alphaR = 2.5;
  nParL  = 10;
  nParR  = 10;    
  isSigmaFixed = false;
  isAlphaFixed = true;
  isNFixed     = true;
  includePsi2s = false;
  // --> exponential
  lambda = 0.7;
  isLambdaFixed = false;

  // minimization
  isChi2Fit = false;

	// defining the ranges for phi
  float phiMin = -TMath::Pi();
  float phiMax = TMath::Pi();
  float rangePhi = phiMax-phiMin;
  float incrementPhi = rangePhi/(float)nPhiBins;
  float actualPhiMin = phiMin;
  float actualPhiMax = phiMin + incrementPhi;

	float binID[3];

  // save the results in a tree
  string fileName = to_string(nPhiBins) + "phi_1pt";
  gSaveFileName = "massFitRes_" + fileName + ".root";
  gSystem->mkdir(gIdentifier.c_str());
  gSystem->mkdir(Form("%s/%s",gIdentifier.c_str(),gNeutronClass.c_str()));
  TFile *saveFile = new TFile(Form("%s/%s/%s",gIdentifier.c_str(),gNeutronClass.c_str(),gSaveFileName.c_str()),"recreate");
	saveFile->Close();
  

	// loop on phi bins
	for(int i=1; i<=nPhiBins; i++){

		binID[0] = i;
    binID[1] = actualPhiMin;
    binID[2] = actualPhiMax;

    gSaveTreeName = "tree_phi" + to_string(i) + "_pt1";

		cout<<"min == "<<binID[1]<<endl;
		cout<<"max == "<<binID[2]<<endl;

		// do the fit
  	doOneDataFit(dataTree, saveFile, binID);

		actualPhiMin = actualPhiMin + incrementPhi;
    actualPhiMax = actualPhiMax + incrementPhi;

	} // end of the loop on the bins

} // end of the main
