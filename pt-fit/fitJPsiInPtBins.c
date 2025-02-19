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

// root headers
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2.h>
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

// RooFit headers
#include <RooGlobalFunc.h>
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

// My headers: needed to save the results in the tree
#include "../library/savedVarInMassFits.h"

using namespace RooFit;

// -----------------------------------------------------------------
// global variables, to be set to drive the fits
// -----------------------------------------------------------------
// --> kinematics
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
  static int nCalls = 0;

  // define variables for the tree
  RooRealVar pt("fPt","pt",minPt,maxPt);
  RooRealVar mass("fM","m_{#mu#mu} (GeV/c^{2})",minMass,maxMass);
  RooRealVar rapidity("fRap","rapidity",minRapidity,maxRapidity);
	RooRealVar phiAverage("fPhiAv","phiAverage",binID[1],binID[2]);
  int znMin = -1;
  int znMax = -1;
  
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
  //RooDataSet inData("inData", "inData", RooArgSet(*neutronID, mass,), Import(*dataTree));
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
  if (isChi2Fit) r = fitData->chi2FitTo(inData,Save(),Extended(kTRUE),Save());
  else           r = fitData->fitTo(inData,Save(),Extended(kTRUE),Save());

  int fitStatus = r->status();
  cout<<"status == "<<fitStatus<<endl;

  // Get the output
  // --> set the plot of the invariant mass distribution
  RooPlot *frame = mass.frame(Title(Form("m_{#mu#mu} mass distr. - %s - #phi %d", gNeutronClass.c_str(), (int)binID[0])));
  int massBins = 50;
  inData.plotOn(frame,Binning(massBins));
  fitData->plotOn(frame, LineColor(kBlack));
  if(!excludeJPsi){
    fitData->plotOn(frame,Name("jpsi"), Components(jpsi), LineColor(kRed+1));
    if (includePsi2s) fitData->plotOn(frame,Name("psi2s"), Components(psi2s), LineColor(kGreen+2));
  }  
  fitData->plotOn(frame,Name("nrBg"), Components(nrBg), LineColor(kBlue+1));
  fitData->paramOn(frame,Layout(0.58, 0.88, 0.8));

  // --> calculate chi^2 and NDOF
  double chi2 = frame->chiSquare();  // Reduced chi^2 (chi^2/ndf)

  // --> do the plot
  TCanvas *c = new TCanvas(Form("m_{#mu#mu} mass distr. - %s - #phi %d", gNeutronClass.c_str(), (int)binID[0]),
                           Form("m_{#mu#mu} mass distr. - %s - #phi %d", gNeutronClass.c_str(), (int)binID[0]),900,600);
  frame->Draw();

  // Add chi2 and NDOF to the canvas
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04); // Set text size
  latex.DrawLatex(gXpos, gYpos, Form("#chi^{2}/ndf = %.2f", chi2));
  latex.DrawLatex(gXpos, gYpos-0.05, Form("# entries = %.d", inData.numEntries()));
  latex.DrawLatex(gXpos, gYpos-0.10, Form("#it{p}_{T} in (%.2f, %.2f) GeV/#it{c}",minPt,maxPt));
  latex.DrawLatex(gXpos, gYpos-0.15, Form("m_{#mu#mu} in (%.2f, %.2f) GeV/#it{c}^{2}",minMass,maxMass));
  latex.DrawLatex(gXpos, gYpos-0.20, Form("Y in (%.2f, %.2f)",minRapidity, maxRapidity));
  latex.DrawLatex(gXpos, gYpos-0.25, Form("neutron class = %s",gNeutronClass.c_str()));

  float degMin = binID[1]*180./TMath::Pi();
  float degMax = binID[2]*180./TMath::Pi();
  latex.DrawLatex(gXpos, gYpos-0.30, Form("#phi ~ in = (%d#circ, %d#circ)", (int)degMin, (int)degMax));


  // --> correlation matrix
  TCanvas *cCM = new TCanvas(Form("CM pt in (%.2f,%.2f) - #phi %d",minPt,maxPt,(int)binID[0]),
			     Form("CM pt in (%.2f,%.2f) - #phi %d",minPt,maxPt,(int)binID[0]),
			     200,100,900, 600);
  TH2* h_CorrM = r->correlationHist();
  h_CorrM->SetMarkerSize(1.2);
  h_CorrM->Draw("zcol,text");

  // --> print some values
  cout << "chi^2 = " << frame->chiSquare() << endl;
  cout << " Entries " << inData.numEntries() << endl;
  cout << " nJpsi = " << nJpsi.getVal() << " ± " << nJpsi.getError() << endl;

  gStyle->SetOptFit(1);

  if(nCalls==0) c->Print(Form("%s/%s/%s.pdf[",gIdentifier.c_str(),gNeutronClass.c_str(),gSaveFileName.c_str()));
  c->Print(Form("%s/%s/%s.pdf",gIdentifier.c_str(),gNeutronClass.c_str(),gSaveFileName.c_str()));
  if(nCalls==gPhiBins-1) c->Print(Form("%s/%s/%s.pdf]",gIdentifier.c_str(),gNeutronClass.c_str(),gSaveFileName.c_str()));
  saveFile = new TFile(Form("%s/%s/%s.root",gIdentifier.c_str(),gNeutronClass.c_str(),gSaveFileName.c_str()),"update");

  createSaveFitTree(gSaveTreeName);

  // fill the variables in the tree
  numJPsi = nJpsi.getVal();
  errNumJPsi = nJpsi.getError();
  numBkg = nBg.getVal();
  errNumBkg = nBg.getError();
  entries = inData.numEntries();
  errEntries = TMath::Sqrt(entries);

  meanPhiAverage = (binID[1] + binID[2])/2;
  //numJPsi = 3;
  saveFitTree->Fill();
  saveFile->Write();
  saveFile->Close();

  //save the number of j/Psi
  ofstream fout;
  fout.open(Form("jpsi-%s.txt",gIdentifier.c_str()), ios::app);
  fout<<numJPsi<<"\t"<<errNumJPsi<<"\t"<<minPt<<"\t"<<maxPt<<endl;

  nCalls++;
}

// -----------------------------------------------------------------
// entry point: set up and call fitting function
void fitJPsiInPtBins(vector<double> &ptFitBinning, int ptId, double massRange[], double rapidityRange[], string identifier)
{
  // put these as input if you need to run on different values
  int nPhiBins = 1;
  string nClass = "noSelection";
  bool isMC = false;
  bool notShow = true;

  if(notShow) gROOT->SetBatch(kTRUE);

  gIdentifier = identifier;
  gNeutronClass = nClass;
  gPhiBins = nPhiBins;

  // get the tree with the data
  TFile *dataFile = NULL;
  TTree *dataTree = NULL;

  if(!isMC){
    dataFile = new TFile("Inputs/merged-data.root");
    dataTree = (TTree *) dataFile->Get("DF_2336518079279520/O2dimu");
  }
  else if(isMC){
    dataFile = new TFile("../MC-test/reco.root");
    dataTree = (TTree *) dataFile->Get("DF_2336518085565599/dimu");
    gNeutronClass = "noSelection";
  }
  

  // chose kinematics
  minPt = ptFitBinning[ptId-1];
  maxPt = ptFitBinning[ptId];
  minMass = massRange[0];
  maxMass = massRange[1];
  minRapidity = rapidityRange[0];
  maxRapidity = rapidityRange[1];

  // chose parameters
  // --> crystal ball
  sigmaL = 0.08;
  sigmaR = 0.07;
  alphaL = 1.2;
  alphaR = 2.5;
  nParL = 10;
  nParR = 10;    
  isSigmaFixed = false;
  isAlphaFixed = true;
  isNFixed = true;
  includePsi2s = false;
  // --> exponential
  lambda = 0.7;
  isLambdaFixed = false;

  //minimization
  isChi2Fit = false;

	//defining the ranges for phi
  float phiMin = -TMath::Pi();
  float phiMax = TMath::Pi();
  float rangePhi = phiMax-phiMin;
  float incrementPhi = rangePhi/(float)nPhiBins;
  float actualPhiMin = phiMin;
  float actualPhiMax = phiMin + incrementPhi;

	float binID[3];

  // save the results in a tree
  string fileName = to_string(nPhiBins) + "phi_1pt_" + to_string(ptId) + "ptId";
  gSaveFileName = "massFitRes_" + fileName;
  gSystem->mkdir(identifier.c_str());
  gSystem->mkdir(Form("%s/%s",identifier.c_str(),gNeutronClass.c_str()));
  TFile *saveFile = new TFile(Form("%s/%s/%s.root",gIdentifier.c_str(),gNeutronClass.c_str(),gSaveFileName.c_str()),"recreate");
	saveFile->Close();
  

	//loop on phi bins
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

	}
  
}
