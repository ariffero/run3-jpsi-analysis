//
// makes a fit to the transverse mass distribution 
//

// -----------------------------------------------------------------
// all headers are defined here
// -----------------------------------------------------------------

// c++ headers
#include <iostream>
#include <fstream>
#include <vector>

// root headers
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1D.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>

// RooFit headers
#include "RooTFnBinding.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooAbsReal.h"
#include "RooGenericPdf.h"
#include "RooAbsArg.h"
#include "RooFormulaVar.h"

using namespace RooFit;

// -----------------------------------------------------------------
// global variables, to be set to drive the fits
// -----------------------------------------------------------------

//----switches
//chose to do chi2 fit
bool gIsChi2fit = false;
//chose to do a binned fit
bool gBinned = true;
//use yy template
bool gYyTemplate = false;

// position of the legend info
float gXpos = 0.17;
float gYpos = 0.85;

// --> number of non-resonant pairs from a mass fit in the pt range
double numMuMu = 5235; 

// --> feed-down parameters from UPC measurement
double f_DC = 0.08; // feed down from coherent psip
double f_DC_err = 0.01; 
double f_DI = f_DC; // feed down from incoherent psip
double f_DI_err = f_DC_err; 

// --> H1 parameters for dissociation from https://arxiv.org/abs/1304.5162
double bH1 = 1.79; 
double nH1 = 3.58; 

// --> coherent pt range (GeV)
double ptCohRangeMin = 0; 
double ptCohRangeMax = 0.25;

vector<double> gPtFitBinning;

// do one data fit
void doOneDataFit(TH1D *dataH, RooDataSet data,
                  RooDataSet jpsiCoh, RooDataSet jpsiIncoh,
                  RooDataSet psi2sFdCoh, RooDataSet psi2sFdIncoh,
                  RooDataSet mumu, RooRealVar pt);

// -----------------------------------------------------------------
// entry point
void ptFit(vector<double> &ptFitBinning, bool yyTemplate, bool binned, bool isChi2Fit, double massRange[], double rapidityRange[], string configName)
{

  gROOT->SetBatch(kFALSE);

  // set passed vales to global variables
  gYyTemplate = yyTemplate;
  gPtFitBinning = ptFitBinning;
  gIsChi2fit = isChi2Fit;
  gBinned = binned;
  gBinned = true;
  
  //cout<<"pt = "<<gPtFitBinning[0]<<" "<<gPtFitBinning[gPtFitBinning.size()-1]<<endl;

  //-----Read the file for the data----------------------
  TFile *fData = NULL;
  TH1D *hData = NULL;
  TTree *dataTree = NULL;
  if(!gYyTemplate){
    string ptDataFile = "ptFitData-" + configName + ".root";
    fData = new TFile(ptDataFile.c_str());
    hData = (TH1D*) fData->Get("hPtFit");
  }
  else if(gYyTemplate){
    fData = new TFile("Inputs/merged-data.root");
    dataTree = (TTree*)fData->Get("DF_2336518079279520/O2dimu");
  }

  //-----Read the files with the MC----------------------
  TFile *fCohJpsi = new TFile("Inputs/jpsi_coh_reco_tree.root");
  TTree *jpsiCohTree = (TTree*)fCohJpsi->Get("DF_2336518085565631/dimu");

  TFile *fIncohJpsi = new TFile("Inputs/jpsi_incoh_reco_tree.root");
  TTree *jpsiIncohTree = (TTree*)fIncohJpsi->Get("DF_2336518085565599/dimu");

  TFile *fCohpsi2sFd = new TFile("Inputs/psi2s_coh_fd_reco_tree.root");
  TTree *psi2sFdCohTree = (TTree*)fCohpsi2sFd->Get("DF_2336518085565631/dimu");

  TFile *fIncohpsi2sFd = new TFile("Inputs/psi2s_incoh_fd_reco_tree.root");
  TTree *psi2sFdIncohTree = (TTree*)fIncohpsi2sFd->Get("DF_2336518085565599/dimu");

  TTree *mumuTree = NULL;
  TFile *fMumu = NULL;
  if(gYyTemplate){
    fMumu = new TFile("Inputs/mumu_mid_reco_tree.root");
    mumuTree = (TTree*)fMumu->Get("DF_2336518081075359/dimu");  
  }
  
  // Construct the RooDataSet
  // Start defining the RooRealVar
  RooRealVar pt("fPt","p_{T} (GeV/c)",gPtFitBinning[0],gPtFitBinning[gPtFitBinning.size()-1]);
  RooRealVar mass("fM","m_{#mu#mu} (GeV/c^{2})", massRange[0], massRange[1]);
  RooRealVar rapidity("fRap","rapidity", rapidityRange[0], rapidityRange[1]);

  RooDataSet jpsiCoh("jpsiCoh", "jpsiCoh", RooArgSet(pt,mass,rapidity), Import(*jpsiCohTree));
  RooDataSet jpsiIncoh("jpsiIncoh", "jpsiIncoh", RooArgSet(pt,mass,rapidity), Import(*jpsiIncohTree));
  RooDataSet psi2sFdCoh("psi2sFdCoh", "psi2sFdCoh", RooArgSet(pt,mass,rapidity), Import(*psi2sFdCohTree));
  RooDataSet psi2sFdIncoh("psi2sFdIncoh", "psi2sFdIncoh", RooArgSet(pt,mass,rapidity), Import(*psi2sFdIncohTree));
  RooDataSet mumu("mumu", "mumu", RooArgSet(pt,mass,rapidity), Import(*mumuTree));

  RooDataSet data("data", "data", RooArgSet(pt,mass,rapidity), Import(*dataTree));

  //------Do the fit-----------------------------------------
  doOneDataFit(hData, data, jpsiCoh, jpsiIncoh, psi2sFdCoh, psi2sFdIncoh, mumu, pt);
}

// do one data fit: implementation
void doOneDataFit(TH1D *hData, RooDataSet data,
                  RooDataSet jpsiCoh, RooDataSet jpsiIncoh,
                  RooDataSet psi2sFdCoh, RooDataSet psi2sFdIncoh,
                  RooDataSet mumu, RooRealVar pt)
{
  //------Create RooDataHist for data--------------------------------------
  RooBinning ptBinning(gPtFitBinning.size()-1,&gPtFitBinning[0]);
  if(gBinned) pt.setBinning(ptBinning);

  RooDataHist *Hist_Data = NULL;
  if(gYyTemplate) Hist_Data = new RooDataHist("Hist_Data","Hist_Data",RooArgSet(pt),data);
  else if(!gYyTemplate) Hist_Data = new RooDataHist("Hist_Data","Hist_Data",RooArgSet(pt),hData);

  //------Create MC pt templates--------------------------------------------
  //--> coherent jpsi
  RooDataHist Hist_MCdata_CohJpsi("Hist_MCdata_CohJpsi","Hist_MCdata_CohJpsi",RooArgSet(pt),jpsiCoh);
  RooHistPdf  PDF_CohJpsi("PDF_CohJpsi","PDF_CohJpsi",pt,Hist_MCdata_CohJpsi,0);
  //--> incoherent jpsi
  RooDataHist Hist_MCdata_IncohJpsi("Hist_MCdata_IncohJpsi","Hist_MCdata_IncohJpsi",RooArgSet(pt),jpsiIncoh);
  RooHistPdf  PDF_IncohJpsi("PDF_IncohJpsi","PDF_IncohJpsi",pt,Hist_MCdata_IncohJpsi,0);
  //--> coherent psip
  RooDataHist Hist_MCdata_CohPsip("Hist_MCdata_CohPsip","Hist_MCdata_CohPsip",RooArgSet(pt),psi2sFdCoh);
  RooHistPdf  PDF_CohPsip("PDF_CohPsip","PDF_CohPsip",pt,Hist_MCdata_CohPsip,0);
  //--> incoherent psip
  RooDataHist Hist_MCdata_IncohPsip("Hist_MCdata_IncohPsip","Hist_MCdata_IncohPsip",RooArgSet(pt),psi2sFdIncoh);
  RooHistPdf  PDF_IncohPsip("PDF_IncohPsip","PDF_IncohPsip",pt,Hist_MCdata_IncohPsip,0);
  //--> yy->mumu
  RooDataHist *Hist_MCdata_Mumu = NULL;
  RooHistPdf  *PDF_Mumu = NULL;
  if(gYyTemplate){
    Hist_MCdata_Mumu = new RooDataHist("Hist_MCdata_Mumu","Hist_MCdata_Mumu",RooArgSet(pt),mumu);
    PDF_Mumu = new RooHistPdf("PDF_Mumu","PDF_Mumu",pt,*Hist_MCdata_Mumu,0);
  }
  
  // ------Create Incoherent Disocitation PDF--------------------------------------------
  //--> Values and formula taken from https://arxiv.org/abs/1304.5162
  RooRealVar b("b","b",bH1, 0.1, 5);
  RooRealVar n("n","n",nH1);
  n.setConstant(kTRUE);
  //RooGenericPdf PDF_dissJpsi("PDF_dissJpsi","pt*pow((1+pow(pt,2)*b/n),-n)",RooArgSet(pt, b, n));
  RooGenericPdf PDF_dissJpsi("PDF_dissJpsi", "@0 * pow((1 + pow(@0,2) * @1 / @2), -@2)", RooArgSet(pt, b, n));

  // ------Create the model--------------------------------------------
  //--> j/psi
  int nEvts = 0;
  if(gYyTemplate) nEvts = data.sumEntries();
  if(!gYyTemplate) nEvts = hData->Integral();

  RooRealVar N_CohJpsi("N_CohJpsi","number of coherent jpsi events",0.5*nEvts,0.3*nEvts,nEvts);
  RooRealVar N_IncohJpsi("N_IncohJps","number of incoherent jpsi events",0.1*nEvts,0,nEvts);
  //--> Psi' feed down fixed by JPsi values * feed down coeficient
  /*
  RooRealVar N_CohPsip("N_CohPsip","number of coherent psip events",0.05*nEvts,0,nEvts); 
  RooRealVar N_IncohPsip("N_IncohPsip","number of incoherent psip events",0.01*nEvts,0,nEvts);
  */
  // code to use the constraint from data
  RooRealVar var_f_D_coh("var_f_D_coh","var_f_D_coh",f_DC); // from UPC analysis
  var_f_D_coh.setError(f_DC_err);
  var_f_D_coh.setConstant(kTRUE);
  RooFormulaVar N_CohPsip("N_CohPsip","@0*@1",RooArgList(N_CohJpsi, var_f_D_coh));
  RooRealVar var_f_D_incoh("var_f_D_incoh","var_f_D_incoh",f_DI); // from UPC analysis
  var_f_D_incoh.setError(f_DI_err);
  var_f_D_incoh.setConstant(kTRUE);
  RooFormulaVar N_IncohPsip("N_IncohPsip","@0*@1",RooArgList(N_IncohJpsi, var_f_D_incoh));
  //--> Disociative
  RooRealVar N_DissJpsi("N_DissJpsi","number of dissociated jpsi events",0.05*nEvts,0,nEvts);
  //--> Muons from yy
  RooRealVar N_Muons("N_Muons","number of yy->mumu events",0.5*nEvts,0.2*nEvts,nEvts);
  //--> Now the model
  RooAddPdf *Pt_fit_func = NULL;
  if(!gYyTemplate){
    Pt_fit_func = new RooAddPdf("Pt_fit_func","Sum of templates",
      RooArgList(PDF_CohJpsi,PDF_IncohJpsi, PDF_CohPsip, PDF_IncohPsip, PDF_dissJpsi),
      RooArgList(N_CohJpsi,N_IncohJpsi, N_CohPsip, N_IncohPsip, N_DissJpsi));
  }
  
  if(gYyTemplate){
    Pt_fit_func = new RooAddPdf("Pt_fit_func","Sum of templates",
      RooArgList(PDF_CohJpsi,PDF_IncohJpsi, PDF_CohPsip, PDF_IncohPsip, PDF_dissJpsi,*PDF_Mumu),
      RooArgList(N_CohJpsi,N_IncohJpsi, N_CohPsip, N_IncohPsip, N_DissJpsi,N_Muons));
  }
  

  // ------Fit model to data--------------------------------------------
  RooFitResult* fit_pt = NULL;
  if(!gIsChi2fit) fit_pt = Pt_fit_func->fitTo(*Hist_Data,Extended(kTRUE),SumW2Error(kTRUE),Save(),Range(""));
  if(gIsChi2fit)  fit_pt = Pt_fit_func->chi2FitTo(*Hist_Data,Extended(kTRUE),SumW2Error(kTRUE),Save(),Range(""));

  // ------plot fit------------------------------------------
  TCanvas *cPt = new TCanvas("cPt","cPt",800,800);
  RooPlot* frame_pt = pt.frame(Title("Pt fit")) ;
  //--> data and fit

  Hist_Data->plotOn(frame_pt,Name("Hist_Data"), Binning(ptBinning),MarkerStyle(20),MarkerSize(0.9),RooFit::DataError(RooAbsData::SumW2));

  //--> full fit
  Pt_fit_func->plotOn(frame_pt,Name("Pt_fit_func"), LineColor(kBlack), LineWidth(2));
  //--> jpsi MC
  Pt_fit_func->plotOn(frame_pt,Name("PDF_CohJpsi"), Components(PDF_CohJpsi),
		     LineStyle(1), LineColor(kBlue), LineWidth(2), Name("coherent j/#psi"));

  Pt_fit_func->plotOn(frame_pt,Name("PDF_IncohJpsi"), Components(PDF_IncohJpsi),
		     LineStyle(1), LineColor(kRed), LineWidth(2), Name("incoherent j/#psi"));
  // --> psip MC
  Pt_fit_func->plotOn(frame_pt,Name("PDF_CohPsip"), Components(PDF_CohPsip),
		     LineStyle(1), LineColor(kOrange-7), LineWidth(2), Name("coherent #psi(2s)"));

  Pt_fit_func->plotOn(frame_pt,Name("PDF_IncohPsip"), Components(PDF_IncohPsip),
		     LineStyle(1), LineColor(kOrange), LineWidth(2), Name("incoherent #psi(2s)"));
  // --> diss
  Pt_fit_func->plotOn(frame_pt,Name("PDF_dissJpsi"), Components(PDF_dissJpsi), 
         LineStyle(1), LineColor(kMagenta), LineWidth(2), Name("dissociative j/#psi"));
  
  if(gYyTemplate){
    Pt_fit_func->plotOn(frame_pt,Name("PDF_Mumu"), Components(*PDF_Mumu), 
         LineStyle(1), LineColor(kGreen), LineWidth(2), Name("Continuum"));
  }
  //Pt_fit_func.paramOn(frame_pt);

  // --> draw
  frame_pt->Draw();
  gPad->SetLogy();

  double chi2 = frame_pt->chiSquare(); 
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04); // Set text size
  latex.DrawLatex(gXpos, gYpos, Form("#chi^{2}/ndf = %.2f", chi2));

  TLegend *legend = new TLegend(0.6, 0.6, 0.85, 0.85); // Adjust position (x1, y1, x2, y2)
  legend->SetBorderSize(0); // No border
  legend->SetFillStyle(0);  // Transparent background

  // Add entries for each component
  legend->AddEntry(frame_pt->findObject("Hist_Data"), "Data", "P");
  legend->AddEntry(frame_pt->findObject("Pt_fit_func"), "Total Fit", "L");
  legend->AddEntry(frame_pt->findObject("coherent j/#psi"), "Coherent J/#psi", "L");
  legend->AddEntry(frame_pt->findObject("incoherent j/#psi"), "Incoherent J/#psi", "L");
  legend->AddEntry(frame_pt->findObject("coherent #psi(2s)"), "Coherent J/#psi from #psi(2S) decay", "L");
  legend->AddEntry(frame_pt->findObject("incoherent #psi(2s)"), "Incoherent J/#psi from #psi(2S) decay", "L");
  legend->AddEntry(frame_pt->findObject("dissociative j/#psi"), "Dissociative J/#psi", "L");
  if(gYyTemplate)legend->AddEntry(frame_pt->findObject("Continuum"), "Continuum #gamma#gamma->#mu#mu", "L");

  // Draw the legend
  legend->Draw();

  // ------plot correlation matrix--------------------------
  TCanvas *cCm = new TCanvas("cCm","cCm",800,0,800,800);
  cCm->SetTopMargin(0.03);
  cCm->SetBottomMargin(0.11);
  cCm->SetRightMargin(0.17);
  cCm->SetLeftMargin(0.15);
  TH2* hCm = fit_pt->correlationHist();
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  hCm->Draw("colz,text");

  // ------compute fI and fD fractions--------------------------
  //--> range where the fraction will be computed
  pt.setRange("ptCohRange",ptCohRangeMin,ptCohRangeMax);

  //--> integrals in this range
  RooAbsReal *iCohJpsi = PDF_CohJpsi.createIntegral(pt,NormSet(pt),Range("ptCohRange")); 
  RooAbsReal *iCohPsip = PDF_CohPsip.createIntegral(pt,NormSet(pt),Range("ptCohRange"));
  RooAbsReal *iIncJpsi = PDF_IncohJpsi.createIntegral(pt,NormSet(pt),Range("ptCohRange"));
  RooAbsReal *iIncPsip = PDF_IncohPsip.createIntegral(pt,NormSet(pt),Range("ptCohRange"));
  RooAbsReal *iDissJpsi = PDF_dissJpsi.createIntegral(pt,NormSet(pt),Range("ptCohRange"));
  RooAbsReal *iContinuum = NULL;
  if(gYyTemplate) iContinuum = PDF_Mumu->createIntegral(pt,NormSet(pt),Range("ptCohRange"));
  //--> events in this range
  double nCohJpsi = iCohJpsi->getVal()*N_CohJpsi.getVal();
  double nCohJpsiErr = iCohJpsi->getVal()*N_CohJpsi.getError();
  double nCohPsip = iCohPsip->getVal()*N_CohPsip.getVal();
  // double nCohPsipErr = iCohPsip->getVal()*N_CohPsip.getError();
  double nIncJpsi = iIncJpsi->getVal()*N_IncohJpsi.getVal();
  double nIncJpsiErr = iIncJpsi->getVal()*N_IncohJpsi.getError();
  double nIncPsip = iIncPsip->getVal()*N_IncohPsip.getVal();
  // double nIncPsipErr = iIncPsip->getVal()*N_IncohPsip.getError();
  double nDissJpsi = iDissJpsi->getVal()*N_DissJpsi.getVal();
  double nDissJpsiErr = iDissJpsi->getVal()*N_DissJpsi.getError();
  //muons from yy
  double nContinuum = 0;
  if(gYyTemplate) nContinuum = iContinuum->getVal()*N_Muons.getVal();

  //--> fI = nInc/nCoh
  double fI = (nIncJpsi+nDissJpsi)/nCohJpsi;
  double fIe = 0; // to be computed taking into account correlations
  //--> fD = nPsip/nCoh
  double fD = (nCohPsip+nIncPsip)/nCohJpsi;
  double fDe = 0; // to be computed taking into account correlations
  // print the values out
  cout << " nCohJpsi events in coherent range = " << nCohJpsi << " ± " << nCohJpsiErr << endl;
  cout << " nCohPsip events in coherent range = " << nCohPsip << endl;
  cout << " nIncJpsi events in coherent range = " << nIncJpsi << " ± " << nIncJpsiErr << endl;
  cout << " nIncPsip events in coherent range = " << nIncPsip << endl;
  cout << " continuum events in coherent range = " << nContinuum << endl;
  cout << " nDissJpsi events in coherent range = " << nDissJpsi << " ± " << nDissJpsiErr << endl;
  cout << " J/Psi yield  in coherent range = " << (nCohJpsi+nCohPsip+nIncJpsi+nIncPsip+nDissJpsi) << endl;
  cout << " Fractions in coherent range: fI = " << fI << " fD = " << fD << endl;
}
