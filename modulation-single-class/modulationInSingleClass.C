//
// makes a fit to the number of (corrected) jpsi as a function of phi
// in a single neutron class (hence it doeas not take into account migrations)
//

// -----------------------------------------------------------------
// all headers are defined here
// -----------------------------------------------------------------

// c++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>

using namespace std;

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TF1.h"
#include "TTreeFormula.h"

// my headers: used to read the trees with the mass fit results
#include "../library/massTreeVars.h"

// -----------------------------------------------------------------
// global variables
// -----------------------------------------------------------------

// vectors with the possible variables that can be plotted
vector<string> possibleVar{"numJPsi", "numJPsiCorr", "entries", "numBkg", "AxE"};
vector<string> possibleErrVar{"errNumJPsi", "errNumJPsiCorr", "errEntries", "errNumBkg", "errAxE"};
// the variable will be chosen using varID
// 0 = # j/psi (not corrected for AxE)
// 1 = # j/psi (AxE-corrected)
// 2 = entries in the tree
// 3 = # bkg (not corrected for AxE)
// 4 = AxE


// do not show the pop up of the canvas
bool notShow = false;

//global variable for neutron classes
string gNeutronClass = "";
//name of the (1st) folder that contains the results
string gIdentifier = "";
// pdf with the results
string pdfFilePath = "";
// for the name of the input files
string massFitFile = "";

// phi bins
int gPhiBins = 0;

// --> switches 
bool fitWithConstant = false;
bool doNotFit = false;


// function to normalize the histo so the mean y =1
void normY1(TH1D *h, int nBins){
  int sum = 0;
  for(int i=0;i<nBins;i++) {
  sum += h->GetBinContent(i+1); 
}
h->Scale(((double) (nBins))/sum);  
}

// function that fill the histogram and the error using the 3 vector passed 
// it fills the histo starting from the bin pStart and ends at pEnd
void fillHistoAndSetError(TH1 *h, vector<float> vX,vector<float> vY,vector<float> vErr, int pStart, int pEnd){
  int j = 1;
  for(int i=pStart; i<=pEnd; i++){
      h->Fill(vX[i],vY[i]);
      h->SetBinError(j,vErr[i]);
      j++;
  }
}

// -----------------------------------------------------------------------------
// function that return the value of the variable requested in the requested event of the tree
float getEventValFromTree(string massFitFile, string ID, string varName, int event){
  // open the file that stores the tree 
  TFile *file = new TFile(massFitFile.c_str());
  // call the function that allow to access the TTree accordingly to the tipe of the data in the file
  createDataTree(file,ID.c_str());

  float temp = -500;

  // use a TTreeFormula to get the right variable
  TTreeFormula getVar("getVar", varName.c_str(), globalTree);
  // get right event 
  globalTree->GetEvent(event);
  // get the chosen variable
  temp = getVar.EvalInstance(0);

  file->Close();

  if(temp != -500) return temp;

  cout<<"error: the requested variable hasn't been found"<<endl;
  return -500;
}

// ----------------------------------------------
// Fit model for the asymmetry
double fitAsymmetry(double *x, double *par){
  
  double phi = x[0];
  double norm = par[0];
  double a1 = par[1];
  double a2 = par[2];
  double a3 = par[3];
  double a4 = par[4];

  double asymmetry = norm*(1 + a1*TMath::Cos(phi) + a2*TMath::Cos(2*phi) + a3*TMath::Cos(3*phi) + a4*TMath::Cos(4*phi));

  return asymmetry;
}

// -------------------------------------------------------------------------------
// function that perform one fit fit and saves the results of the cos(2phi) amplitude and its err
void onePhiFit(int event, int varID){

  // number of times that the function is called
  static int calls = 0;

  // declaration of the vector used for the phi fit
  vector<float> vec_deltaPhi;
  vector<float> vec_var;
  vector<float> vec_errVar;

  // ID of the tree -> used to move along the phi bins
  string ID;

  // fill the vector
  for(int j=1; j<=gPhiBins; j++){
    ID = "tree_phi"+to_string(j)+"_pt"+to_string(1);
    //filling the vectors
    vec_deltaPhi.push_back(getEventValFromTree(massFitFile, ID, "meanPhiAverage", event));
    vec_var.push_back(getEventValFromTree(massFitFile, ID, possibleVar[varID], event));
    vec_errVar.push_back(getEventValFromTree(massFitFile, ID, possibleErrVar[varID], event));
  }

  // check that all the vectors have the same size
  cout<<"sanity check: "<<vec_deltaPhi.size()<<" "<<vec_var.size()<<" "<<vec_errVar.size()<<endl;

  // create the single histo with the var vs phi
  float xMin = -TMath::Pi();
  float xMax = TMath::Pi();
 
  TH1D *hVar = new TH1D("hVar",Form("%s modulation - %s",possibleVar[varID].c_str(), gNeutronClass.c_str()), gPhiBins, xMin, xMax);
  hVar->GetXaxis()->SetTitle("#phi");
  hVar->GetYaxis()->SetTitle(possibleVar[varID].c_str());
  //hVar->GetYaxis()->SetTitle("J/#psi yield (normalized)");

  // change the labels on the x-axis
  hVar->GetXaxis()->SetNdivisions(-504);
  hVar->GetXaxis()->ChangeLabelByValue(-TMath::Pi(),-1,-1,-1,-1,-1,"-#pi");
  hVar->GetXaxis()->ChangeLabelByValue(-TMath::Pi()/2,-1,-1,-1,-1,-1,"-#pi/2");
  hVar->GetXaxis()->ChangeLabelByValue(0,-1,-1,-1,-1,-1,"0");
  hVar->GetXaxis()->ChangeLabelByValue(TMath::Pi()/2,-1,-1,-1,-1,-1,"#pi/2");
  hVar->GetXaxis()->ChangeLabelByValue(TMath::Pi(),-1,-1,-1,-1,-1,"#pi");

  // fill the histo
  fillHistoAndSetError(hVar, vec_deltaPhi, vec_var, vec_errVar, 0, gPhiBins-1);

  // prepare a canvas
  TCanvas *cPhi = new TCanvas("cPhi","cPhi",1920,1080);

  // if the fit is not request save just the histo with the distribution
  // before normalizing it
  if(doNotFit){
    hVar->DrawClone("ep");
    cPhi->Print(pdfFilePath.c_str());
    return;
  }

  // normalize the histo
  if(varID!=4) normY1(hVar, gPhiBins); // do not normalize for the AxE

  // check the normalization
  hVar->DrawClone("ep");

  //saving the fits in a pdf file
  gErrorIgnoreLevel = kWarning;

  // --> make the fit

  // print fit info on the canvas
  gStyle->SetOptFit(1);
  // do not print stat info
  gStyle->SetOptStat(0);

  // prepare the fitting function
  int fitParameters = 5;
  if(fitWithConstant) fitParameters = 1;

  TF1 *asymmetryModel = new TF1("asymmetryModel", fitAsymmetry, TMath::Pi(), TMath::Pi(), fitParameters);
  asymmetryModel->SetLineColor(kRed);
  asymmetryModel->SetLineWidth(3);
  asymmetryModel->SetLineStyle(1);
  // setting the parameters
  asymmetryModel->SetParNames("norm", "a_{1}", "a_{2}", "a_{3}", "a_{4}");
  // norm
  asymmetryModel->FixParameter(0,1);
  // a1
  asymmetryModel->SetParameter(1,0);
  // a2
  asymmetryModel->SetParameter(2,0);
  // a3
  asymmetryModel->FixParameter(3,0);
  // a4
  asymmetryModel->FixParameter(4,0);

  if(fitWithConstant){
    asymmetryModel->FixParameter(0,1);
    asymmetryModel->FixParameter(1,0);
    asymmetryModel->FixParameter(2,0);
    asymmetryModel->FixParameter(3,0);
    asymmetryModel->FixParameter(4,0);
  }

  TFitResultPtr phiPtr = 0;

  // do the fit
  phiPtr = hVar->Fit(asymmetryModel,"RM+SEI");
  // chi2
  double chi2 = asymmetryModel->GetChisquare();
  double ndf = asymmetryModel->GetNDF();
  
  // print the correlation matrix
  phiPtr->GetCorrelationMatrix().Print();

  // print the chi2
  cout<<"chi2/ndf = "<<chi2/ndf<<" ("<<chi2<<"/"<<ndf<<")"<<endl;

  // adjust the y axis
  hVar->GetYaxis()->SetRangeUser(hVar->GetMinimum()/1.15,hVar->GetMaximum()*1.15);

  // save the result on a pdf
  cPhi->Print(pdfFilePath.c_str());

}

// --------------------------------------------------------------------------------
// Entry point: set up and call the fitting function
void modulationInSingleClass(string neutronClass = "0n0n", int varID = 1, string identifier = "jPsi", int phiBins = 12){
  
  // set the input values to the global variables
  gIdentifier = identifier;
  gNeutronClass = neutronClass;
  gPhiBins = phiBins;

  // read the trees that contains the results of the mass fit
  // strings used to read the input files
  string massFitFolder = "../mass-fits";
  string totBinID = to_string(phiBins)+"phi_1pt";
  massFitFile = Form("%s/%s/%s/massFitRes_%s.root",massFitFolder.c_str(),gIdentifier.c_str(),neutronClass.c_str(),totBinID.c_str());

  cout<<massFitFile<<endl;

  //check if the files exist
  if(!filesystem::exists(massFitFile)){
    cout<<"The file: "<<massFitFile<<" does not exist"<<endl;
    return;
  }

  // check that the vectors with the possible vars have the same size
  if(possibleErrVar.size() != possibleVar.size()){
    cout<<"The vectors of possible vars have different size! Bye."<<endl;
    return;
  }

  pdfFilePath = "mod-results-" + to_string(phiBins) + "-" + gIdentifier + "-" + possibleVar[varID] + "-" + gNeutronClass + ".pdf";

  // number of events, the structure allows for multiple events for syst study
  int event = 0; // fix it to one for the moment

  onePhiFit(event, varID);
}