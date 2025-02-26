#ifndef MASSTREEVARS_H
#define MASSTREEVARS_H

  // TTree globl pionter
  TTree *globalTree = NULL;

  // fit results
  // j/psi
  float numJPsi;
  float errNumJPsi;
  // background
  float numBkg;
  float errNumBkg;
  // entries
  float entries;
  float errEntries;
  // AxECohJPsi
  float AxECohJPsi;
  float errAxECohJPsi;
  // corrected j/psi
  float numJPsiCorr;
  float errNumJPsiCorr;
  // AxEMumuMid
  float AxEMumuMid;
  float errAxEMumuMid;
  // corrected background
  float numBkgCorr;
  float errNumBkgCorr;
  // mean pt and phi of the bin
  float meanPt;
  float meanPhiAverage;

  void createDataTree(TFile *globalFile, string treePath){

    // get the tree
    globalTree = (TTree*)globalFile->Get(treePath.c_str());
    
    // defining the branches of the tree
    // j/psi
    globalTree->SetBranchAddress("numJPsi",&numJPsi);
    globalTree->SetBranchAddress("errNumJPsi",&errNumJPsi);
    // background
    globalTree->SetBranchAddress("numBkg",&numBkg);
    globalTree->SetBranchAddress("errNumBkg",&errNumBkg);
    // AxECohJPsi
    globalTree->SetBranchAddress("AxECohJPsi",&AxECohJPsi);
    globalTree->SetBranchAddress("errAxECohJPsi",&errAxECohJPsi);
    // j/psi corrected by AxE = j/psi /AxECohJPsi
    globalTree->SetBranchAddress("numJPsiCorr",&numJPsiCorr);
    globalTree->SetBranchAddress("errNumJPsiCorr",&errNumJPsiCorr);
    // AxEMumuMid
    globalTree->SetBranchAddress("AxEMumuMid",&AxEMumuMid);
    globalTree->SetBranchAddress("errAxEMumuMid",&errAxEMumuMid);
    // background corrected by AxE = bkg / AxEMumuMid
    globalTree->SetBranchAddress("numBkgCorr",&numBkgCorr);
    globalTree->SetBranchAddress("errNumBkgCorr",&errNumBkgCorr);

    // entries
    globalTree->SetBranchAddress("entries",&entries);
    globalTree->SetBranchAddress("errEntries",&errEntries);

    globalTree->SetBranchAddress("meanPt",&meanPt);
    globalTree->SetBranchAddress("meanPhiAverage",&meanPhiAverage);
  }
  
#endif