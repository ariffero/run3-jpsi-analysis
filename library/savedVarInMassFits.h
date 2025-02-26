#ifndef SAVEDVARINMASSFITS_H
#define SAVEDVARINMASSFITS_H

  // TTree globl pionter
  TTree *saveFitTree = NULL;

  // fit results
  // j/psi
  float numJPsi;
  float errNumJPsi;
  //background
  float numBkg;
  float errNumBkg;
  // entries
  float entries;
  float errEntries;
  // AxECohJPsi = AxE for coh j/psi process
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

  void createSaveFitTree(string name){
    saveFitTree = new TTree (name.c_str(),name.c_str());

    // defining the branches of the tree
    // j/psi
    saveFitTree->Branch("numJPsi",&numJPsi,"numJPsi/F");
    saveFitTree->Branch("errNumJPsi",&errNumJPsi,"errNumJPsi/F");
    // background
    saveFitTree->Branch("numBkg",&numBkg,"numBkg/F");
    saveFitTree->Branch("errNumBkg",&errNumBkg,"errNumBkg/F");
    // AxECohJPsi
    saveFitTree->Branch("AxECohJPsi",&AxECohJPsi,"AxECohJPsi/F");
    saveFitTree->Branch("errAxECohJPsi",&errAxECohJPsi,"errAxECohJPsi/F");
    // j/psi corrected by AxECohJPsi = j/psi /AxECohJPsi
    saveFitTree->Branch("numJPsiCorr",&numJPsiCorr,"numJPsiCorr/F");
    saveFitTree->Branch("errNumJPsiCorr",&errNumJPsiCorr,"errNumJPsiCorr/F");
    // AxEMumuMid
    saveFitTree->Branch("AxEMumuMid",&AxEMumuMid,"AxEMumuMid/F");
    saveFitTree->Branch("errAxEMumuMid",&errAxEMumuMid,"errAxEMumuMid/F");
    // background corrected by AxE = bkg / AxEMumuMid
    saveFitTree->Branch("numBkgCorr",&numBkgCorr,"numBkgCorr/F");
    saveFitTree->Branch("errNumBkgCorr",&errNumBkgCorr,"errNumBkgCorr/F");

    // entries
    saveFitTree->Branch("entries",&entries,"entries/F");
    saveFitTree->Branch("errEntries",&errEntries,"errEntries/F");

    saveFitTree->Branch("meanPt",&meanPt,"meanPt/F");
    saveFitTree->Branch("meanPhiAverage",&meanPhiAverage,"meanPhiAverage/F");
  }
#endif