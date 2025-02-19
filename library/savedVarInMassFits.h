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

  // ->mean pt and phi of the bin
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
    // entries
    saveFitTree->Branch("entries",&entries,"entries/F");
    saveFitTree->Branch("errEntries",&errEntries,"errEntries/F");

    saveFitTree->Branch("meanPt",&meanPt,"meanPt/F");
    saveFitTree->Branch("meanPhiAverage",&meanPhiAverage,"meanPhiAverage/F");
  }
#endif