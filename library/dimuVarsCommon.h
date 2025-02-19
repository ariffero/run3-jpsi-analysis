// variables saved in data and reco MC trees

#ifndef DIMU_VARS_COMMON_H
#define DIMU_VARS_COMMON_H

  // DATA ------------------------------------------------------------

  // TTree globl pionter
  TTree *globalTree = NULL;

  // Declaration of leaf types
  Float_t         fM;
  Float_t         fPt;
  Float_t         fRap;
  Float_t         fPhi;
  Float_t         fPhiAv;
  Float_t         fPhiCh;
  Float_t         fPtp;
  Float_t         fEtap;
  Float_t         fPhip;
  Float_t         fPtn;
  Float_t         fEtan;
  Float_t         fPhin;
  Int_t           fNclass;

  void createDataTree(TFile *globalFile, string treePath){
    globalTree = (TTree*)globalFile->Get(treePath.c_str());

    globalTree->SetBranchAddress("fM", &fM);
    globalTree->SetBranchAddress("fPt", &fPt);
    globalTree->SetBranchAddress("fRap", &fRap);
    globalTree->SetBranchAddress("fPhi", &fPhi);
    globalTree->SetBranchAddress("fPhiAv", &fPhiAv);
    globalTree->SetBranchAddress("fPhiCh", &fPhiCh);
    globalTree->SetBranchAddress("fPtp", &fPtp);
    globalTree->SetBranchAddress("fEtap", &fEtap);
    globalTree->SetBranchAddress("fPhip", &fPhip);
    globalTree->SetBranchAddress("fPtn", &fPtn);
    globalTree->SetBranchAddress("fEtan", &fEtan);
    globalTree->SetBranchAddress("fPhin", &fPhin);
    globalTree->SetBranchAddress("fNclass", &fNclass);
  }
  // ----------------------------------------------------------------

  // RECO MC --------------------------------------------------------
  // Leaf types added in reco MC
  Float_t         fGenPt;
  Float_t         fGenRap;
  Float_t         fGenPhi;
  Float_t         fGenPhiAv;
  Float_t         fGenPhiCh;
  Float_t         fGenPtp;
  Float_t         fGenEtap;
  Float_t         fGenPhip;
  Float_t         fGenPtn;
  Float_t         fGenEtan;
  Float_t         fGenPhin;

  void createRecoTree(TFile *globalFile, string treePath){
    globalTree = (TTree*)globalFile->Get(treePath.c_str());

    globalTree->SetBranchAddress("fM",     &fM);
    globalTree->SetBranchAddress("fPt",    &fPt);
    globalTree->SetBranchAddress("fRap",   &fRap);
    globalTree->SetBranchAddress("fPhi",   &fPhi);
    globalTree->SetBranchAddress("fPhiAv", &fPhiAv);
    globalTree->SetBranchAddress("fPhiCh", &fPhiCh);
    globalTree->SetBranchAddress("fPtp",   &fPtp);
    globalTree->SetBranchAddress("fEtap",  &fEtap);
    globalTree->SetBranchAddress("fPhip",  &fPhip);
    globalTree->SetBranchAddress("fPtn",   &fPtn);
    globalTree->SetBranchAddress("fEtan",  &fEtan);
    globalTree->SetBranchAddress("fPhin",  &fPhin);

    globalTree->SetBranchAddress("fGenPt",    &fGenPt);
    globalTree->SetBranchAddress("fGenRap",   &fGenRap);
    globalTree->SetBranchAddress("fGenPhi",   &fGenPhi);
    globalTree->SetBranchAddress("fGenPtp",   &fGenPtp);
    globalTree->SetBranchAddress("fGenEtap",  &fGenEtap);
    globalTree->SetBranchAddress("fGenPhip",  &fGenPhip);
    globalTree->SetBranchAddress("fGenPtn",   &fGenPtn);
    globalTree->SetBranchAddress("fGenEtan",  &fGenEtan);
    globalTree->SetBranchAddress("fGenPhin",  &fGenPhin);
  }
  // ----------------------------------------------------------------

#endif