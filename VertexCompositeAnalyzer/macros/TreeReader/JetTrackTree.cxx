#ifndef JetTrackTree_cxx
#define JetTrackTree_cxx
#include "JetTrackTree.hxx"
#include <iostream>

using std::cout;
using std::endl;

JetTrackTree::JetTrackTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("v0ana_4497.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("v0ana_4497.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("v0ana_4497.root:/jetanalyzer");
      dir->GetObject("JetTrackTree",tree);

   }
   Init(tree);
}

JetTrackTree::~JetTrackTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetTrackTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetTrackTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

Long64_t JetTrackTree::GetEntries()
{
  if (!fChain) return -1;
  return fChain->GetEntries();
}

Long64_t JetTrackTree::GetEntriesFast()
{
  if (!fChain) return -1;
  return fChain->GetEntriesFast();
}

void JetTrackTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   _xVtx = 0;
   _nTracksVtx = 0;
   _ptSumVtx = 0;
   _trkPt = 0;
   _trkEta = 0;
   _trkPhi = 0;
   _trkCharge = 0;
   _trkPDFId = 0;
   _trkNHits = 0;
   _highPurity = 0;
   _trkAssociatedVtxIndx = 0;
   _jetNumDaughters = 0;
   _jetEta = 0;
   _jetPt = 0;
   _jetPhi = 0;
   _jetTheta = 0;
   _jetMass = 0;
   _muonMultiplicity = 0;
   _chargedMultiplicity = 0;
   _dau_pt_sum = 0;
   _dau_chg = 0;
   _dau_pid = 0;
   _dau_vref = 0;
   _dau_pt = 0;
   _dau_eta = 0;
   _dau_phi = 0;
   _dau_theta = 0;
   _dau_vz = 0;
   _dau_vy = 0;
   _dau_vx = 0;
   _dau_vrefz = 0;
   _dau_vrefy = 0;
   _dau_vrefx = 0;
   _dau_vp_difZ = 0;
   _dau_vp_difY = 0;
   _dau_vp_difX = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nRun", &_nRun, &b_nRun);
   fChain->SetBranchAddress("nEv", &_nEv, &b_nRun);
   fChain->SetBranchAddress("nLumi", &_nLumi, &b_nLumi);
   fChain->SetBranchAddress("xVtx", &_xVtx, &b_xVtx);
   fChain->SetBranchAddress("nTracksVtx", &_nTracksVtx, &b_nTracksVtx);
   fChain->SetBranchAddress("ptSumVtx", &_ptSumVtx, &b_ptSumVtx);
   fChain->SetBranchAddress("trkPt", &_trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkEta", &_trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", &_trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("trkCharge", &_trkCharge, &b_trkCharge);
   fChain->SetBranchAddress("trkPDFId", &_trkPDFId, &b_trkPDFId);
   fChain->SetBranchAddress("trkNHits", &_trkNHits, &b_trkNHits);
   fChain->SetBranchAddress("highPurity", &_highPurity, &b_highPurity);
   fChain->SetBranchAddress("trkAssociatedVtxIndx", &_trkAssociatedVtxIndx, &b_trkAssociatedVtxIndx);
   fChain->SetBranchAddress("jetNumDaughters", &_jetNumDaughters, &b_jetNumDaughters);
   fChain->SetBranchAddress("jetEta", &_jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPt", &_jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetPhi", &_jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetTheta", &_jetTheta, &b_jetTheta);
   fChain->SetBranchAddress("jetMass", &_jetMass, &b_jetMass);
   fChain->SetBranchAddress("muonMultiplicity", &_muonMultiplicity, &b_muonMultiplicity);
   fChain->SetBranchAddress("chargedMultiplicity", &_chargedMultiplicity, &b_chargedMultiplicity);
   fChain->SetBranchAddress("jetN", &_jetN, &b_jetN);
   fChain->SetBranchAddress("dau_pt_sum", &_dau_pt_sum, &b_dau_pt_sum);
   fChain->SetBranchAddress("dau_chg", &_dau_chg, &b_dau_chg);
   fChain->SetBranchAddress("dau_pid", &_dau_pid, &b_dau_pid);
   fChain->SetBranchAddress("dau_vref", &_dau_vref, &b_dau_vref);
   fChain->SetBranchAddress("dau_pt", &_dau_pt, &b_dau_pt);
   fChain->SetBranchAddress("dau_eta", &_dau_eta, &b_dau_eta);
   fChain->SetBranchAddress("dau_phi", &_dau_phi, &b_dau_phi);
   fChain->SetBranchAddress("dau_theta", &_dau_theta, &b_dau_theta);
   fChain->SetBranchAddress("dau_vz", &_dau_vz, &b_dau_vz);
   fChain->SetBranchAddress("dau_vy", &_dau_vy, &b_dau_vy);
   fChain->SetBranchAddress("dau_vx", &_dau_vx, &b_dau_vx);
   fChain->SetBranchAddress("dau_vrefz", &_dau_vrefz, &b_dau_vrefz);
   fChain->SetBranchAddress("dau_vrefy", &_dau_vrefy, &b_dau_vrefy);
   fChain->SetBranchAddress("dau_vrefx", &_dau_vrefx, &b_dau_vrefx);
   fChain->SetBranchAddress("dau_vp_difZ", &_dau_vp_difZ, &b_dau_vp_difZ);
   fChain->SetBranchAddress("dau_vp_difY", &_dau_vp_difY, &b_dau_vp_difY);
   fChain->SetBranchAddress("dau_vp_difX", &_dau_vp_difX, &b_dau_vp_difX);
   fChain->SetBranchAddress("didHLTFire", &_didHLTFire, &b_didHLTFire);
   Notify();
}

Bool_t JetTrackTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetTrackTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetTrackTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JetTrackTree_cxx
