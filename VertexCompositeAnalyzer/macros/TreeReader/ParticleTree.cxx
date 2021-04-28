#ifndef ParticleTree_cxx
#define ParticleTree_cxx
#include "ParticleTree.hxx"
#include <iostream>

using std::cout;
using std::endl;

ParticleTree::ParticleTree(TTree *tree) : fChain(tree) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    fChain == nullptr;
    return;
  }
  Init(tree);
}

ParticleTree::~ParticleTree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  fChain = nullptr;
}

Int_t ParticleTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t ParticleTree::GetEntries()
{
  if (!fChain) return -1;
  return fChain->GetEntries();
}

Long64_t ParticleTree::GetEntriesFast()
{
  if (!fChain) return -1;
  return fChain->GetEntriesFast();
}

Long64_t ParticleTree::LoadTree(Long64_t entry)
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

void ParticleTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  _evtSel = 0;
  _cand_charge = 0;
  _cand_pdgId = 0;
  _cand_status = 0;
  _cand_trkIdx = 0;
  _cand_angle2D = 0;
  _cand_angle3D = 0;
  _cand_dca = 0;
  _cand_decayLength2D = 0;
  _cand_decayLength3D = 0;
  _cand_decayLengthError2D = 0;
  _cand_decayLengthError3D = 0;
  _cand_eta = 0;
  _cand_mass = 0;
  _cand_p = 0;
  _cand_pT = 0;
  _cand_phi = 0;
  _cand_pseudoDecayLengthError2D = 0;
  _cand_pseudoDecayLengthError3D = 0;
  _cand_vtxChi2 = 0;
  _cand_vtxProb = 0;
  _cand_y = 0;
  _cand_dauIdx = 0;
  _cand_momIdx = 0;
  _cand_etaDau = 0;
  _cand_massDau = 0;
  _cand_pTDau = 0;
  _cand_phiDau = 0;
  _trk_isHP = 0;
  _trk_nHit = 0;
  _trk_dEdx_dedxHarmonic2 = 0;
  _trk_dEdx_dedxPixelHarmonic2 = 0;
  _trk_nChi2 = 0;
  _trk_pTErr = 0;
  _trk_xyDCASignificance = 0;
  _trk_zDCASignificance = 0;
  _trk_candIdx = 0;

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("nPV", &_nPV, &b_nPV);
  fChain->SetBranchAddress("BXNb", &_BXNb, &b_BXNb);
//  fChain->SetBranchAddress("Ntrkoffline", &_Ntrkoffline, &b_Ntrkoffline);
  fChain->SetBranchAddress("EventNb", &_EventNb, &b_EventNb);
  fChain->SetBranchAddress("LSNb", &_LSNb, &b_LSNb);
  fChain->SetBranchAddress("RunNb", &_RunNb, &b_RunNb);
//  fChain->SetBranchAddress("HFsumETMinus", &_HFsumETMinus, &b_HFsumETMinus);
//  fChain->SetBranchAddress("HFsumETPlus", &_HFsumETPlus, &b_HFsumETPlus);
//  fChain->SetBranchAddress("Npixel", &_Npixel, &b_Npixel);
//  fChain->SetBranchAddress("ZDCMinus", &_ZDCMinus, &b_ZDCMinus);
//  fChain->SetBranchAddress("ZDCPlus", &_ZDCPlus, &b_ZDCPlus);
  fChain->SetBranchAddress("bestvtxX", &_bestvtxX, &b_bestvtxX);
  fChain->SetBranchAddress("bestvtxY", &_bestvtxY, &b_bestvtxY);
  fChain->SetBranchAddress("bestvtxZ", &_bestvtxZ, &b_bestvtxZ);
  fChain->SetBranchAddress("evtSel", &_evtSel, &b_evtSel);
  fChain->SetBranchAddress("cand_charge", &_cand_charge, &b_cand_charge);
  fChain->SetBranchAddress("cand_pdgId", &_cand_pdgId, &b_cand_pdgId);
  fChain->SetBranchAddress("cand_status", &_cand_status, &b_cand_status);
  fChain->SetBranchAddress("cand_trkIdx", &_cand_trkIdx, &b_cand_trkIdx);
  fChain->SetBranchAddress("cand_angle2D", &_cand_angle2D, &b_cand_angle2D);
  fChain->SetBranchAddress("cand_angle3D", &_cand_angle3D, &b_cand_angle3D);
  fChain->SetBranchAddress("cand_dca", &_cand_dca, &b_cand_dca);
  fChain->SetBranchAddress("cand_decayLength2D", &_cand_decayLength2D, &b_cand_decayLength2D);
  fChain->SetBranchAddress("cand_decayLength3D", &_cand_decayLength3D, &b_cand_decayLength3D);
  fChain->SetBranchAddress("cand_decayLengthError2D", &_cand_decayLengthError2D, &b_cand_decayLengthError2D);
  fChain->SetBranchAddress("cand_decayLengthError3D", &_cand_decayLengthError3D, &b_cand_decayLengthError3D);
  fChain->SetBranchAddress("cand_eta", &_cand_eta, &b_cand_eta);
  fChain->SetBranchAddress("cand_mass", &_cand_mass, &b_cand_mass);
  fChain->SetBranchAddress("cand_p", &_cand_p, &b_cand_p);
  fChain->SetBranchAddress("cand_pT", &_cand_pT, &b_cand_pT);
  fChain->SetBranchAddress("cand_phi", &_cand_phi, &b_cand_phi);
  fChain->SetBranchAddress("cand_pseudoDecayLengthError2D", &_cand_pseudoDecayLengthError2D, &b_cand_pseudoDecayLengthError2D);
  fChain->SetBranchAddress("cand_pseudoDecayLengthError3D", &_cand_pseudoDecayLengthError3D, &b_cand_pseudoDecayLengthError3D);
  fChain->SetBranchAddress("cand_vtxChi2", &_cand_vtxChi2, &b_cand_vtxChi2);
  fChain->SetBranchAddress("cand_vtxProb", &_cand_vtxProb, &b_cand_vtxProb);
  fChain->SetBranchAddress("cand_y", &_cand_y, &b_cand_y);
  fChain->SetBranchAddress("cand_dauIdx", &_cand_dauIdx, &b_cand_dauIdx);
  fChain->SetBranchAddress("cand_momIdx", &_cand_momIdx, &b_cand_momIdx);
  fChain->SetBranchAddress("cand_etaDau", &_cand_etaDau, &b_cand_etaDau);
  fChain->SetBranchAddress("cand_massDau", &_cand_massDau, &b_cand_massDau);
  fChain->SetBranchAddress("cand_pTDau", &_cand_pTDau, &b_cand_pTDau);
  fChain->SetBranchAddress("cand_phiDau", &_cand_phiDau, &b_cand_phiDau);
  fChain->SetBranchAddress("trk_isHP", &_trk_isHP, &b_trk_isHP);
  fChain->SetBranchAddress("trk_nHit", &_trk_nHit, &b_trk_nHit);
//  fChain->SetBranchAddress("trk_dEdx_dedxHarmonic2", &_trk_dEdx_dedxHarmonic2, &b_trk_dEdx_dedxHarmonic2);
//  fChain->SetBranchAddress("trk_dEdx_dedxPixelHarmonic2", &_trk_dEdx_dedxPixelHarmonic2, &b_trk_dEdx_dedxPixelHarmonic2);
  fChain->SetBranchAddress("trk_nChi2", &_trk_nChi2, &b_trk_nChi2);
  fChain->SetBranchAddress("trk_pTErr", &_trk_pTErr, &b_trk_pTErr);
  fChain->SetBranchAddress("trk_xyDCASignificance", &_trk_xyDCASignificance, &b_trk_xyDCASignificance);
  fChain->SetBranchAddress("trk_zDCASignificance", &_trk_zDCASignificance, &b_trk_zDCASignificance);
  fChain->SetBranchAddress("trk_candIdx", &_trk_candIdx, &b_trk_candIdx);
  Notify();
}

Bool_t ParticleTree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void ParticleTree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t ParticleTree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef ParticleTree_cxx
