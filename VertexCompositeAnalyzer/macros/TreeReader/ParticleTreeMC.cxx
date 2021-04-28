#ifndef ParticleTreeMC_cxx
#define ParticleTreeMC_cxx
#include <iostream>
#include "ParticleTreeMC.hxx"

using std::cout;
using std::endl;

ParticleTreeMC::ParticleTreeMC(TTree *tree) : fChain(tree)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    fChain == nullptr;
    return;
  }
  ParticleTree::Init(tree);
  Init(tree);
}

ParticleTreeMC::~ParticleTreeMC()
{
}

void ParticleTreeMC::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  _cand_matchGEN = 0;
  _cand_momMatchGEN = 0;
  _cand_genPdgId = 0;
  _cand_isSwap = 0;
  _cand_genIdx = 0;
  _cand_momMatchIdx = 0;
  _gen_charge = 0;
  _gen_pdgId = 0;
  _gen_status = 0;
  _gen_statusBit = 0;
  _gen_angle2D = 0;
  _gen_angle3D = 0;
  _gen_decayLength2D = 0;
  _gen_decayLength3D = 0;
  _gen_eta = 0;
  _gen_mass = 0;
  _gen_p = 0;
  _gen_pT = 0;
  _gen_phi = 0;
  _gen_y = 0;
  _gen_dauIdx = 0;
  _gen_momIdx = 0;
  _gen_candIdx = 0;

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("Ntrkgen", &_Ntrkgen, &b_Ntrkgen);
  fChain->SetBranchAddress("genWeight", &_genWeight, &b_genWeight);
  fChain->SetBranchAddress("pTHat", &_pTHat, &b_pTHat);
  fChain->SetBranchAddress("cand_matchGEN", &_cand_matchGEN, &b_cand_matchGEN);
  fChain->SetBranchAddress("cand_momMatchGEN", &_cand_momMatchGEN, &b_cand_momMatchGEN);
  fChain->SetBranchAddress("cand_genPdgId", &_cand_genPdgId, &b_cand_genPdgId);
  fChain->SetBranchAddress("cand_isSwap", &_cand_isSwap, &b_cand_isSwap);
  fChain->SetBranchAddress("cand_genIdx", &_cand_genIdx, &b_cand_genIdx);
  fChain->SetBranchAddress("cand_momMatchIdx", &_cand_momMatchIdx, &b_cand_momMatchIdx);
  fChain->SetBranchAddress("gen_charge", &_gen_charge, &b_gen_charge);
  fChain->SetBranchAddress("gen_pdgId", &_gen_pdgId, &b_gen_pdgId);
  fChain->SetBranchAddress("gen_status", &_gen_status, &b_gen_status);
  fChain->SetBranchAddress("gen_statusBit", &_gen_statusBit, &b_gen_statusBit);
  fChain->SetBranchAddress("gen_angle2D", &_gen_angle2D, &b_gen_angle2D);
  fChain->SetBranchAddress("gen_angle3D", &_gen_angle3D, &b_gen_angle3D);
  fChain->SetBranchAddress("gen_decayLength2D", &_gen_decayLength2D, &b_gen_decayLength2D);
  fChain->SetBranchAddress("gen_decayLength3D", &_gen_decayLength3D, &b_gen_decayLength3D);
  fChain->SetBranchAddress("gen_eta", &_gen_eta, &b_gen_eta);
  fChain->SetBranchAddress("gen_mass", &_gen_mass, &b_gen_mass);
  fChain->SetBranchAddress("gen_p", &_gen_p, &b_gen_p);
  fChain->SetBranchAddress("gen_pT", &_gen_pT, &b_gen_pT);
  fChain->SetBranchAddress("gen_phi", &_gen_phi, &b_gen_phi);
  fChain->SetBranchAddress("gen_y", &_gen_y, &b_gen_y);
  fChain->SetBranchAddress("gen_dauIdx", &_gen_dauIdx, &b_gen_dauIdx);
  fChain->SetBranchAddress("gen_momIdx", &_gen_momIdx, &b_gen_momIdx);
  fChain->SetBranchAddress("gen_candIdx", &_gen_candIdx, &b_gen_candIdx);
}
#endif
