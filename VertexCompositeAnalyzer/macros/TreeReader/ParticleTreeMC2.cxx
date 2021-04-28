#ifndef ParticleTreeMC2_cxx
#define ParticleTreeMC2_cxx
#include <iostream>
#include "ParticleTreeMC2.hxx"

using std::cout;
using std::endl;

ParticleTreeMC2::ParticleTreeMC2(TTree *tree) : fChain(tree)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    fChain == nullptr;
    return;
  }
  ParticleTree::Init(tree);
  ParticleTreeMC::Init(tree);
  Init(tree);
}

ParticleTreeMC2::~ParticleTreeMC2()
{
}

void ParticleTreeMC2::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  _gen_trackParameters = 0;
  _trk_trackParameters = 0;
  _trk_trackCovariance = 0;

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("gen_trackParameters", &_gen_trackParameters, &b_gen_trackParameters);
  fChain->SetBranchAddress("trk_trackParameters", &_trk_trackParameters, &b_trk_trackParameters);
  fChain->SetBranchAddress("trk_trackCovariance", &_trk_trackCovariance, &b_trk_trackCovariance);
}
#endif
