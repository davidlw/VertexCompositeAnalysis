//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 18 17:08:18 2021 by ROOT version 6.06/01
// from TTree ParticleTree/ParticleTree
// found on file: /eos/cms/store/group/phys_heavyions/yousen/RiceHIN/OpenHF2020_LamCKsP/MC/LambdaC-KsPr_LCpT-0p9_pPb-EmbEPOS_8p16_Pythia8/PA8TeV_pPb_LamCKsP0p9_pT0p9to6p1_y1p1_MC_Training_20210217/Merged_PA8TeV_pPb_LamCKsP0p9_pT0p9to6p1_y1p1_MC_Training_20210217.root
//////////////////////////////////////////////////////////

#ifndef ParticleTreeMC2_hxx
#define ParticleTreeMC2_hxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

#include "ParticleTreeMC.hxx"

class ParticleTreeMC2 : public ParticleTreeMC {
  public :
    ParticleTreeMC2(TTree *tree=0);
    virtual ~ParticleTreeMC2();

    std::vector<std::vector<float>>& gen_trackParameters() const { return *_gen_trackParameters; }
    std::vector<std::vector<float>>& trk_trackParameters() const { return *_trk_trackParameters; }
    std::vector<std::vector<float>>& trk_trackCovariance() const { return *_trk_trackCovariance; }

  private:
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    std::vector<std::vector<float> > *_gen_trackParameters;
    std::vector<std::vector<float> > *_trk_trackParameters;
    std::vector<std::vector<float> > *_trk_trackCovariance;

    // List of branches
    TBranch        *b_gen_trackParameters;   //!
    TBranch        *b_trk_trackParameters;   //!
    TBranch        *b_trk_trackCovariance;   //!

    virtual void     Init(TTree *tree);
};

#endif
