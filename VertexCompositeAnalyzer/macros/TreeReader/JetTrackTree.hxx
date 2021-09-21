//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 30 15:00:08 2021 by ROOT version 6.18/04
// from TTree JetTrackTree/v1
// found on file: v0ana_4497.root
//////////////////////////////////////////////////////////

#ifndef JetTrackTree_h
#define JetTrackTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class JetTrackTree {
public :
   JetTrackTree(TTree *tree=0);
   virtual ~JetTrackTree();

   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Long64_t GetEntries();
   virtual Long64_t GetEntriesFast();

   Int_t           nRun() const { return _nRun; }
   Int_t           nEv() const { return _nEv; }
   Int_t           nLumi() const { return _nLumi; }
   std::vector<float>&   xVtx() const { return *_xVtx; }
   std::vector<int>&     nTracksVtx() const { return *_nTracksVtx; }
   std::vector<float>&   ptSumVtx() const { return *_ptSumVtx; }
   std::vector<float>&   trkPt() const { return *_trkPt; }
   std::vector<float>&   trkEta() const { return *_trkEta; }
   std::vector<float>&   trkPhi() const { return *_trkPhi; }
   std::vector<char>&    trkCharge() const { return *_trkCharge; }
   std::vector<int>&     trkPDFId() const { return *_trkPDFId; }
   std::vector<char>&    trkNHits() const { return *_trkNHits; }
   std::vector<bool>&    highPurity() const { return *_highPurity; }
   std::vector<int>&     trkAssociatedVtxIndx() const { return *_trkAssociatedVtxIndx; }
   std::vector<int>&     jetNumDaughters() const { return *_jetNumDaughters; }
   std::vector<float>&   jetEta() const { return *_jetEta; }
   std::vector<float>&   jetPt() const { return *_jetPt; }
   std::vector<float>&   jetPhi() const { return *_jetPhi; }
   std::vector<float>&   jetTheta() const { return *_jetTheta; }
   std::vector<float>&   jetMass() const { return *_jetMass; }
   std::vector<int>&     muonMultiplicity() const { return *_muonMultiplicity; }
   std::vector<int>&     chargedMultiplicity() const { return *_chargedMultiplicity; }
   Int_t           jetN() const { return _jetN; };
   std::vector<float>&   dau_pt_sum() const { return *_dau_pt_sum; }
   std::vector<std::vector<int> >& dau_chg() const { return *_dau_chg; }
   std::vector<std::vector<int> >& dau_pid() const { return *_dau_pid; }
   std::vector<std::vector<unsigned int> >& dau_vref() const { return *_dau_vref; }
   std::vector<std::vector<float> >& dau_pt() const { return *_dau_pt; }
   std::vector<std::vector<float> >& dau_eta() const { return *_dau_eta; }
   std::vector<std::vector<float> >& dau_phi() const { return *_dau_phi; }
   std::vector<std::vector<float> >& dau_theta() const { return *_dau_theta; }
   std::vector<std::vector<float> >& dau_vz() const { return *_dau_vz; }
   std::vector<std::vector<float> >& dau_vy() const { return *_dau_vy; }
   std::vector<std::vector<float> >& dau_vx() const { return *_dau_vx; }
   std::vector<std::vector<float> >& dau_vrefz() const { return *_dau_vrefz; }
   std::vector<std::vector<float> >& dau_vrefy() const { return *_dau_vrefy; }
   std::vector<std::vector<float> >& dau_vrefx() const { return *_dau_vrefx; }
   std::vector<std::vector<float> >& dau_vp_difZ() const { return *_dau_vp_difZ; }
   std::vector<std::vector<float> >& dau_vp_difY() const { return *_dau_vp_difY; }
   std::vector<std::vector<float> >& dau_vp_difX() const { return *_dau_vp_difX; }
   Bool_t          didHLTFire() const { return _didHLTFire; }

  protected:
    virtual void     Init(TTree *tree);

  private:
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           _nRun;
   Int_t           _nEv;
   Int_t           _nLumi;
   std::vector<float>   *_xVtx;
   std::vector<int>     *_nTracksVtx;
   std::vector<float>   *_ptSumVtx;
   std::vector<float>   *_trkPt;
   std::vector<float>   *_trkEta;
   std::vector<float>   *_trkPhi;
   std::vector<char>    *_trkCharge;
   std::vector<int>     *_trkPDFId;
   std::vector<char>    *_trkNHits;
   std::vector<bool>    *_highPurity;
   std::vector<int>     *_trkAssociatedVtxIndx;
   std::vector<int>     *_jetNumDaughters;
   std::vector<float>   *_jetEta;
   std::vector<float>   *_jetPt;
   std::vector<float>   *_jetPhi;
   std::vector<float>   *_jetTheta;
   std::vector<float>   *_jetMass;
   std::vector<int>     *_muonMultiplicity;
   std::vector<int>     *_chargedMultiplicity;
   Int_t           _jetN;
   std::vector<float>   *_dau_pt_sum;
   std::vector<std::vector<int> > *_dau_chg;
   std::vector<std::vector<int> > *_dau_pid;
   std::vector<std::vector<unsigned int> > *_dau_vref;
   std::vector<std::vector<float> > *_dau_pt;
   std::vector<std::vector<float> > *_dau_eta;
   std::vector<std::vector<float> > *_dau_phi;
   std::vector<std::vector<float> > *_dau_theta;
   std::vector<std::vector<float> > *_dau_vz;
   std::vector<std::vector<float> > *_dau_vy;
   std::vector<std::vector<float> > *_dau_vx;
   std::vector<std::vector<float> > *_dau_vrefz;
   std::vector<std::vector<float> > *_dau_vrefy;
   std::vector<std::vector<float> > *_dau_vrefx;
   std::vector<std::vector<float> > *_dau_vp_difZ;
   std::vector<std::vector<float> > *_dau_vp_difY;
   std::vector<std::vector<float> > *_dau_vp_difX;
   Bool_t          _didHLTFire;

   // List of branches
   TBranch        *b_nRun;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_nTracksVtx;   //!
   TBranch        *b_ptSumVtx;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkPDFId;   //!
   TBranch        *b_trkNHits;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_trkAssociatedVtxIndx;   //!
   TBranch        *b_jetNumDaughters;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetTheta;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_muonMultiplicity;   //!
   TBranch        *b_chargedMultiplicity;   //!
   TBranch        *b_jetN;   //!
   TBranch        *b_dau_pt_sum;   //!
   TBranch        *b_dau_chg;   //!
   TBranch        *b_dau_pid;   //!
   TBranch        *b_dau_vref;   //!
   TBranch        *b_dau_pt;   //!
   TBranch        *b_dau_eta;   //!
   TBranch        *b_dau_phi;   //!
   TBranch        *b_dau_theta;   //!
   TBranch        *b_dau_vz;   //!
   TBranch        *b_dau_vy;   //!
   TBranch        *b_dau_vx;   //!
   TBranch        *b_dau_vrefz;   //!
   TBranch        *b_dau_vrefy;   //!
   TBranch        *b_dau_vrefx;   //!
   TBranch        *b_dau_vp_difZ;   //!
   TBranch        *b_dau_vp_difY;   //!
   TBranch        *b_dau_vp_difX;   //!
   TBranch        *b_didHLTFire;   //!
};

#endif
