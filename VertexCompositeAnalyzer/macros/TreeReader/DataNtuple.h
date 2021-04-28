//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 27 10:03:42 2020 by ROOT version 6.06/01
// from TTree ParticleNTuple/ParticleNTuple
// found on file: /storage1/osg/stage_out/store/user/yousen/RiceHIN/pPb2016/LamCTree/VertexCompositeTree_PAHM6_lamCDATA_20200822/PAHighMultiplicity6/VertexCompositeTree_PAHM6_lamCDATA_20200822/200822_201043/0001/lambdacana_data_1000.root
//////////////////////////////////////////////////////////

#ifdef ParticleNTuple_h
#define ParticleNTuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ParticleNTuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Bool_t          cand_dau0_matchTRG0;
   Bool_t          cand_dau0_matchTRG1;
   Bool_t          cand_dau1_matchTRG0;
   Bool_t          cand_dau1_matchTRG1;
   Bool_t          cand_gdau0_matchTRG0;
   Bool_t          cand_gdau0_matchTRG1;
   Bool_t          cand_gdau1_matchTRG0;
   Bool_t          cand_gdau1_matchTRG1;
   Bool_t          evtSel0;
   Bool_t          evtSel1;
   Bool_t          evtSel2;
   Bool_t          evtSel3;
   Bool_t          evtSel4;
   Bool_t          evtSel5;
   Bool_t          passHLT0;
   Bool_t          passHLT1;
   Bool_t          passHLT2;
   Bool_t          passHLT3;
   Bool_t          passHLT4;
   Bool_t          passHLT5;
   Bool_t          passHLT6;
   Bool_t          passHLTPrescaler0;
   Bool_t          passHLTPrescaler1;
   Bool_t          passHLTPrescaler2;
   Bool_t          passHLTPrescaler3;
   Bool_t          passHLTPrescaler4;
   Bool_t          passHLTPrescaler5;
   Bool_t          passHLTPrescaler6;
   Bool_t          passL10;
   Bool_t          passL11;
   Bool_t          passL12;
   Bool_t          passL13;
   Bool_t          passL14;
   Bool_t          passL15;
   Bool_t          passL16;
   Bool_t          passL1Prescaler0;
   Bool_t          passL1Prescaler1;
   Bool_t          passL1Prescaler2;
   Bool_t          passL1Prescaler3;
   Bool_t          passL1Prescaler4;
   Bool_t          passL1Prescaler5;
   Bool_t          passL1Prescaler6;
   Bool_t          trk_dau0_isHP;
   Bool_t          trk_dau1_isHP;
   Bool_t          trk_gdau0_isHP;
   Bool_t          trk_gdau1_isHP;
   Bool_t          trk_isHP;
   Bool_t          validPrescale0;
   Bool_t          validPrescale1;
   Bool_t          validPrescale2;
   Bool_t          validPrescale3;
   Bool_t          validPrescale4;
   Bool_t          validPrescale5;
   Bool_t          validPrescale6;
   UChar_t         cand_dau0_status;
   UChar_t         cand_dau1_status;
   UChar_t         cand_gdau0_status;
   UChar_t         cand_gdau1_status;
   UChar_t         cand_status;
   UChar_t         hltPDs0;
   UChar_t         hltPDs1;
   UChar_t         hltPDs2;
   UChar_t         hltPDs3;
   UChar_t         hltPDs4;
   UChar_t         hltPDs5;
   UChar_t         hltPDs6;
   UChar_t         nPV;
   Char_t          cand_charge;
   Char_t          cand_dau0_charge;
   Char_t          cand_dau1_charge;
   Char_t          cand_gdau0_charge;
   Char_t          cand_gdau1_charge;
   UShort_t        BXNb;
   UShort_t        hltPrescale0;
   UShort_t        hltPrescale1;
   UShort_t        hltPrescale2;
   UShort_t        hltPrescale3;
   UShort_t        hltPrescale4;
   UShort_t        hltPrescale5;
   UShort_t        hltPrescale6;
   UShort_t        l1Prescale0;
   UShort_t        l1Prescale1;
   UShort_t        l1Prescale2;
   UShort_t        l1Prescale3;
   UShort_t        l1Prescale4;
   UShort_t        l1Prescale5;
   UShort_t        l1Prescale6;
   UShort_t        trk_dau0_nHit;
   UShort_t        trk_dau1_nHit;
   UShort_t        trk_gdau0_nHit;
   UShort_t        trk_gdau1_nHit;
   UShort_t        trk_nHit;
   UInt_t          EventNb;
   UInt_t          LSNb;
   UInt_t          RunNb;
   Int_t           cand_dau0_pdgId;
   Int_t           cand_dau1_pdgId;
   Int_t           cand_gdau0_pdgId;
   Int_t           cand_gdau1_pdgId;
   Int_t           cand_pdgId;
   Float_t         HFsumETMinus;
   Float_t         HFsumETPlus;
   Float_t         Npixel;
   Float_t         Ntrkoffline;
   Float_t         ZDCMinus;
   Float_t         ZDCPlus;
   Float_t         bestvtxX;
   Float_t         bestvtxY;
   Float_t         bestvtxZ;
   Float_t         cand_angle2D;
   Float_t         cand_angle3D;
   Float_t         cand_dau0_angle2D;
   Float_t         cand_dau0_angle3D;
   Float_t         cand_dau0_dca;
   Float_t         cand_dau0_decayLength2D;
   Float_t         cand_dau0_decayLength3D;
   Float_t         cand_dau0_decayLengthError2D;
   Float_t         cand_dau0_decayLengthError3D;
   Float_t         cand_dau0_eta;
   Float_t         cand_dau0_etaDau0;
   Float_t         cand_dau0_etaDau1;
   Float_t         cand_dau0_mass;
   Float_t         cand_dau0_massDau0;
   Float_t         cand_dau0_massDau1;
   Float_t         cand_dau0_p;
   Float_t         cand_dau0_pT;
   Float_t         cand_dau0_pTDau0;
   Float_t         cand_dau0_pTDau1;
   Float_t         cand_dau0_phi;
   Float_t         cand_dau0_phiDau0;
   Float_t         cand_dau0_phiDau1;
   Float_t         cand_dau0_pseudoDecayLengthError2D;
   Float_t         cand_dau0_pseudoDecayLengthError3D;
   Float_t         cand_dau0_vtxChi2;
   Float_t         cand_dau0_vtxProb;
   Float_t         cand_dau0_y;
   Float_t         cand_dau1_angle2D;
   Float_t         cand_dau1_angle3D;
   Float_t         cand_dau1_dca;
   Float_t         cand_dau1_decayLength2D;
   Float_t         cand_dau1_decayLength3D;
   Float_t         cand_dau1_decayLengthError2D;
   Float_t         cand_dau1_decayLengthError3D;
   Float_t         cand_dau1_eta;
   Float_t         cand_dau1_mass;
   Float_t         cand_dau1_p;
   Float_t         cand_dau1_pT;
   Float_t         cand_dau1_phi;
   Float_t         cand_dau1_pseudoDecayLengthError2D;
   Float_t         cand_dau1_pseudoDecayLengthError3D;
   Float_t         cand_dau1_vtxChi2;
   Float_t         cand_dau1_vtxProb;
   Float_t         cand_dau1_y;
   Float_t         cand_dca;
   Float_t         cand_decayLength2D;
   Float_t         cand_decayLength3D;
   Float_t         cand_decayLengthError2D;
   Float_t         cand_decayLengthError3D;
   Float_t         cand_eta;
   Float_t         cand_etaDau0;
   Float_t         cand_etaDau1;
   Float_t         cand_gdau0_angle2D;
   Float_t         cand_gdau0_angle3D;
   Float_t         cand_gdau0_dca;
   Float_t         cand_gdau0_decayLength2D;
   Float_t         cand_gdau0_decayLength3D;
   Float_t         cand_gdau0_decayLengthError2D;
   Float_t         cand_gdau0_decayLengthError3D;
   Float_t         cand_gdau0_eta;
   Float_t         cand_gdau0_mass;
   Float_t         cand_gdau0_p;
   Float_t         cand_gdau0_pT;
   Float_t         cand_gdau0_phi;
   Float_t         cand_gdau0_pseudoDecayLengthError2D;
   Float_t         cand_gdau0_pseudoDecayLengthError3D;
   Float_t         cand_gdau0_vtxChi2;
   Float_t         cand_gdau0_vtxProb;
   Float_t         cand_gdau0_y;
   Float_t         cand_gdau1_angle2D;
   Float_t         cand_gdau1_angle3D;
   Float_t         cand_gdau1_dca;
   Float_t         cand_gdau1_decayLength2D;
   Float_t         cand_gdau1_decayLength3D;
   Float_t         cand_gdau1_decayLengthError2D;
   Float_t         cand_gdau1_decayLengthError3D;
   Float_t         cand_gdau1_eta;
   Float_t         cand_gdau1_mass;
   Float_t         cand_gdau1_p;
   Float_t         cand_gdau1_pT;
   Float_t         cand_gdau1_phi;
   Float_t         cand_gdau1_pseudoDecayLengthError2D;
   Float_t         cand_gdau1_pseudoDecayLengthError3D;
   Float_t         cand_gdau1_vtxChi2;
   Float_t         cand_gdau1_vtxProb;
   Float_t         cand_gdau1_y;
   Float_t         cand_mass;
   Float_t         cand_massDau0;
   Float_t         cand_massDau1;
   Float_t         cand_p;
   Float_t         cand_pT;
   Float_t         cand_pTDau0;
   Float_t         cand_pTDau1;
   Float_t         cand_phi;
   Float_t         cand_phiDau0;
   Float_t         cand_phiDau1;
   Float_t         cand_pseudoDecayLengthError2D;
   Float_t         cand_pseudoDecayLengthError3D;
   Float_t         cand_vtxChi2;
   Float_t         cand_vtxProb;
   Float_t         cand_y;
   Float_t         pileup;
   Float_t         rawInstLumi;
   Float_t         trk_dau0_nChi2;
   Float_t         trk_dau0_pTErr;
   Float_t         trk_dau0_xyDCASignificance;
   Float_t         trk_dau0_zDCASignificance;
   Float_t         trk_dau1_nChi2;
   Float_t         trk_dau1_pTErr;
   Float_t         trk_dau1_xyDCASignificance;
   Float_t         trk_dau1_zDCASignificance;
   Float_t         trk_gdau0_nChi2;
   Float_t         trk_gdau0_pTErr;
   Float_t         trk_gdau0_xyDCASignificance;
   Float_t         trk_gdau0_zDCASignificance;
   Float_t         trk_gdau1_nChi2;
   Float_t         trk_gdau1_pTErr;
   Float_t         trk_gdau1_xyDCASignificance;
   Float_t         trk_gdau1_zDCASignificance;
   Float_t         trk_nChi2;
   Float_t         trk_pTErr;
   Float_t         trk_xyDCASignificance;
   Float_t         trk_zDCASignificance;

   // List of branches
   TBranch        *b_cand_dau0_matchTRG0;   //!
   TBranch        *b_cand_dau0_matchTRG1;   //!
   TBranch        *b_cand_dau1_matchTRG0;   //!
   TBranch        *b_cand_dau1_matchTRG1;   //!
   TBranch        *b_cand_gdau0_matchTRG0;   //!
   TBranch        *b_cand_gdau0_matchTRG1;   //!
   TBranch        *b_cand_gdau1_matchTRG0;   //!
   TBranch        *b_cand_gdau1_matchTRG1;   //!
   TBranch        *b_evtSel0;   //!
   TBranch        *b_evtSel1;   //!
   TBranch        *b_evtSel2;   //!
   TBranch        *b_evtSel3;   //!
   TBranch        *b_evtSel4;   //!
   TBranch        *b_evtSel5;   //!
   TBranch        *b_passHLT0;   //!
   TBranch        *b_passHLT1;   //!
   TBranch        *b_passHLT2;   //!
   TBranch        *b_passHLT3;   //!
   TBranch        *b_passHLT4;   //!
   TBranch        *b_passHLT5;   //!
   TBranch        *b_passHLT6;   //!
   TBranch        *b_passHLTPrescaler0;   //!
   TBranch        *b_passHLTPrescaler1;   //!
   TBranch        *b_passHLTPrescaler2;   //!
   TBranch        *b_passHLTPrescaler3;   //!
   TBranch        *b_passHLTPrescaler4;   //!
   TBranch        *b_passHLTPrescaler5;   //!
   TBranch        *b_passHLTPrescaler6;   //!
   TBranch        *b_passL10;   //!
   TBranch        *b_passL11;   //!
   TBranch        *b_passL12;   //!
   TBranch        *b_passL13;   //!
   TBranch        *b_passL14;   //!
   TBranch        *b_passL15;   //!
   TBranch        *b_passL16;   //!
   TBranch        *b_passL1Prescaler0;   //!
   TBranch        *b_passL1Prescaler1;   //!
   TBranch        *b_passL1Prescaler2;   //!
   TBranch        *b_passL1Prescaler3;   //!
   TBranch        *b_passL1Prescaler4;   //!
   TBranch        *b_passL1Prescaler5;   //!
   TBranch        *b_passL1Prescaler6;   //!
   TBranch        *b_trk_dau0_isHP;   //!
   TBranch        *b_trk_dau1_isHP;   //!
   TBranch        *b_trk_gdau0_isHP;   //!
   TBranch        *b_trk_gdau1_isHP;   //!
   TBranch        *b_trk_isHP;   //!
   TBranch        *b_validPrescale0;   //!
   TBranch        *b_validPrescale1;   //!
   TBranch        *b_validPrescale2;   //!
   TBranch        *b_validPrescale3;   //!
   TBranch        *b_validPrescale4;   //!
   TBranch        *b_validPrescale5;   //!
   TBranch        *b_validPrescale6;   //!
   TBranch        *b_cand_dau0_status;   //!
   TBranch        *b_cand_dau1_status;   //!
   TBranch        *b_cand_gdau0_status;   //!
   TBranch        *b_cand_gdau1_status;   //!
   TBranch        *b_cand_status;   //!
   TBranch        *b_hltPDs0;   //!
   TBranch        *b_hltPDs1;   //!
   TBranch        *b_hltPDs2;   //!
   TBranch        *b_hltPDs3;   //!
   TBranch        *b_hltPDs4;   //!
   TBranch        *b_hltPDs5;   //!
   TBranch        *b_hltPDs6;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_cand_charge;   //!
   TBranch        *b_cand_dau0_charge;   //!
   TBranch        *b_cand_dau1_charge;   //!
   TBranch        *b_cand_gdau0_charge;   //!
   TBranch        *b_cand_gdau1_charge;   //!
   TBranch        *b_BXNb;   //!
   TBranch        *b_hltPrescale0;   //!
   TBranch        *b_hltPrescale1;   //!
   TBranch        *b_hltPrescale2;   //!
   TBranch        *b_hltPrescale3;   //!
   TBranch        *b_hltPrescale4;   //!
   TBranch        *b_hltPrescale5;   //!
   TBranch        *b_hltPrescale6;   //!
   TBranch        *b_l1Prescale0;   //!
   TBranch        *b_l1Prescale1;   //!
   TBranch        *b_l1Prescale2;   //!
   TBranch        *b_l1Prescale3;   //!
   TBranch        *b_l1Prescale4;   //!
   TBranch        *b_l1Prescale5;   //!
   TBranch        *b_l1Prescale6;   //!
   TBranch        *b_trk_dau0_nHit;   //!
   TBranch        *b_trk_dau1_nHit;   //!
   TBranch        *b_trk_gdau0_nHit;   //!
   TBranch        *b_trk_gdau1_nHit;   //!
   TBranch        *b_trk_nHit;   //!
   TBranch        *b_EventNb;   //!
   TBranch        *b_LSNb;   //!
   TBranch        *b_RunNb;   //!
   TBranch        *b_cand_dau0_pdgId;   //!
   TBranch        *b_cand_dau1_pdgId;   //!
   TBranch        *b_cand_gdau0_pdgId;   //!
   TBranch        *b_cand_gdau1_pdgId;   //!
   TBranch        *b_cand_pdgId;   //!
   TBranch        *b_HFsumETMinus;   //!
   TBranch        *b_HFsumETPlus;   //!
   TBranch        *b_Npixel;   //!
   TBranch        *b_Ntrkoffline;   //!
   TBranch        *b_ZDCMinus;   //!
   TBranch        *b_ZDCPlus;   //!
   TBranch        *b_bestvtxX;   //!
   TBranch        *b_bestvtxY;   //!
   TBranch        *b_bestvtxZ;   //!
   TBranch        *b_cand_angle2D;   //!
   TBranch        *b_cand_angle3D;   //!
   TBranch        *b_cand_dau0_angle2D;   //!
   TBranch        *b_cand_dau0_angle3D;   //!
   TBranch        *b_cand_dau0_dca;   //!
   TBranch        *b_cand_dau0_decayLength2D;   //!
   TBranch        *b_cand_dau0_decayLength3D;   //!
   TBranch        *b_cand_dau0_decayLengthError2D;   //!
   TBranch        *b_cand_dau0_decayLengthError3D;   //!
   TBranch        *b_cand_dau0_eta;   //!
   TBranch        *b_cand_dau0_etaDau0;   //!
   TBranch        *b_cand_dau0_etaDau1;   //!
   TBranch        *b_cand_dau0_mass;   //!
   TBranch        *b_cand_dau0_massDau0;   //!
   TBranch        *b_cand_dau0_massDau1;   //!
   TBranch        *b_cand_dau0_p;   //!
   TBranch        *b_cand_dau0_pT;   //!
   TBranch        *b_cand_dau0_pTDau0;   //!
   TBranch        *b_cand_dau0_pTDau1;   //!
   TBranch        *b_cand_dau0_phi;   //!
   TBranch        *b_cand_dau0_phiDau0;   //!
   TBranch        *b_cand_dau0_phiDau1;   //!
   TBranch        *b_cand_dau0_pseudoDecayLengthError2D;   //!
   TBranch        *b_cand_dau0_pseudoDecayLengthError3D;   //!
   TBranch        *b_cand_dau0_vtxChi2;   //!
   TBranch        *b_cand_dau0_vtxProb;   //!
   TBranch        *b_cand_dau0_y;   //!
   TBranch        *b_cand_dau1_angle2D;   //!
   TBranch        *b_cand_dau1_angle3D;   //!
   TBranch        *b_cand_dau1_dca;   //!
   TBranch        *b_cand_dau1_decayLength2D;   //!
   TBranch        *b_cand_dau1_decayLength3D;   //!
   TBranch        *b_cand_dau1_decayLengthError2D;   //!
   TBranch        *b_cand_dau1_decayLengthError3D;   //!
   TBranch        *b_cand_dau1_eta;   //!
   TBranch        *b_cand_dau1_mass;   //!
   TBranch        *b_cand_dau1_p;   //!
   TBranch        *b_cand_dau1_pT;   //!
   TBranch        *b_cand_dau1_phi;   //!
   TBranch        *b_cand_dau1_pseudoDecayLengthError2D;   //!
   TBranch        *b_cand_dau1_pseudoDecayLengthError3D;   //!
   TBranch        *b_cand_dau1_vtxChi2;   //!
   TBranch        *b_cand_dau1_vtxProb;   //!
   TBranch        *b_cand_dau1_y;   //!
   TBranch        *b_cand_dca;   //!
   TBranch        *b_cand_decayLength2D;   //!
   TBranch        *b_cand_decayLength3D;   //!
   TBranch        *b_cand_decayLengthError2D;   //!
   TBranch        *b_cand_decayLengthError3D;   //!
   TBranch        *b_cand_eta;   //!
   TBranch        *b_cand_etaDau0;   //!
   TBranch        *b_cand_etaDau1;   //!
   TBranch        *b_cand_gdau0_angle2D;   //!
   TBranch        *b_cand_gdau0_angle3D;   //!
   TBranch        *b_cand_gdau0_dca;   //!
   TBranch        *b_cand_gdau0_decayLength2D;   //!
   TBranch        *b_cand_gdau0_decayLength3D;   //!
   TBranch        *b_cand_gdau0_decayLengthError2D;   //!
   TBranch        *b_cand_gdau0_decayLengthError3D;   //!
   TBranch        *b_cand_gdau0_eta;   //!
   TBranch        *b_cand_gdau0_mass;   //!
   TBranch        *b_cand_gdau0_p;   //!
   TBranch        *b_cand_gdau0_pT;   //!
   TBranch        *b_cand_gdau0_phi;   //!
   TBranch        *b_cand_gdau0_pseudoDecayLengthError2D;   //!
   TBranch        *b_cand_gdau0_pseudoDecayLengthError3D;   //!
   TBranch        *b_cand_gdau0_vtxChi2;   //!
   TBranch        *b_cand_gdau0_vtxProb;   //!
   TBranch        *b_cand_gdau0_y;   //!
   TBranch        *b_cand_gdau1_angle2D;   //!
   TBranch        *b_cand_gdau1_angle3D;   //!
   TBranch        *b_cand_gdau1_dca;   //!
   TBranch        *b_cand_gdau1_decayLength2D;   //!
   TBranch        *b_cand_gdau1_decayLength3D;   //!
   TBranch        *b_cand_gdau1_decayLengthError2D;   //!
   TBranch        *b_cand_gdau1_decayLengthError3D;   //!
   TBranch        *b_cand_gdau1_eta;   //!
   TBranch        *b_cand_gdau1_mass;   //!
   TBranch        *b_cand_gdau1_p;   //!
   TBranch        *b_cand_gdau1_pT;   //!
   TBranch        *b_cand_gdau1_phi;   //!
   TBranch        *b_cand_gdau1_pseudoDecayLengthError2D;   //!
   TBranch        *b_cand_gdau1_pseudoDecayLengthError3D;   //!
   TBranch        *b_cand_gdau1_vtxChi2;   //!
   TBranch        *b_cand_gdau1_vtxProb;   //!
   TBranch        *b_cand_gdau1_y;   //!
   TBranch        *b_cand_mass;   //!
   TBranch        *b_cand_massDau0;   //!
   TBranch        *b_cand_massDau1;   //!
   TBranch        *b_cand_p;   //!
   TBranch        *b_cand_pT;   //!
   TBranch        *b_cand_pTDau0;   //!
   TBranch        *b_cand_pTDau1;   //!
   TBranch        *b_cand_phi;   //!
   TBranch        *b_cand_phiDau0;   //!
   TBranch        *b_cand_phiDau1;   //!
   TBranch        *b_cand_pseudoDecayLengthError2D;   //!
   TBranch        *b_cand_pseudoDecayLengthError3D;   //!
   TBranch        *b_cand_vtxChi2;   //!
   TBranch        *b_cand_vtxProb;   //!
   TBranch        *b_cand_y;   //!
   TBranch        *b_pileup;   //!
   TBranch        *b_rawInstLumi;   //!
   TBranch        *b_trk_dau0_nChi2;   //!
   TBranch        *b_trk_dau0_pTErr;   //!
   TBranch        *b_trk_dau0_xyDCASignificance;   //!
   TBranch        *b_trk_dau0_zDCASignificance;   //!
   TBranch        *b_trk_dau1_nChi2;   //!
   TBranch        *b_trk_dau1_pTErr;   //!
   TBranch        *b_trk_dau1_xyDCASignificance;   //!
   TBranch        *b_trk_dau1_zDCASignificance;   //!
   TBranch        *b_trk_gdau0_nChi2;   //!
   TBranch        *b_trk_gdau0_pTErr;   //!
   TBranch        *b_trk_gdau0_xyDCASignificance;   //!
   TBranch        *b_trk_gdau0_zDCASignificance;   //!
   TBranch        *b_trk_gdau1_nChi2;   //!
   TBranch        *b_trk_gdau1_pTErr;   //!
   TBranch        *b_trk_gdau1_xyDCASignificance;   //!
   TBranch        *b_trk_gdau1_zDCASignificance;   //!
   TBranch        *b_trk_nChi2;   //!
   TBranch        *b_trk_pTErr;   //!
   TBranch        *b_trk_xyDCASignificance;   //!
   TBranch        *b_trk_zDCASignificance;   //!

   ParticleNTuple(TTree *tree=0);
   virtual ~ParticleNTuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t GetEntries();
   virtual Long64_t GetEntriesFast();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

ParticleNTuple::ParticleNTuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/storage1/osg/stage_out/store/user/yousen/RiceHIN/pPb2016/LamCTree/VertexCompositeTree_PAHM6_lamCDATA_20200822/PAHighMultiplicity6/VertexCompositeTree_PAHM6_lamCDATA_20200822/200822_201043/0001/lambdacana_data_1000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/storage1/osg/stage_out/store/user/yousen/RiceHIN/pPb2016/LamCTree/VertexCompositeTree_PAHM6_lamCDATA_20200822/PAHighMultiplicity6/VertexCompositeTree_PAHM6_lamCDATA_20200822/200822_201043/0001/lambdacana_data_1000.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/storage1/osg/stage_out/store/user/yousen/RiceHIN/pPb2016/LamCTree/VertexCompositeTree_PAHM6_lamCDATA_20200822/PAHighMultiplicity6/VertexCompositeTree_PAHM6_lamCDATA_20200822/200822_201043/0001/lambdacana_data_1000.root:/lambdacAna_mc");
      dir->GetObject("ParticleNTuple",tree);

   }
   Init(tree);
}

ParticleNTuple::~ParticleNTuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ParticleNTuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t ParticleNTuple::GetEntries()
{
  if (!fChain) return 0;
  return fChain->GetEntries();
}

Long64_t ParticleNTuple::GetEntriesFast()
{
  if (!fChain) return 0;
  return fChain->GetEntriesFast();
}

Long64_t ParticleNTuple::LoadTree(Long64_t entry)
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

void ParticleNTuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("cand_dau0_matchTRG0", &cand_dau0_matchTRG0, &b_cand_dau0_matchTRG0);
   fChain->SetBranchAddress("cand_dau0_matchTRG1", &cand_dau0_matchTRG1, &b_cand_dau0_matchTRG1);
   fChain->SetBranchAddress("cand_dau1_matchTRG0", &cand_dau1_matchTRG0, &b_cand_dau1_matchTRG0);
   fChain->SetBranchAddress("cand_dau1_matchTRG1", &cand_dau1_matchTRG1, &b_cand_dau1_matchTRG1);
   fChain->SetBranchAddress("cand_gdau0_matchTRG0", &cand_gdau0_matchTRG0, &b_cand_gdau0_matchTRG0);
   fChain->SetBranchAddress("cand_gdau0_matchTRG1", &cand_gdau0_matchTRG1, &b_cand_gdau0_matchTRG1);
   fChain->SetBranchAddress("cand_gdau1_matchTRG0", &cand_gdau1_matchTRG0, &b_cand_gdau1_matchTRG0);
   fChain->SetBranchAddress("cand_gdau1_matchTRG1", &cand_gdau1_matchTRG1, &b_cand_gdau1_matchTRG1);
   fChain->SetBranchAddress("evtSel0", &evtSel0, &b_evtSel0);
   fChain->SetBranchAddress("evtSel1", &evtSel1, &b_evtSel1);
   fChain->SetBranchAddress("evtSel2", &evtSel2, &b_evtSel2);
   fChain->SetBranchAddress("evtSel3", &evtSel3, &b_evtSel3);
   fChain->SetBranchAddress("evtSel4", &evtSel4, &b_evtSel4);
   fChain->SetBranchAddress("evtSel5", &evtSel5, &b_evtSel5);
   fChain->SetBranchAddress("passHLT0", &passHLT0, &b_passHLT0);
   fChain->SetBranchAddress("passHLT1", &passHLT1, &b_passHLT1);
   fChain->SetBranchAddress("passHLT2", &passHLT2, &b_passHLT2);
   fChain->SetBranchAddress("passHLT3", &passHLT3, &b_passHLT3);
   fChain->SetBranchAddress("passHLT4", &passHLT4, &b_passHLT4);
   fChain->SetBranchAddress("passHLT5", &passHLT5, &b_passHLT5);
   fChain->SetBranchAddress("passHLT6", &passHLT6, &b_passHLT6);
   fChain->SetBranchAddress("passHLTPrescaler0", &passHLTPrescaler0, &b_passHLTPrescaler0);
   fChain->SetBranchAddress("passHLTPrescaler1", &passHLTPrescaler1, &b_passHLTPrescaler1);
   fChain->SetBranchAddress("passHLTPrescaler2", &passHLTPrescaler2, &b_passHLTPrescaler2);
   fChain->SetBranchAddress("passHLTPrescaler3", &passHLTPrescaler3, &b_passHLTPrescaler3);
   fChain->SetBranchAddress("passHLTPrescaler4", &passHLTPrescaler4, &b_passHLTPrescaler4);
   fChain->SetBranchAddress("passHLTPrescaler5", &passHLTPrescaler5, &b_passHLTPrescaler5);
   fChain->SetBranchAddress("passHLTPrescaler6", &passHLTPrescaler6, &b_passHLTPrescaler6);
   fChain->SetBranchAddress("passL10", &passL10, &b_passL10);
   fChain->SetBranchAddress("passL11", &passL11, &b_passL11);
   fChain->SetBranchAddress("passL12", &passL12, &b_passL12);
   fChain->SetBranchAddress("passL13", &passL13, &b_passL13);
   fChain->SetBranchAddress("passL14", &passL14, &b_passL14);
   fChain->SetBranchAddress("passL15", &passL15, &b_passL15);
   fChain->SetBranchAddress("passL16", &passL16, &b_passL16);
   fChain->SetBranchAddress("passL1Prescaler0", &passL1Prescaler0, &b_passL1Prescaler0);
   fChain->SetBranchAddress("passL1Prescaler1", &passL1Prescaler1, &b_passL1Prescaler1);
   fChain->SetBranchAddress("passL1Prescaler2", &passL1Prescaler2, &b_passL1Prescaler2);
   fChain->SetBranchAddress("passL1Prescaler3", &passL1Prescaler3, &b_passL1Prescaler3);
   fChain->SetBranchAddress("passL1Prescaler4", &passL1Prescaler4, &b_passL1Prescaler4);
   fChain->SetBranchAddress("passL1Prescaler5", &passL1Prescaler5, &b_passL1Prescaler5);
   fChain->SetBranchAddress("passL1Prescaler6", &passL1Prescaler6, &b_passL1Prescaler6);
   fChain->SetBranchAddress("trk_dau0_isHP", &trk_dau0_isHP, &b_trk_dau0_isHP);
   fChain->SetBranchAddress("trk_dau1_isHP", &trk_dau1_isHP, &b_trk_dau1_isHP);
   fChain->SetBranchAddress("trk_gdau0_isHP", &trk_gdau0_isHP, &b_trk_gdau0_isHP);
   fChain->SetBranchAddress("trk_gdau1_isHP", &trk_gdau1_isHP, &b_trk_gdau1_isHP);
   fChain->SetBranchAddress("trk_isHP", &trk_isHP, &b_trk_isHP);
   fChain->SetBranchAddress("validPrescale0", &validPrescale0, &b_validPrescale0);
   fChain->SetBranchAddress("validPrescale1", &validPrescale1, &b_validPrescale1);
   fChain->SetBranchAddress("validPrescale2", &validPrescale2, &b_validPrescale2);
   fChain->SetBranchAddress("validPrescale3", &validPrescale3, &b_validPrescale3);
   fChain->SetBranchAddress("validPrescale4", &validPrescale4, &b_validPrescale4);
   fChain->SetBranchAddress("validPrescale5", &validPrescale5, &b_validPrescale5);
   fChain->SetBranchAddress("validPrescale6", &validPrescale6, &b_validPrescale6);
   fChain->SetBranchAddress("cand_dau0_status", &cand_dau0_status, &b_cand_dau0_status);
   fChain->SetBranchAddress("cand_dau1_status", &cand_dau1_status, &b_cand_dau1_status);
   fChain->SetBranchAddress("cand_gdau0_status", &cand_gdau0_status, &b_cand_gdau0_status);
   fChain->SetBranchAddress("cand_gdau1_status", &cand_gdau1_status, &b_cand_gdau1_status);
   fChain->SetBranchAddress("cand_status", &cand_status, &b_cand_status);
   fChain->SetBranchAddress("hltPDs0", &hltPDs0, &b_hltPDs0);
   fChain->SetBranchAddress("hltPDs1", &hltPDs1, &b_hltPDs1);
   fChain->SetBranchAddress("hltPDs2", &hltPDs2, &b_hltPDs2);
   fChain->SetBranchAddress("hltPDs3", &hltPDs3, &b_hltPDs3);
   fChain->SetBranchAddress("hltPDs4", &hltPDs4, &b_hltPDs4);
   fChain->SetBranchAddress("hltPDs5", &hltPDs5, &b_hltPDs5);
   fChain->SetBranchAddress("hltPDs6", &hltPDs6, &b_hltPDs6);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("cand_charge", &cand_charge, &b_cand_charge);
   fChain->SetBranchAddress("cand_dau0_charge", &cand_dau0_charge, &b_cand_dau0_charge);
   fChain->SetBranchAddress("cand_dau1_charge", &cand_dau1_charge, &b_cand_dau1_charge);
   fChain->SetBranchAddress("cand_gdau0_charge", &cand_gdau0_charge, &b_cand_gdau0_charge);
   fChain->SetBranchAddress("cand_gdau1_charge", &cand_gdau1_charge, &b_cand_gdau1_charge);
   fChain->SetBranchAddress("BXNb", &BXNb, &b_BXNb);
   fChain->SetBranchAddress("hltPrescale0", &hltPrescale0, &b_hltPrescale0);
   fChain->SetBranchAddress("hltPrescale1", &hltPrescale1, &b_hltPrescale1);
   fChain->SetBranchAddress("hltPrescale2", &hltPrescale2, &b_hltPrescale2);
   fChain->SetBranchAddress("hltPrescale3", &hltPrescale3, &b_hltPrescale3);
   fChain->SetBranchAddress("hltPrescale4", &hltPrescale4, &b_hltPrescale4);
   fChain->SetBranchAddress("hltPrescale5", &hltPrescale5, &b_hltPrescale5);
   fChain->SetBranchAddress("hltPrescale6", &hltPrescale6, &b_hltPrescale6);
   fChain->SetBranchAddress("l1Prescale0", &l1Prescale0, &b_l1Prescale0);
   fChain->SetBranchAddress("l1Prescale1", &l1Prescale1, &b_l1Prescale1);
   fChain->SetBranchAddress("l1Prescale2", &l1Prescale2, &b_l1Prescale2);
   fChain->SetBranchAddress("l1Prescale3", &l1Prescale3, &b_l1Prescale3);
   fChain->SetBranchAddress("l1Prescale4", &l1Prescale4, &b_l1Prescale4);
   fChain->SetBranchAddress("l1Prescale5", &l1Prescale5, &b_l1Prescale5);
   fChain->SetBranchAddress("l1Prescale6", &l1Prescale6, &b_l1Prescale6);
   fChain->SetBranchAddress("trk_dau0_nHit", &trk_dau0_nHit, &b_trk_dau0_nHit);
   fChain->SetBranchAddress("trk_dau1_nHit", &trk_dau1_nHit, &b_trk_dau1_nHit);
   fChain->SetBranchAddress("trk_gdau0_nHit", &trk_gdau0_nHit, &b_trk_gdau0_nHit);
   fChain->SetBranchAddress("trk_gdau1_nHit", &trk_gdau1_nHit, &b_trk_gdau1_nHit);
   fChain->SetBranchAddress("trk_nHit", &trk_nHit, &b_trk_nHit);
   fChain->SetBranchAddress("EventNb", &EventNb, &b_EventNb);
   fChain->SetBranchAddress("LSNb", &LSNb, &b_LSNb);
   fChain->SetBranchAddress("RunNb", &RunNb, &b_RunNb);
   fChain->SetBranchAddress("cand_dau0_pdgId", &cand_dau0_pdgId, &b_cand_dau0_pdgId);
   fChain->SetBranchAddress("cand_dau1_pdgId", &cand_dau1_pdgId, &b_cand_dau1_pdgId);
   fChain->SetBranchAddress("cand_gdau0_pdgId", &cand_gdau0_pdgId, &b_cand_gdau0_pdgId);
   fChain->SetBranchAddress("cand_gdau1_pdgId", &cand_gdau1_pdgId, &b_cand_gdau1_pdgId);
   fChain->SetBranchAddress("cand_pdgId", &cand_pdgId, &b_cand_pdgId);
   fChain->SetBranchAddress("HFsumETMinus", &HFsumETMinus, &b_HFsumETMinus);
   fChain->SetBranchAddress("HFsumETPlus", &HFsumETPlus, &b_HFsumETPlus);
   fChain->SetBranchAddress("Npixel", &Npixel, &b_Npixel);
   fChain->SetBranchAddress("Ntrkoffline", &Ntrkoffline, &b_Ntrkoffline);
   fChain->SetBranchAddress("ZDCMinus", &ZDCMinus, &b_ZDCMinus);
   fChain->SetBranchAddress("ZDCPlus", &ZDCPlus, &b_ZDCPlus);
   fChain->SetBranchAddress("bestvtxX", &bestvtxX, &b_bestvtxX);
   fChain->SetBranchAddress("bestvtxY", &bestvtxY, &b_bestvtxY);
   fChain->SetBranchAddress("bestvtxZ", &bestvtxZ, &b_bestvtxZ);
   fChain->SetBranchAddress("cand_angle2D", &cand_angle2D, &b_cand_angle2D);
   fChain->SetBranchAddress("cand_angle3D", &cand_angle3D, &b_cand_angle3D);
   fChain->SetBranchAddress("cand_dau0_angle2D", &cand_dau0_angle2D, &b_cand_dau0_angle2D);
   fChain->SetBranchAddress("cand_dau0_angle3D", &cand_dau0_angle3D, &b_cand_dau0_angle3D);
   fChain->SetBranchAddress("cand_dau0_dca", &cand_dau0_dca, &b_cand_dau0_dca);
   fChain->SetBranchAddress("cand_dau0_decayLength2D", &cand_dau0_decayLength2D, &b_cand_dau0_decayLength2D);
   fChain->SetBranchAddress("cand_dau0_decayLength3D", &cand_dau0_decayLength3D, &b_cand_dau0_decayLength3D);
   fChain->SetBranchAddress("cand_dau0_decayLengthError2D", &cand_dau0_decayLengthError2D, &b_cand_dau0_decayLengthError2D);
   fChain->SetBranchAddress("cand_dau0_decayLengthError3D", &cand_dau0_decayLengthError3D, &b_cand_dau0_decayLengthError3D);
   fChain->SetBranchAddress("cand_dau0_eta", &cand_dau0_eta, &b_cand_dau0_eta);
   fChain->SetBranchAddress("cand_dau0_etaDau0", &cand_dau0_etaDau0, &b_cand_dau0_etaDau0);
   fChain->SetBranchAddress("cand_dau0_etaDau1", &cand_dau0_etaDau1, &b_cand_dau0_etaDau1);
   fChain->SetBranchAddress("cand_dau0_mass", &cand_dau0_mass, &b_cand_dau0_mass);
   fChain->SetBranchAddress("cand_dau0_massDau0", &cand_dau0_massDau0, &b_cand_dau0_massDau0);
   fChain->SetBranchAddress("cand_dau0_massDau1", &cand_dau0_massDau1, &b_cand_dau0_massDau1);
   fChain->SetBranchAddress("cand_dau0_p", &cand_dau0_p, &b_cand_dau0_p);
   fChain->SetBranchAddress("cand_dau0_pT", &cand_dau0_pT, &b_cand_dau0_pT);
   fChain->SetBranchAddress("cand_dau0_pTDau0", &cand_dau0_pTDau0, &b_cand_dau0_pTDau0);
   fChain->SetBranchAddress("cand_dau0_pTDau1", &cand_dau0_pTDau1, &b_cand_dau0_pTDau1);
   fChain->SetBranchAddress("cand_dau0_phi", &cand_dau0_phi, &b_cand_dau0_phi);
   fChain->SetBranchAddress("cand_dau0_phiDau0", &cand_dau0_phiDau0, &b_cand_dau0_phiDau0);
   fChain->SetBranchAddress("cand_dau0_phiDau1", &cand_dau0_phiDau1, &b_cand_dau0_phiDau1);
   fChain->SetBranchAddress("cand_dau0_pseudoDecayLengthError2D", &cand_dau0_pseudoDecayLengthError2D, &b_cand_dau0_pseudoDecayLengthError2D);
   fChain->SetBranchAddress("cand_dau0_pseudoDecayLengthError3D", &cand_dau0_pseudoDecayLengthError3D, &b_cand_dau0_pseudoDecayLengthError3D);
   fChain->SetBranchAddress("cand_dau0_vtxChi2", &cand_dau0_vtxChi2, &b_cand_dau0_vtxChi2);
   fChain->SetBranchAddress("cand_dau0_vtxProb", &cand_dau0_vtxProb, &b_cand_dau0_vtxProb);
   fChain->SetBranchAddress("cand_dau0_y", &cand_dau0_y, &b_cand_dau0_y);
   fChain->SetBranchAddress("cand_dau1_angle2D", &cand_dau1_angle2D, &b_cand_dau1_angle2D);
   fChain->SetBranchAddress("cand_dau1_angle3D", &cand_dau1_angle3D, &b_cand_dau1_angle3D);
   fChain->SetBranchAddress("cand_dau1_dca", &cand_dau1_dca, &b_cand_dau1_dca);
   fChain->SetBranchAddress("cand_dau1_decayLength2D", &cand_dau1_decayLength2D, &b_cand_dau1_decayLength2D);
   fChain->SetBranchAddress("cand_dau1_decayLength3D", &cand_dau1_decayLength3D, &b_cand_dau1_decayLength3D);
   fChain->SetBranchAddress("cand_dau1_decayLengthError2D", &cand_dau1_decayLengthError2D, &b_cand_dau1_decayLengthError2D);
   fChain->SetBranchAddress("cand_dau1_decayLengthError3D", &cand_dau1_decayLengthError3D, &b_cand_dau1_decayLengthError3D);
   fChain->SetBranchAddress("cand_dau1_eta", &cand_dau1_eta, &b_cand_dau1_eta);
   fChain->SetBranchAddress("cand_dau1_mass", &cand_dau1_mass, &b_cand_dau1_mass);
   fChain->SetBranchAddress("cand_dau1_p", &cand_dau1_p, &b_cand_dau1_p);
   fChain->SetBranchAddress("cand_dau1_pT", &cand_dau1_pT, &b_cand_dau1_pT);
   fChain->SetBranchAddress("cand_dau1_phi", &cand_dau1_phi, &b_cand_dau1_phi);
   fChain->SetBranchAddress("cand_dau1_pseudoDecayLengthError2D", &cand_dau1_pseudoDecayLengthError2D, &b_cand_dau1_pseudoDecayLengthError2D);
   fChain->SetBranchAddress("cand_dau1_pseudoDecayLengthError3D", &cand_dau1_pseudoDecayLengthError3D, &b_cand_dau1_pseudoDecayLengthError3D);
   fChain->SetBranchAddress("cand_dau1_vtxChi2", &cand_dau1_vtxChi2, &b_cand_dau1_vtxChi2);
   fChain->SetBranchAddress("cand_dau1_vtxProb", &cand_dau1_vtxProb, &b_cand_dau1_vtxProb);
   fChain->SetBranchAddress("cand_dau1_y", &cand_dau1_y, &b_cand_dau1_y);
   fChain->SetBranchAddress("cand_dca", &cand_dca, &b_cand_dca);
   fChain->SetBranchAddress("cand_decayLength2D", &cand_decayLength2D, &b_cand_decayLength2D);
   fChain->SetBranchAddress("cand_decayLength3D", &cand_decayLength3D, &b_cand_decayLength3D);
   fChain->SetBranchAddress("cand_decayLengthError2D", &cand_decayLengthError2D, &b_cand_decayLengthError2D);
   fChain->SetBranchAddress("cand_decayLengthError3D", &cand_decayLengthError3D, &b_cand_decayLengthError3D);
   fChain->SetBranchAddress("cand_eta", &cand_eta, &b_cand_eta);
   fChain->SetBranchAddress("cand_etaDau0", &cand_etaDau0, &b_cand_etaDau0);
   fChain->SetBranchAddress("cand_etaDau1", &cand_etaDau1, &b_cand_etaDau1);
   fChain->SetBranchAddress("cand_gdau0_angle2D", &cand_gdau0_angle2D, &b_cand_gdau0_angle2D);
   fChain->SetBranchAddress("cand_gdau0_angle3D", &cand_gdau0_angle3D, &b_cand_gdau0_angle3D);
   fChain->SetBranchAddress("cand_gdau0_dca", &cand_gdau0_dca, &b_cand_gdau0_dca);
   fChain->SetBranchAddress("cand_gdau0_decayLength2D", &cand_gdau0_decayLength2D, &b_cand_gdau0_decayLength2D);
   fChain->SetBranchAddress("cand_gdau0_decayLength3D", &cand_gdau0_decayLength3D, &b_cand_gdau0_decayLength3D);
   fChain->SetBranchAddress("cand_gdau0_decayLengthError2D", &cand_gdau0_decayLengthError2D, &b_cand_gdau0_decayLengthError2D);
   fChain->SetBranchAddress("cand_gdau0_decayLengthError3D", &cand_gdau0_decayLengthError3D, &b_cand_gdau0_decayLengthError3D);
   fChain->SetBranchAddress("cand_gdau0_eta", &cand_gdau0_eta, &b_cand_gdau0_eta);
   fChain->SetBranchAddress("cand_gdau0_mass", &cand_gdau0_mass, &b_cand_gdau0_mass);
   fChain->SetBranchAddress("cand_gdau0_p", &cand_gdau0_p, &b_cand_gdau0_p);
   fChain->SetBranchAddress("cand_gdau0_pT", &cand_gdau0_pT, &b_cand_gdau0_pT);
   fChain->SetBranchAddress("cand_gdau0_phi", &cand_gdau0_phi, &b_cand_gdau0_phi);
   fChain->SetBranchAddress("cand_gdau0_pseudoDecayLengthError2D", &cand_gdau0_pseudoDecayLengthError2D, &b_cand_gdau0_pseudoDecayLengthError2D);
   fChain->SetBranchAddress("cand_gdau0_pseudoDecayLengthError3D", &cand_gdau0_pseudoDecayLengthError3D, &b_cand_gdau0_pseudoDecayLengthError3D);
   fChain->SetBranchAddress("cand_gdau0_vtxChi2", &cand_gdau0_vtxChi2, &b_cand_gdau0_vtxChi2);
   fChain->SetBranchAddress("cand_gdau0_vtxProb", &cand_gdau0_vtxProb, &b_cand_gdau0_vtxProb);
   fChain->SetBranchAddress("cand_gdau0_y", &cand_gdau0_y, &b_cand_gdau0_y);
   fChain->SetBranchAddress("cand_gdau1_angle2D", &cand_gdau1_angle2D, &b_cand_gdau1_angle2D);
   fChain->SetBranchAddress("cand_gdau1_angle3D", &cand_gdau1_angle3D, &b_cand_gdau1_angle3D);
   fChain->SetBranchAddress("cand_gdau1_dca", &cand_gdau1_dca, &b_cand_gdau1_dca);
   fChain->SetBranchAddress("cand_gdau1_decayLength2D", &cand_gdau1_decayLength2D, &b_cand_gdau1_decayLength2D);
   fChain->SetBranchAddress("cand_gdau1_decayLength3D", &cand_gdau1_decayLength3D, &b_cand_gdau1_decayLength3D);
   fChain->SetBranchAddress("cand_gdau1_decayLengthError2D", &cand_gdau1_decayLengthError2D, &b_cand_gdau1_decayLengthError2D);
   fChain->SetBranchAddress("cand_gdau1_decayLengthError3D", &cand_gdau1_decayLengthError3D, &b_cand_gdau1_decayLengthError3D);
   fChain->SetBranchAddress("cand_gdau1_eta", &cand_gdau1_eta, &b_cand_gdau1_eta);
   fChain->SetBranchAddress("cand_gdau1_mass", &cand_gdau1_mass, &b_cand_gdau1_mass);
   fChain->SetBranchAddress("cand_gdau1_p", &cand_gdau1_p, &b_cand_gdau1_p);
   fChain->SetBranchAddress("cand_gdau1_pT", &cand_gdau1_pT, &b_cand_gdau1_pT);
   fChain->SetBranchAddress("cand_gdau1_phi", &cand_gdau1_phi, &b_cand_gdau1_phi);
   fChain->SetBranchAddress("cand_gdau1_pseudoDecayLengthError2D", &cand_gdau1_pseudoDecayLengthError2D, &b_cand_gdau1_pseudoDecayLengthError2D);
   fChain->SetBranchAddress("cand_gdau1_pseudoDecayLengthError3D", &cand_gdau1_pseudoDecayLengthError3D, &b_cand_gdau1_pseudoDecayLengthError3D);
   fChain->SetBranchAddress("cand_gdau1_vtxChi2", &cand_gdau1_vtxChi2, &b_cand_gdau1_vtxChi2);
   fChain->SetBranchAddress("cand_gdau1_vtxProb", &cand_gdau1_vtxProb, &b_cand_gdau1_vtxProb);
   fChain->SetBranchAddress("cand_gdau1_y", &cand_gdau1_y, &b_cand_gdau1_y);
   fChain->SetBranchAddress("cand_mass", &cand_mass, &b_cand_mass);
   fChain->SetBranchAddress("cand_massDau0", &cand_massDau0, &b_cand_massDau0);
   fChain->SetBranchAddress("cand_massDau1", &cand_massDau1, &b_cand_massDau1);
   fChain->SetBranchAddress("cand_p", &cand_p, &b_cand_p);
   fChain->SetBranchAddress("cand_pT", &cand_pT, &b_cand_pT);
   fChain->SetBranchAddress("cand_pTDau0", &cand_pTDau0, &b_cand_pTDau0);
   fChain->SetBranchAddress("cand_pTDau1", &cand_pTDau1, &b_cand_pTDau1);
   fChain->SetBranchAddress("cand_phi", &cand_phi, &b_cand_phi);
   fChain->SetBranchAddress("cand_phiDau0", &cand_phiDau0, &b_cand_phiDau0);
   fChain->SetBranchAddress("cand_phiDau1", &cand_phiDau1, &b_cand_phiDau1);
   fChain->SetBranchAddress("cand_pseudoDecayLengthError2D", &cand_pseudoDecayLengthError2D, &b_cand_pseudoDecayLengthError2D);
   fChain->SetBranchAddress("cand_pseudoDecayLengthError3D", &cand_pseudoDecayLengthError3D, &b_cand_pseudoDecayLengthError3D);
   fChain->SetBranchAddress("cand_vtxChi2", &cand_vtxChi2, &b_cand_vtxChi2);
   fChain->SetBranchAddress("cand_vtxProb", &cand_vtxProb, &b_cand_vtxProb);
   fChain->SetBranchAddress("cand_y", &cand_y, &b_cand_y);
   fChain->SetBranchAddress("pileup", &pileup, &b_pileup);
   fChain->SetBranchAddress("rawInstLumi", &rawInstLumi, &b_rawInstLumi);
   fChain->SetBranchAddress("trk_dau0_nChi2", &trk_dau0_nChi2, &b_trk_dau0_nChi2);
   fChain->SetBranchAddress("trk_dau0_pTErr", &trk_dau0_pTErr, &b_trk_dau0_pTErr);
   fChain->SetBranchAddress("trk_dau0_xyDCASignificance", &trk_dau0_xyDCASignificance, &b_trk_dau0_xyDCASignificance);
   fChain->SetBranchAddress("trk_dau0_zDCASignificance", &trk_dau0_zDCASignificance, &b_trk_dau0_zDCASignificance);
   fChain->SetBranchAddress("trk_dau1_nChi2", &trk_dau1_nChi2, &b_trk_dau1_nChi2);
   fChain->SetBranchAddress("trk_dau1_pTErr", &trk_dau1_pTErr, &b_trk_dau1_pTErr);
   fChain->SetBranchAddress("trk_dau1_xyDCASignificance", &trk_dau1_xyDCASignificance, &b_trk_dau1_xyDCASignificance);
   fChain->SetBranchAddress("trk_dau1_zDCASignificance", &trk_dau1_zDCASignificance, &b_trk_dau1_zDCASignificance);
   fChain->SetBranchAddress("trk_gdau0_nChi2", &trk_gdau0_nChi2, &b_trk_gdau0_nChi2);
   fChain->SetBranchAddress("trk_gdau0_pTErr", &trk_gdau0_pTErr, &b_trk_gdau0_pTErr);
   fChain->SetBranchAddress("trk_gdau0_xyDCASignificance", &trk_gdau0_xyDCASignificance, &b_trk_gdau0_xyDCASignificance);
   fChain->SetBranchAddress("trk_gdau0_zDCASignificance", &trk_gdau0_zDCASignificance, &b_trk_gdau0_zDCASignificance);
   fChain->SetBranchAddress("trk_gdau1_nChi2", &trk_gdau1_nChi2, &b_trk_gdau1_nChi2);
   fChain->SetBranchAddress("trk_gdau1_pTErr", &trk_gdau1_pTErr, &b_trk_gdau1_pTErr);
   fChain->SetBranchAddress("trk_gdau1_xyDCASignificance", &trk_gdau1_xyDCASignificance, &b_trk_gdau1_xyDCASignificance);
   fChain->SetBranchAddress("trk_gdau1_zDCASignificance", &trk_gdau1_zDCASignificance, &b_trk_gdau1_zDCASignificance);
   fChain->SetBranchAddress("trk_nChi2", &trk_nChi2, &b_trk_nChi2);
   fChain->SetBranchAddress("trk_pTErr", &trk_pTErr, &b_trk_pTErr);
   fChain->SetBranchAddress("trk_xyDCASignificance", &trk_xyDCASignificance, &b_trk_xyDCASignificance);
   fChain->SetBranchAddress("trk_zDCASignificance", &trk_zDCASignificance, &b_trk_zDCASignificance);
   Notify();
}

Bool_t ParticleNTuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ParticleNTuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ParticleNTuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ParticleNTuple_cxx
