#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

void copyTree()
{
  TChain* chain = new TChain("jpsiana_wrongsign/VertexCompositeNtuple");
  chain->Add("./training_samples/Merged_PbpData_JPsi_TrainingTree_HLT185_all_v3.root");
//  chain->Add("/mnt/hadoop/cms/store/user/davidlw/PAHighMultiplicity1/pPb_Tree_D0_training_run285505_HLT185_v7/180828_103545/0000/d0ana_training_2.root");
//  chain->Add("/mnt/hadoop/cms/store/user/davidlw/PAHighMultiplicity1/pPb_Tree_D0_training_run285505_HLT185_v7/180828_103545/0000/d0ana_training_3.root");
//  chain->Add("/mnt/hadoop/cms/store/user/davidlw/PAHighMultiplicity1/pPb_Tree_D0_training_run285505_HLT185_v7/180828_103545/0000/d0ana_training_4.root");
//  chain->Add("/mnt/hadoop/cms/store/user/davidlw/PAHighMultiplicity1/pPb_Tree_D0_training_run285505_HLT185_v7/180828_103545/0000/d0ana_training_5.root");
//  chain->Add("/mnt/hadoop/cms/store/user/davidlw/PAHighMultiplicity1/pPb_Tree_D0_training_run285505_HLT185_v7/180828_103545/0000/d0ana_training_6.root");
//  chain->Add("/mnt/hadoop/cms/store/user/davidlw/PAHighMultiplicity1/pPb_Tree_D0_training_run285505_HLT185_v7/180828_103545/0000/d0ana_training_7.root");
//  chain->Add("/mnt/hadoop/cms/store/user/davidlw/PAHighMultiplicity1/pPb_Tree_D0_training_run285505_HLT185_v7/180828_103545/0000/d0ana_training_8.root");

  TFile* outFile = TFile::Open("./training_samples/Merged_PbpData_JPsi_TrainingTree_WS_HLT185_v3.root","RECREATE");
  TTree* copyTree = chain->CopyTree("Ntrkoffline>=185", "");
  outFile->Write();
  outFile->Close();
}
