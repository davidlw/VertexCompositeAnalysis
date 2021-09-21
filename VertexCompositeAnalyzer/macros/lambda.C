#include <iostream>

#include "TreeReader/ParticleTreeMC.cxx"
#include "TreeReader/ParticleTree.cxx"
#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "THashList.h"
#include "TH1D.h"
#include "TH2D.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "Math/GenVector/VectorUtil.h"
#include "TMath.h"

using namespace std; 

void lambda()
{
  using PtEtaPhiM_t = ROOT::Math::PtEtaPhiM4D<double>;

  std::unique_ptr<TH2D> hKsRecoAll(new TH2D("hKsRecoAll", "All channels, Ks p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hKsRecoAll_mass(new TH2D("hKsRecoAll_mass", "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55));

  std::unique_ptr<TH2D> hLamRecoAll(new TH2D("hLamRecoAll", "All channels, Lam p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hLamRecoAll_mass(new TH2D("hLamRecoAll_mass", "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16));

  std::unique_ptr<TH2D> hXiRecoAll(new TH2D("hXiRecoAll", "All channels, Xi p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hXiRecoAll_mass(new TH2D("hXiRecoAll_mass", "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39));

  std::unique_ptr<TH2D> hOmRecoAll(new TH2D("hOmRecoAll", "All channels, Om p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hOmRecoAll_mass(new TH2D("hOmRecoAll_mass", "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73));

  const TString& inputList = "filelist_lambda.txt";
  const TString& treeDir_ks = "kshortana";
  const TString& treeDir_lam = "lambdaana";
  const TString& treeDir_xi = "xiana";
  const TString& treeDir_om = "omegaana";

  TFileCollection tf("tf", "", inputList);

  TChain t_ks(treeDir_ks+"/ParticleTree");
  t_ks.AddFileInfoList(tf.GetList());
  ParticleTree p_ks(&t_ks);

  TChain t_lam(treeDir_lam+"/ParticleTree");
  t_lam.AddFileInfoList(tf.GetList());
  ParticleTree p_lam(&t_lam);

  TChain t_xi(treeDir_xi+"/ParticleTree");
  t_xi.AddFileInfoList(tf.GetList());
  ParticleTree p_xi(&t_xi);

  TChain t_om(treeDir_om+"/ParticleTree");
  t_om.AddFileInfoList(tf.GetList());
  ParticleTree p_om(&t_om);

  auto nentries = p_ks.GetEntries();
  //auto nentries = 1000L;

  cout << "Tree lambdaana/ParticleTree in " << inputList
       << " has " << nentries << " entries." << endl;;

  for (Long64_t ientry=0; ientry<nentries; ientry++) {

    if( !(ientry % 10000) ) std::cout<<"Processing "<<ientry<<" events"<<std::endl;

    p_ks.GetEntry(ientry);
    p_lam.GetEntry(ientry);
    p_xi.GetEntry(ientry);
    p_om.GetEntry(ientry);

    auto recosize_ks = p_ks.cand_mass().size();
    auto recosize_lam = p_lam.cand_mass().size();
    auto recosize_xi = p_xi.cand_mass().size();
    auto recosize_om = p_om.cand_mass().size();

    //if(gensize == 0 || recosize == 0) continue;

    auto pdgId_ks = p_ks.cand_pdgId();
    auto pdgId_lam = p_lam.cand_pdgId();
    auto pdgId_xi = p_xi.cand_pdgId();
    auto pdgId_om = p_om.cand_pdgId();

    auto dauIdxEvt_ks = p_ks.cand_dauIdx();
    auto dauIdxEvt_lam = p_lam.cand_dauIdx();
    auto dauIdxEvt_xi = p_xi.cand_dauIdx();
    auto dauIdxEvt_om = p_om.cand_dauIdx();

    for (size_t ireco=0; ireco<recosize_ks; ireco++) {
       if( fabs(pdgId_ks[ireco]) != 310 ) continue;
       if(p_ks.cand_decayLength3D()[ireco]/p_ks.cand_decayLengthError3D()[ireco]<5) continue;
       if(cos(p_ks.cand_angle3D()[ireco])<0.999) continue;

       hKsRecoAll->Fill(p_ks.cand_pT()[ireco], p_ks.cand_y()[ireco]);
       hKsRecoAll_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);
    }
      
    for (size_t ireco=0; ireco<recosize_lam; ireco++) {
       if( fabs(pdgId_lam[ireco]) != 3122 ) continue;
       if(p_lam.cand_decayLength3D()[ireco]/p_lam.cand_decayLengthError3D()[ireco]<5) continue;
       if(cos(p_lam.cand_angle3D()[ireco])<0.999) continue;

       hLamRecoAll->Fill(p_lam.cand_pT()[ireco], p_lam.cand_y()[ireco]);
       hLamRecoAll_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
    }

    for (size_t ireco=0; ireco<recosize_xi; ireco++) {
       if( fabs(pdgId_xi[ireco]) != 3312 ) continue;
       if(p_xi.cand_decayLength3D()[ireco]/p_xi.cand_decayLengthError3D()[ireco]<3) continue;
       if(cos(p_xi.cand_angle3D()[ireco])<0.9998) continue;

       auto dauIdx = dauIdxEvt_xi.at(ireco);
       if(p_xi.cand_decayLength3D()[dauIdx[1]]/p_xi.cand_decayLengthError3D()[dauIdx[1]]<12) continue;

       auto trkIdx0 = p_xi.cand_trkIdx().at(dauIdx.at(0));
       if(p_xi.trk_xyDCASignificance().at(trkIdx0)<5) continue;
       if(p_xi.trk_zDCASignificance().at(trkIdx0)<5) continue;

       auto gdauIdx = dauIdxEvt_xi.at(dauIdx[1]);
       auto gtrkIdx0 = p_xi.cand_trkIdx().at(gdauIdx.at(0));
       auto gtrkIdx1 = p_xi.cand_trkIdx().at(gdauIdx.at(1));
       if(p_xi.trk_xyDCASignificance().at(gtrkIdx0)<3) continue;
       if(p_xi.trk_zDCASignificance().at(gtrkIdx0)<3) continue;
       if(p_xi.trk_xyDCASignificance().at(gtrkIdx1)<3) continue;
       if(p_xi.trk_zDCASignificance().at(gtrkIdx1)<3) continue;

       hXiRecoAll->Fill(p_xi.cand_pT()[ireco], p_xi.cand_y()[ireco]);
       hXiRecoAll_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
    }

    for (size_t ireco=0; ireco<recosize_om; ireco++) {
       if( fabs(pdgId_om[ireco]) != 3334 ) continue;
       if(p_om.cand_decayLength3D()[ireco]/p_om.cand_decayLengthError3D()[ireco]<2) continue;
       if(cos(p_om.cand_angle3D()[ireco])<0.9998) continue;

       auto dauIdx = dauIdxEvt_om.at(ireco);
       if(p_om.cand_decayLength3D()[dauIdx[1]]/p_om.cand_decayLengthError3D()[dauIdx[1]]<10) continue;

       auto trkIdx0 = p_om.cand_trkIdx().at(dauIdx.at(0));
       if(p_om.trk_xyDCASignificance().at(trkIdx0)<4) continue;
       if(p_om.trk_zDCASignificance().at(trkIdx0)<4) continue;

       auto gdauIdx = dauIdxEvt_om.at(dauIdx[1]);
       auto gtrkIdx0 = p_om.cand_trkIdx().at(gdauIdx.at(0));
       auto gtrkIdx1 = p_om.cand_trkIdx().at(gdauIdx.at(1));
       if(p_om.trk_xyDCASignificance().at(gtrkIdx0)<2) continue;
       if(p_om.trk_zDCASignificance().at(gtrkIdx0)<2) continue;
       if(p_om.trk_xyDCASignificance().at(gtrkIdx1)<2) continue;
       if(p_om.trk_zDCASignificance().at(gtrkIdx1)<2) continue;

       hOmRecoAll->Fill(p_om.cand_pT()[ireco], p_om.cand_y()[ireco]);
       hOmRecoAll_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
    }
  }

  TFile* fout = new TFile("lambda_output.root","recreate");

  hKsRecoAll->Write();
  hKsRecoAll_mass->Write();

  hLamRecoAll->Write();
  hLamRecoAll_mass->Write();

  hXiRecoAll->Write();
  hXiRecoAll_mass->Write();

  hOmRecoAll->Write();
  hOmRecoAll_mass->Write();

  fout->Close();
}
