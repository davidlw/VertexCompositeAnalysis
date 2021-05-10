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

void test()
{
  using PtEtaPhiM_t = ROOT::Math::PtEtaPhiM4D<double>;

  std::unique_ptr<TH2D> hKsRecoAll(new TH2D("hKsRecoAll", "All channels, Ks p3;pT (GeV);y;",
          40, 0, 8, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hKsRecoAll_mass(new TH2D("hKsRecoAll_mass", "All channels, Ks p3;pT (GeV);mass (GeV);",
          40, 0, 8, 80, 0.44, 0.55));

  std::unique_ptr<TH2D> hKsRecoDaus(new TH2D("hKsRecoDaus", "All daughters, Ks p3;pT (GeV);pT (GeV;",
          40, 0, 8, 40, 0, 8));

  std::unique_ptr<TH2D> hLamRecoAll(new TH2D("hLamRecoAll", "All channels, Lam p3;pT (GeV);y;",
          40, 0, 8, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hLamRecoAll_mass(new TH2D("hLamRecoAll_mass", "All channels, Lam p3;pT (GeV);mass (GeV);",
          40, 0, 8, 80, 1.08, 1.16));

  std::unique_ptr<TH2D> hLamRecoDaus(new TH2D("hLamRecoDaus", "All daughters, Lam p3;pT (GeV);pT (GeV;",
          40, 0, 8, 40, 0, 8));

  std::unique_ptr<TH2D> hXiRecoAll(new TH2D("hXiRecoAll", "All channels, Xi p3;pT (GeV);y;",
          40, 0, 8, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hXiRecoAll_mass(new TH2D("hXiRecoAll_mass", "All channels, Xi p3;pT (GeV);mass (GeV);",
          40, 0, 8, 80, 1.26, 1.39));

  std::unique_ptr<TH2D> hOmRecoAll(new TH2D("hOmRecoAll", "All channels, Om p3;pT (GeV);y;",
          40, 0, 8, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hOmRecoAll_mass(new TH2D("hOmRecoAll_mass", "All channels, Om p3;pT (GeV);mass (GeV);",
          40, 0, 8, 80, 1.61, 1.73));

  const TString& inputList = "filelist.txt";
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

/*
    auto charge = p.cand_charge();
    auto gen_pdgId = p.gen_pdgId();
    auto gen_charge = p.gen_charge();

    auto dauIdxEvt = p.cand_dauIdx();
    auto gen_dauIdxEvt = p.gen_dauIdx();

cout<<"get here"<<endl;
*/
    for (size_t ireco=0; ireco<recosize_ks; ireco++) {
       if( fabs(pdgId_ks[ireco]) != 310 ) continue;
       hKsRecoAll->Fill(p_ks.cand_pT()[ireco], p_ks.cand_y()[ireco]);
       hKsRecoAll_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);

//       auto dauIdx = dauIdxEvt_ks.at(ireco);
// cout<<dauIdx[0]<<" "<<dauIdx[1]<<endl;

//       hKsRecoDaus->Fill(p_ks.cand_pT()[dauIdx[0]],p_ks.cand_pT()[dauIdx[1]]);
    }
      
    for (size_t ireco=0; ireco<recosize_lam; ireco++) {
       if( fabs(pdgId_lam[ireco]) != 3122 ) continue;
       hLamRecoAll->Fill(p_lam.cand_pT()[ireco], p_lam.cand_y()[ireco]);
       hLamRecoAll_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);

//       auto dauIdx = dauIdxEvt_lam.at(ireco);
//       hLamRecoDaus->Fill(p_lam.cand_pT()[dauIdx[0]],p_lam.cand_pT()[dauIdx[1]]);
    }

    for (size_t ireco=0; ireco<recosize_xi; ireco++) {
       if( fabs(pdgId_xi[ireco]) != 3312 ) continue;
       hXiRecoAll->Fill(p_xi.cand_pT()[ireco], p_xi.cand_y()[ireco]);
       hXiRecoAll_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
    }

    for (size_t ireco=0; ireco<recosize_om; ireco++) {
       if( fabs(pdgId_om[ireco]) != 3334 ) continue;
       hOmRecoAll->Fill(p_om.cand_pT()[ireco], p_om.cand_y()[ireco]);
       hOmRecoAll_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
    }

/*
    for (size_t ireco=0; ireco<recosize; ireco++) {
      // begin Ks
      if (pdgId[ireco] == 310) {
        hKsRecoAll->Fill(p.cand_pT()[ireco], p.cand_y()[ireco]);
        hKsRecoAll_mass->Fill(p.cand_pT()[ireco], p.cand_mass()[ireco]);
      }
      if (pdgId[ireco] == 3122) {
        hLamRecoAll->Fill(p.cand_pT()[ireco], p.cand_y()[ireco]);
        hLamRecoAll_mass->Fill(p.cand_pT()[ireco], p.cand_mass()[ireco]);
      }
      if (pdgId[ireco] == 3312) {
        hXiRecoAll->Fill(p.cand_pT()[ireco], p.cand_y()[ireco]);
        hXiRecoAll_mass->Fill(p.cand_pT()[ireco], p.cand_mass()[ireco]);
      }
      if (pdgId[ireco] == 3334) {
        hOmRecoAll->Fill(p.cand_pT()[ireco], p.cand_y()[ireco]);
        hOmRecoAll_mass->Fill(p.cand_pT()[ireco], p.cand_mass()[ireco]);
      }
*/
	// reco daughters
/* 
	bool matchGEN = false;
	bool breakMatchGEN = false;
	bool breakMatch0Pi = false;
	bool breakMatch1Pi = false;
	bool breakMatch2Pi = false;

	auto dauIdx = dauIdxEvt.at(ireco);

	PtEtaPhiM_t recoPi0(
			    p.cand_pT()[dauIdx[0]],
			    p.cand_eta()[dauIdx[0]],
			    p.cand_phi()[dauIdx[0]],
			    p.cand_mass()[dauIdx[0]]
			    );
	PtEtaPhiM_t recoPi1(
			    p.cand_pT()[dauIdx[1]],
			    p.cand_eta()[dauIdx[1]],
			    p.cand_phi()[dauIdx[1]],
			    p.cand_mass()[dauIdx[1]]
			    );
	auto chargePi0 = charge[dauIdx[0]];
	auto chargePi1 = charge[dauIdx[1]];
	auto pdgIdPi0 = pdgId[dauIdx[0]];
	auto pdgIdPi1 = pdgId[dauIdx[1]];
	//cout << "pdgIdPi0: " << pdgIdPi0 << endl;
	//cout << "pdgIdPi1: " << pdgIdPi1 << endl;
	//cout << "pdgId Pi0 " << chargePi0 * pdgIdPi0 << endl;
	//cout << "pdgId Pi1 " << chargePi1 * pdgIdPi1 << endl;
	PtEtaPhiM_t recoKs (
			    p.cand_pT()[ireco],
			    p.cand_eta()[ireco],
			    p.cand_phi()[ireco],
			    p.cand_mass()[ireco]
			    );
*/
//        cout<<p.cand_pT()[ireco]<<" "<<p.cand_eta()[ireco]<<" "<<p.cand_phi()[ireco]<<" "<<p.cand_mass()[ireco]<<endl;

/*
	// reco-gen matching
	for (size_t igen=0; igen<gensize; igen++) {
			     
			     p.gen_pT()[igen],
			     p.gen_eta()[igen],
			     p.gen_phi()[igen],
			     p.gen_mass()[igen]

          hKsRecoMatched->Fill(p.cand_pT()[ireco], p.cand_y()[ireco]);
          hKsRecoMatched_mass->Fill(p.cand_pT()[ireco], p.cand_mass()[ireco]);

      } // end Ks
    }
*/

  }

  TFile* fout = new TFile("output.root","recreate");

  hKsRecoAll->Write();
  hKsRecoAll_mass->Write();
  hKsRecoDaus->Write();

  hLamRecoAll->Write();
  hLamRecoAll_mass->Write();
  hLamRecoDaus->Write();

  hXiRecoAll->Write();
  hXiRecoAll_mass->Write();

  hOmRecoAll->Write();
  hOmRecoAll_mass->Write();

  fout->Close();
}
