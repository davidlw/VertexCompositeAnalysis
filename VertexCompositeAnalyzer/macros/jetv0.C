#include <iostream>

#include "TreeReader/ParticleTreeMC.cxx"
#include "TreeReader/ParticleTree.cxx"
#include "TreeReader/JetTrackTree.cxx"
#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "THashList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "Math/GenVector/VectorUtil.h"
#include "TMath.h"

using namespace std; 

double nMultBin[] = {0,5,10,20,40,60,80,1000};
int getNMultBin( double nMult)
{
  int bin=-1;
  for(int i=0;i<7;i++) 
    { if(nMult>=nMultBin[i] && nMult<nMultBin[i+1]) bin = i; }

  return bin;
}

void jetv0(const TString& inputList = "filelist_v0.txt")
{
  using PtEtaPhiM_t = ROOT::Math::PtEtaPhiM4D<double>;

  std::unique_ptr<TH2D> hKsRecoAll(new TH2D("hKsRecoAll", "All channels, Ks p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hKsRecoAll_mass(new TH2D("hKsRecoAll_mass", "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55));

  std::unique_ptr<TH2D> hJet4KsRecoAll_mass(new TH2D("hJet4KsRecoAll_mass", "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55));

  std::unique_ptr<TH2D> hJet8KsRecoAll_mass(new TH2D("hJet8KsRecoAll_mass", "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55));

  std::unique_ptr<TH2D> hLamRecoAll(new TH2D("hLamRecoAll", "All channels, Lam p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hLamRecoAll_mass(new TH2D("hLamRecoAll_mass", "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16));

  std::unique_ptr<TH2D> hJet4LamRecoAll_mass(new TH2D("hJet4LamRecoAll_mass", "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16));

  std::unique_ptr<TH2D> hJet8LamRecoAll_mass(new TH2D("hJet8LamRecoAll_mass", "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16));

  std::unique_ptr<TH2D> hXiRecoAll(new TH2D("hXiRecoAll", "All channels, Xi p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hXiRecoAll_mass(new TH2D("hXiRecoAll_mass", "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39));

  std::unique_ptr<TH2D> hJet4XiRecoAll_mass(new TH2D("hJet4XiRecoAll_mass", "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39));

  std::unique_ptr<TH2D> hJet8XiRecoAll_mass(new TH2D("hJet8XiRecoAll_mass", "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39));

  std::unique_ptr<TH2D> hOmRecoAll(new TH2D("hOmRecoAll", "All channels, Om p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hOmRecoAll_mass(new TH2D("hOmRecoAll_mass", "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73));

  std::unique_ptr<TH2D> hJet4OmRecoAll_mass(new TH2D("hJet4OmRecoAll_mass", "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73));

  std::unique_ptr<TH2D> hJet8OmRecoAll_mass(new TH2D("hJet8OmRecoAll_mass", "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73));

  std::unique_ptr<TH2D> hJetEtaPt(new TH2D("hJetEtaPt", ";pT (GeV);#eta;",
          100, 0, 1000, 10, -2.5, 2.5));

  std::unique_ptr<TH1D> hJetMult(new TH1D("hJetMult", ";N_{trk}", 400, 0, 400));

  TH2D* hJet4KsRecoMult_mass[100];
  TH2D* hJet8KsRecoMult_mass[100];
  TH2D* hJet4LamRecoMult_mass[100];
  TH2D* hJet8LamRecoMult_mass[100];
  TH2D* hJet4XiRecoMult_mass[100];
  TH2D* hJet8XiRecoMult_mass[100];
  TH2D* hJet4OmRecoMult_mass[100];
  TH2D* hJet8OmRecoMult_mass[100];
  double nMultBin[] = {0,5,10,20,40,60,80,200};

  for(int i=0;i<7;i++)
  {
    hJet4KsRecoMult_mass[i] = new TH2D(Form("hJet4KsRecoMult_mass_%d",i), "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55);
    hJet8KsRecoMult_mass[i] = new TH2D(Form("hJet8KsRecoMult_mass_%d",i), "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55);

    hJet4LamRecoMult_mass[i] = new TH2D(Form("hJet4LamRecoMult_mass_%d",i), "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16);
    hJet8LamRecoMult_mass[i] = new TH2D(Form("hJet8LamRecoMult_mass_%d",i), "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16);

    hJet4XiRecoMult_mass[i] = new TH2D(Form("hJet4XiRecoMult_mass_%d",i), "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39);
    hJet8XiRecoMult_mass[i] = new TH2D(Form("hJet8XiRecoMult_mass_%d",i), "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39);

    hJet4OmRecoMult_mass[i] = new TH2D(Form("hJet4OmRecoMult_mass_%d",i), "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73);
    hJet8OmRecoMult_mass[i] = new TH2D(Form("hJet8OmRecoMult_mass_%d",i), "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73);
  }
 
  const TString& treeDir_ks = "kshortana";
  const TString& treeDir_lam = "lambdaana";
  const TString& treeDir_xi = "xiana";
  const TString& treeDir_om = "omegaana";
  const TString& treeDir_antilam = "antilambdaana";
  const TString& treeDir_antixi = "antixiana";
  const TString& treeDir_antiom = "antiomegaana";
  const TString& treeDir_jet = "jetanalyzer";

  TFileCollection tf("tf", "", inputList);

  TChain t_ks(treeDir_ks+"/ParticleTree");
  t_ks.AddFileInfoList(tf.GetList());
  ParticleTree p_ks(&t_ks);

  TChain t_antilam(treeDir_antilam+"/ParticleTree");
  t_antilam.AddFileInfoList(tf.GetList());
  ParticleTree p_antilam(&t_antilam);

  TChain t_antixi(treeDir_antixi+"/ParticleTree");
  t_antixi.AddFileInfoList(tf.GetList());
  ParticleTree p_antixi(&t_antixi);

  TChain t_antiom(treeDir_antiom+"/ParticleTree");
  t_antiom.AddFileInfoList(tf.GetList());
  ParticleTree p_antiom(&t_antiom);

  TChain t_lam(treeDir_lam+"/ParticleTree");
  t_lam.AddFileInfoList(tf.GetList());
//  t_lam.Add(&t_antilam);
  ParticleTree p_lam(&t_lam);

  TChain t_xi(treeDir_xi+"/ParticleTree");
  t_xi.AddFileInfoList(tf.GetList());
//  t_xi.Add(&t_antixi);
  ParticleTree p_xi(&t_xi);

  TChain t_om(treeDir_om+"/ParticleTree");
  t_om.AddFileInfoList(tf.GetList());
//  t_om.Add(&t_antiom);
  ParticleTree p_om(&t_om);

  TChain t_jet(treeDir_jet+"/trackTree");
  t_jet.AddFileInfoList(tf.GetList());
  JetTrackTree p_jet(&t_jet);

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
    p_antilam.GetEntry(ientry);
    p_antixi.GetEntry(ientry);
    p_antiom.GetEntry(ientry);
    p_jet.GetEntry(ientry);

    auto recosize_ks = p_ks.cand_mass().size();
    auto recosize_lam = p_lam.cand_mass().size();
    auto recosize_xi = p_xi.cand_mass().size();
    auto recosize_om = p_om.cand_mass().size();
    auto recosize_antilam = p_antilam.cand_mass().size();
    auto recosize_antixi = p_antixi.cand_mass().size();
    auto recosize_antiom = p_antiom.cand_mass().size();
    auto recosize_jet = p_jet.jetPt().size();

    //if(gensize == 0 || recosize == 0) continue;

    auto pdgId_ks = p_ks.cand_pdgId();
    auto pdgId_lam = p_lam.cand_pdgId();
    auto pdgId_xi = p_xi.cand_pdgId();
    auto pdgId_om = p_om.cand_pdgId();
    auto pdgId_antilam = p_antilam.cand_pdgId();
    auto pdgId_antixi = p_antixi.cand_pdgId();
    auto pdgId_antiom = p_antiom.cand_pdgId();

    auto dauIdxEvt_ks = p_ks.cand_dauIdx();
    auto dauIdxEvt_lam = p_lam.cand_dauIdx();
    auto dauIdxEvt_xi = p_xi.cand_dauIdx();
    auto dauIdxEvt_om = p_om.cand_dauIdx();
    auto dauIdxEvt_antilam = p_antilam.cand_dauIdx();
    auto dauIdxEvt_antixi = p_antixi.cand_dauIdx();
    auto dauIdxEvt_antiom = p_antiom.cand_dauIdx();

    double jetpt=0;
    double jeteta=0;
    double jetphi=0; 
    int    jetmult=0;
    for (size_t ireco=0; ireco<recosize_jet; ireco++) { 
      double jetpttmp = p_jet.jetPt()[ireco];
      if(jetpttmp>jetpt)
      {
        jetpt   = jetpttmp; 
        jeteta  = p_jet.jetEta()[ireco];
        jetphi  = p_jet.jetPhi()[ireco];
        jetmult = p_jet.chargedMultiplicity()[ireco];
      }
    }
    hJetEtaPt->Fill(jetpt,jeteta);
    hJetMult->Fill(jetmult);

    TVector3 jet;
    jet.SetPtEtaPhi(jetpt,jeteta,jetphi);

    if(jetpt>500)
    {
      int multBin = getNMultBin(jetmult);

      for (size_t ireco=0; ireco<recosize_ks; ireco++) {
         if( fabs(pdgId_ks[ireco]) != 310 ) continue;
         if(p_ks.cand_decayLength3D()[ireco]/p_ks.cand_decayLengthError3D()[ireco]<5) continue;
         if(cos(p_ks.cand_angle3D()[ireco])<0.99995) continue;

         hKsRecoAll->Fill(p_ks.cand_pT()[ireco], p_ks.cand_y()[ireco]);
         hKsRecoAll_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);

         TVector3 ks;
         ks.SetPtEtaPhi(p_ks.cand_pT()[ireco],p_ks.cand_eta()[ireco],p_ks.cand_phi()[ireco]);
         if(ks.DeltaR(jet)<0.4) 
         {
           hJet4KsRecoAll_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);
           hJet4KsRecoMult_mass[multBin]->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);
         }
         if(ks.DeltaR(jet)<0.8)
         {
           hJet8KsRecoAll_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);
           hJet8KsRecoMult_mass[multBin]->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);
         }
      }

      for (size_t ireco=0; ireco<recosize_lam; ireco++) {
         if( fabs(pdgId_lam[ireco]) != 3122 ) continue;

         if(p_lam.cand_decayLength3D()[ireco]/p_lam.cand_decayLengthError3D()[ireco]<5) continue;
         if(cos(p_lam.cand_angle3D()[ireco])<0.99995) continue;

         auto dauIdx = dauIdxEvt_lam.at(ireco);

         auto trkIdx0 = p_lam.cand_trkIdx().at(dauIdx.at(0));
         auto trkIdx1 = p_lam.cand_trkIdx().at(dauIdx.at(1));
/*
         if(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx0),2))<3) continue;
         if(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx0),2))>45) continue;
         if(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx1),2))<5) continue;
         if(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx1),2))>35) continue;

*/
/*
         if(fabs(p_lam.trk_xyDCASignificance().at(trkIdx0))<5.) continue;
         if(fabs(p_lam.trk_zDCASignificance().at(trkIdx0))<5. || fabs(p_lam.trk_zDCASignificance().at(trkIdx0))>40) continue;
         if(fabs(p_lam.trk_xyDCASignificance().at(trkIdx1))<5.) continue;
         if(fabs(p_lam.trk_zDCASignificance().at(trkIdx1))<5. || fabs(p_lam.trk_zDCASignificance().at(trkIdx1))>30.) continue;
*/
         if(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx0),2))<5) continue;
         if(fabs(p_lam.trk_zDCASignificance().at(trkIdx0))>40) continue;
         if(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx1),2))<5) continue;
         if(fabs(p_lam.trk_zDCASignificance().at(trkIdx1))>30) continue;

/*
         int chgDau0 = p_lam.cand_charge()[dauIdx[0]];
         int chgDau1 = p_lam.cand_charge()[dauIdx[1]];
         double massDau0 = p_lam.cand_mass()[dauIdx[0]];
         double massDau1 = p_lam.cand_mass()[dauIdx[1]];
         auto pdgDau0 = pdgId_lam[dauIdx[0]];
         auto pdgDau1 = pdgId_lam[dauIdx[1]];
*/
         hLamRecoAll->Fill(p_lam.cand_pT()[ireco], p_lam.cand_y()[ireco]);
         hLamRecoAll_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
         
         TVector3 lam;
         lam.SetPtEtaPhi(p_lam.cand_pT()[ireco],p_lam.cand_eta()[ireco],p_lam.cand_phi()[ireco]);
         if(lam.DeltaR(jet)<0.4)
         {
           hJet4LamRecoAll_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
           hJet4LamRecoMult_mass[multBin]->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
         }
         if(lam.DeltaR(jet)<0.8)
         {
           hJet8LamRecoAll_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
           hJet8LamRecoMult_mass[multBin]->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
         }
      }

      for (size_t ireco=0; ireco<recosize_antilam; ireco++) {
         if( fabs(pdgId_antilam[ireco]) != 3122 ) continue;

         if(p_antilam.cand_decayLength3D()[ireco]/p_antilam.cand_decayLengthError3D()[ireco]<5) continue;
         if(cos(p_antilam.cand_angle3D()[ireco])<0.99995) continue;

         auto dauIdx = dauIdxEvt_antilam.at(ireco);

         auto trkIdx0 = p_antilam.cand_trkIdx().at(dauIdx.at(0));
         auto trkIdx1 = p_antilam.cand_trkIdx().at(dauIdx.at(1));
/*
         if(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx0),2))<3) continue;
         if(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx0),2))>45) continue;
         if(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx1),2))<5) continue;
         if(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx1),2))>35) continue;
*/
/*
         if(fabs(p_antilam.trk_xyDCASignificance().at(trkIdx0))<2.) continue;
         if(fabs(p_antilam.trk_zDCASignificance().at(trkIdx0))<2.) continue;
         if(fabs(p_antilam.trk_xyDCASignificance().at(trkIdx1))<2.) continue;
         if(fabs(p_antilam.trk_zDCASignificance().at(trkIdx1))<2.) continue;
*/
         if(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx0),2))<5) continue;
         if(fabs(p_antilam.trk_zDCASignificance().at(trkIdx0))>40) continue;
         if(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx1),2))<5) continue;
         if(fabs(p_antilam.trk_zDCASignificance().at(trkIdx1))>30) continue;

/*
         int chgDau0 = p_antilam.cand_charge()[dauIdx[0]];
         int chgDau1 = p_antilam.cand_charge()[dauIdx[1]];
         double massDau0 = p_antilam.cand_mass()[dauIdx[0]];
         double massDau1 = p_antilam.cand_mass()[dauIdx[1]];
         auto pdgDau0 = pdgId_antilam[dauIdx[0]];
         auto pdgDau1 = pdgId_antilam[dauIdx[1]];
*/
         hLamRecoAll->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_y()[ireco]);
         hLamRecoAll_mass->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);

         TVector3 antilam;
         antilam.SetPtEtaPhi(p_antilam.cand_pT()[ireco],p_antilam.cand_eta()[ireco],p_antilam.cand_phi()[ireco]);
         if(antilam.DeltaR(jet)<0.4)
         {
           hJet4LamRecoAll_mass->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);
           hJet4LamRecoMult_mass[multBin]->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);
         }
         if(antilam.DeltaR(jet)<0.8)
         {
           hJet8LamRecoAll_mass->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);
           hJet8LamRecoMult_mass[multBin]->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);
         }
      }
     
      for (size_t ireco=0; ireco<recosize_xi; ireco++) {
         if( fabs(pdgId_xi[ireco]) != 3312 ) continue;
         if(p_xi.cand_decayLength3D()[ireco]/p_xi.cand_decayLengthError3D()[ireco]<3) continue;
         if(cos(p_xi.cand_angle3D()[ireco])<0.99995) continue;

         auto dauIdx = dauIdxEvt_xi.at(ireco);
         if(p_xi.cand_decayLength3D()[dauIdx[1]]/p_xi.cand_decayLengthError3D()[dauIdx[1]]<15) continue;
         if(fabs(p_xi.cand_mass()[dauIdx[1]]-1.116)>0.005) continue;

         auto trkIdx0 = p_xi.cand_trkIdx().at(dauIdx.at(0));
         if(sqrt(pow(p_xi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(trkIdx0),2))<4) continue;
//         if(sqrt(pow(p_xi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(trkIdx0),2))>50) continue;
         if(fabs(p_xi.trk_zDCASignificance().at(trkIdx0))>40) continue;

         auto gdauIdx = dauIdxEvt_xi.at(dauIdx[1]);
         auto gtrkIdx0 = p_xi.cand_trkIdx().at(gdauIdx.at(0));
         auto gtrkIdx1 = p_xi.cand_trkIdx().at(gdauIdx.at(1));

         if(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx0),2))<15) continue;
//         if(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx0),2))>60) continue;
         if(fabs(p_xi.trk_zDCASignificance().at(gtrkIdx0))>50) continue;
         if(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx1),2))<10) continue;
//         if(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx1),2))>45) continue;
         if(fabs(p_xi.trk_zDCASignificance().at(gtrkIdx1))>35) continue;

         hXiRecoAll->Fill(p_xi.cand_pT()[ireco], p_xi.cand_y()[ireco]);
         hXiRecoAll_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);

         TVector3 xi;
         xi.SetPtEtaPhi(p_xi.cand_pT()[ireco],p_xi.cand_eta()[ireco],p_xi.cand_phi()[ireco]);
         if(xi.DeltaR(jet)<0.4)
         {
           hJet4XiRecoAll_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
           hJet4XiRecoMult_mass[multBin]->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
         }
         if(xi.DeltaR(jet)<0.8)
         {
           hJet8XiRecoAll_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
           hJet8XiRecoMult_mass[multBin]->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
         }
      }

      for (size_t ireco=0; ireco<recosize_antixi; ireco++) {
         if( fabs(pdgId_antixi[ireco]) != 3312 ) continue;
         if(p_antixi.cand_decayLength3D()[ireco]/p_antixi.cand_decayLengthError3D()[ireco]<3) continue;
         if(cos(p_antixi.cand_angle3D()[ireco])<0.99995) continue;

         auto dauIdx = dauIdxEvt_antixi.at(ireco);
         if(p_antixi.cand_decayLength3D()[dauIdx[1]]/p_antixi.cand_decayLengthError3D()[dauIdx[1]]<15) continue;
         if(fabs(p_antixi.cand_mass()[dauIdx[1]]-1.116)>0.005) continue;

         auto trkIdx0 = p_antixi.cand_trkIdx().at(dauIdx.at(0));
         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(trkIdx0),2))<4) continue;
//         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(trkIdx0),2))>50) continue;
         if(fabs(p_antixi.trk_zDCASignificance().at(trkIdx0))>40) continue;

         auto gdauIdx = dauIdxEvt_antixi.at(dauIdx[1]);
         auto gtrkIdx0 = p_antixi.cand_trkIdx().at(gdauIdx.at(0));
         auto gtrkIdx1 = p_antixi.cand_trkIdx().at(gdauIdx.at(1));

         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx0),2))<15) continue;
//         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx0),2))>60) continue;
         if(fabs(p_antixi.trk_zDCASignificance().at(gtrkIdx0))>50) continue;
         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx1),2))<10) continue;
//         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx1),2))>45) continue;
         if(fabs(p_antixi.trk_zDCASignificance().at(gtrkIdx1))>35) continue;

         hXiRecoAll->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_y()[ireco]);
         hXiRecoAll_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);

         TVector3 antixi;
         antixi.SetPtEtaPhi(p_antixi.cand_pT()[ireco],p_antixi.cand_eta()[ireco],p_antixi.cand_phi()[ireco]);
         if(antixi.DeltaR(jet)<0.4)
         {
           hJet4XiRecoAll_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);
           hJet4XiRecoMult_mass[multBin]->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);
         }
         if(antixi.DeltaR(jet)<0.8)
         {
           hJet8XiRecoAll_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);
           hJet8XiRecoMult_mass[multBin]->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);
         }
      }

      for (size_t ireco=0; ireco<recosize_om; ireco++) {
         if( fabs(pdgId_om[ireco]) != 3334 ) continue;
         if(p_om.cand_decayLength3D()[ireco]/p_om.cand_decayLengthError3D()[ireco]<2) continue;
         if(cos(p_om.cand_angle3D()[ireco])<0.99995) continue;

         auto dauIdx = dauIdxEvt_om.at(ireco);
         if(p_om.cand_decayLength3D()[dauIdx[1]]/p_om.cand_decayLengthError3D()[dauIdx[1]]<10) continue;
         if(fabs(p_om.cand_mass()[dauIdx[1]]-1.116)>0.005) continue;

         auto trkIdx0 = p_om.cand_trkIdx().at(dauIdx.at(0));
         if(sqrt(pow(p_om.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_om.trk_zDCASignificance().at(trkIdx0),2))<4) continue;
         if(sqrt(pow(p_om.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_om.trk_zDCASignificance().at(trkIdx0),2))>60) continue;

         auto gdauIdx = dauIdxEvt_om.at(dauIdx[1]);
         auto gtrkIdx0 = p_om.cand_trkIdx().at(gdauIdx.at(0));
         auto gtrkIdx1 = p_om.cand_trkIdx().at(gdauIdx.at(1));

         if(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx0),2))<8) continue;
         if(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx0),2))>60) continue;
         if(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx1),2))<8) continue;
         if(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx1),2))>60) continue;

         hOmRecoAll->Fill(p_om.cand_pT()[ireco], p_om.cand_y()[ireco]);
         hOmRecoAll_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);

         TVector3 om;
         om.SetPtEtaPhi(p_om.cand_pT()[ireco],p_om.cand_eta()[ireco],p_om.cand_phi()[ireco]);
         if(om.DeltaR(jet)<0.4)
         {
           hJet4OmRecoAll_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
           hJet4OmRecoMult_mass[multBin]->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
         }
         if(om.DeltaR(jet)<0.8)
         {
           hJet8OmRecoAll_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
           hJet8OmRecoMult_mass[multBin]->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
         }
      }

      for (size_t ireco=0; ireco<recosize_antiom; ireco++) {
         if( fabs(pdgId_antiom[ireco]) != 3334 ) continue;
         if(p_antiom.cand_decayLength3D()[ireco]/p_antiom.cand_decayLengthError3D()[ireco]<2) continue;
         if(cos(p_antiom.cand_angle3D()[ireco])<0.99995) continue;

         auto dauIdx = dauIdxEvt_antiom.at(ireco);
         if(p_antiom.cand_decayLength3D()[dauIdx[1]]/p_antiom.cand_decayLengthError3D()[dauIdx[1]]<10) continue;
         if(fabs(p_antiom.cand_mass()[dauIdx[1]]-1.116)>0.005) continue;

         auto trkIdx0 = p_antiom.cand_trkIdx().at(dauIdx.at(0));
         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(trkIdx0),2))<4) continue;
         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(trkIdx0),2))>60) continue;

         auto gdauIdx = dauIdxEvt_antiom.at(dauIdx[1]);
         auto gtrkIdx0 = p_antiom.cand_trkIdx().at(gdauIdx.at(0));
         auto gtrkIdx1 = p_antiom.cand_trkIdx().at(gdauIdx.at(1));

         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx0),2))<8) continue;
         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx0),2))>60) continue;
         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx1),2))<8) continue;
         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx1),2))>60) continue;

         hOmRecoAll->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_y()[ireco]);
         hOmRecoAll_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);

         TVector3 antiom;
         antiom.SetPtEtaPhi(p_antiom.cand_pT()[ireco],p_antiom.cand_eta()[ireco],p_antiom.cand_phi()[ireco]);
         if(antiom.DeltaR(jet)<0.4)
         {
           hJet4OmRecoAll_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);
           hJet4OmRecoMult_mass[multBin]->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);
         }
         if(antiom.DeltaR(jet)<0.8)
         {
           hJet8OmRecoAll_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);
           hJet8OmRecoMult_mass[multBin]->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);
         }
      }
    }
  }

  TFile* fout = new TFile(Form("outputs/jetv0_output_%s.root",inputList.Data()),"recreate");

  hKsRecoAll->Write();
  hKsRecoAll_mass->Write();
  hJet4KsRecoAll_mass->Write();
  hJet8KsRecoAll_mass->Write();

  hLamRecoAll->Write();
  hLamRecoAll_mass->Write();
  hJet4LamRecoAll_mass->Write();
  hJet8LamRecoAll_mass->Write();

  hXiRecoAll->Write();
  hXiRecoAll_mass->Write();
  hJet4XiRecoAll_mass->Write();
  hJet8XiRecoAll_mass->Write();

  hOmRecoAll->Write();
  hOmRecoAll_mass->Write();
  hJet4OmRecoAll_mass->Write();
  hJet8OmRecoAll_mass->Write();

  for(int i=0;i<7;i++)
  {
    hJet4KsRecoMult_mass[i]->Write();
    hJet8KsRecoMult_mass[i]->Write();
    hJet4LamRecoMult_mass[i]->Write();
    hJet8LamRecoMult_mass[i]->Write();
    hJet4XiRecoMult_mass[i]->Write();
    hJet8XiRecoMult_mass[i]->Write();
    hJet4OmRecoMult_mass[i]->Write();
    hJet8OmRecoMult_mass[i]->Write();
  }
  hJetEtaPt->Write();
  hJetMult->Write();
  
  fout->Close();
}
