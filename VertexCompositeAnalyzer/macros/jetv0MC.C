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

void jetv0MC(const TString& inputList = "filelist_v0.txt")
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

  std::unique_ptr<TH2D> hKsRecoMatch(new TH2D("hKsRecoMatch", "All channels, Ks p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hKsRecoMatch_mass(new TH2D("hKsRecoMatch_mass", "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55));

  std::unique_ptr<TH2D> hJet4KsRecoMatch_mass(new TH2D("hJet4KsRecoMatch_mass", "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55));

  std::unique_ptr<TH2D> hJet8KsRecoMatch_mass(new TH2D("hJet8KsRecoMatch_mass", "All channels, Ks p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 0.44, 0.55));

  std::unique_ptr<TH2D> hLamRecoMatch(new TH2D("hLamRecoMatch", "All channels, Lam p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hLamRecoMatch_mass(new TH2D("hLamRecoMatch_mass", "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16));

  std::unique_ptr<TH2D> hJet4LamRecoMatch_mass(new TH2D("hJet4LamRecoMatch_mass", "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16));

  std::unique_ptr<TH2D> hJet8LamRecoMatch_mass(new TH2D("hJet8LamRecoMatch_mass", "All channels, Lam p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.08, 1.16));

  std::unique_ptr<TH2D> hXiRecoMatch(new TH2D("hXiRecoMatch", "All channels, Xi p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hXiRecoMatch_mass(new TH2D("hXiRecoMatch_mass", "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39));

  std::unique_ptr<TH2D> hJet4XiRecoMatch_mass(new TH2D("hJet4XiRecoMatch_mass", "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39));

  std::unique_ptr<TH2D> hJet8XiRecoMatch_mass(new TH2D("hJet8XiRecoMatch_mass", "All channels, Xi p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.26, 1.39));

  std::unique_ptr<TH2D> hOmRecoMatch(new TH2D("hOmRecoMatch", "All channels, Om p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hOmRecoMatch_mass(new TH2D("hOmRecoMatch_mass", "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73));

  std::unique_ptr<TH2D> hJet4OmRecoMatch_mass(new TH2D("hJet4OmRecoMatch_mass", "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73));

  std::unique_ptr<TH2D> hJet8OmRecoMatch_mass(new TH2D("hJet8OmRecoMatch_mass", "All channels, Om p3;pT (GeV);mass (GeV);",
          400, 0, 80, 80, 1.61, 1.73));


  std::unique_ptr<TH2D> hKsRecoAll_trkdca1(new TH2D("hKsRecoAll_trkdca1", "All channels, ks p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoAll_trkdca2(new TH2D("hKsRecoAll_trkdca2", "All channels, ks p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoAll_trkdcaxy1(new TH2D("hKsRecoAll_trkdcaxy1", "All channels, ks p3;dca xy;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoAll_trkdcaxy2(new TH2D("hKsRecoAll_trkdcaxy2", "All channels, ks p3;dca xy;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoAll_trkdcaz1(new TH2D("hKsRecoAll_trkdcaz1", "All channels, ks p3;dca z;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoAll_trkdcaz2(new TH2D("hKsRecoAll_trkdcaz2", "All channels, ks p3;dca z;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoAll_angle3D(new TH2D("hKsRecoAll_angle3D", "All channels, ks p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoAll_dl3D(new TH2D("hKsRecoAll_dl3D", "All channels, ks p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoMatch_trkdca1(new TH2D("hKsRecoMatch_trkdca1", "All channels, ks p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoMatch_trkdca2(new TH2D("hKsRecoMatch_trkdca2", "All channels, ks p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoMatch_trkdcaxy1(new TH2D("hKsRecoMatch_trkdcaxy1", "All channels, ks p3;dca xy;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoMatch_trkdcaxy2(new TH2D("hKsRecoMatch_trkdcaxy2", "All channels, ks p3;dca xy;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoMatch_trkdcaz1(new TH2D("hKsRecoMatch_trkdcaz1", "All channels, ks p3;dca z;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoMatch_trkdcaz2(new TH2D("hKsRecoMatch_trkdcaz2", "All channels, ks p3;dca z;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoMatch_angle3D(new TH2D("hKsRecoMatch_angle3D", "All channels, ks p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hKsRecoMatch_dl3D(new TH2D("hKsRecoMatch_dl3D", "All channels, ks p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));

  std::unique_ptr<TH2D> hLamRecoAll_trkdca1(new TH2D("hLamRecoAll_trkdca1", "All channels, lam p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoAll_trkdca2(new TH2D("hLamRecoAll_trkdca2", "All channels, lam p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoAll_trkdcaxy1(new TH2D("hLamRecoAll_trkdcaxy1", "All channels, lam p3;dca xy;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoAll_trkdcaxy2(new TH2D("hLamRecoAll_trkdcaxy2", "All channels, lam p3;dca xy;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoAll_trkdcaz1(new TH2D("hLamRecoAll_trkdcaz1", "All channels, lam p3;dca z;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoAll_trkdcaz2(new TH2D("hLamRecoAll_trkdcaz2", "All channels, lam p3;dca z;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoAll_angle3D(new TH2D("hLamRecoAll_angle3D", "All channels, lam p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoAll_dl3D(new TH2D("hLamRecoAll_dl3D", "All channels, lam p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoMatch_trkdca1(new TH2D("hLamRecoMatch_trkdca1", "All channels, lam p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoMatch_trkdca2(new TH2D("hLamRecoMatch_trkdca2", "All channels, lam p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoMatch_trkdcaxy1(new TH2D("hLamRecoMatch_trkdcaxy1", "All channels, lam p3;dca xy;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoMatch_trkdcaxy2(new TH2D("hLamRecoMatch_trkdcaxy2", "All channels, lam p3;dca xy;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoMatch_trkdcaz1(new TH2D("hLamRecoMatch_trkdcaz1", "All channels, lam p3;dca z;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoMatch_trkdcaz2(new TH2D("hLamRecoMatch_trkdcaz2", "All channels, lam p3;dca z;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoMatch_angle3D(new TH2D("hLamRecoMatch_angle3D", "All channels, lam p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hLamRecoMatch_dl3D(new TH2D("hLamRecoMatch_dl3D", "All channels, lam p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));

  std::unique_ptr<TH2D> hXiRecoAll_trkdca1(new TH2D("hXiRecoAll_trkdca1", "All channels, xi p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoAll_gtrkdca1(new TH2D("hXiRecoAll_gtrkdca1", "All channels, xi p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoAll_gtrkdca2(new TH2D("hXiRecoAll_gtrkdca2", "All channels, xi p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoAll_lamangle3D(new TH2D("hXiRecoAll_lamangle3D", "All channels, xi p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoAll_lamdl3D(new TH2D("hXiRecoAll_lamdl3D", "All channels, xi p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoAll_angle3D(new TH2D("hXiRecoAll_angle3D", "All channels, xi p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoAll_dl3D(new TH2D("hXiRecoAll_dl3D", "All channels, xi p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoMatch_trkdca1(new TH2D("hXiRecoMatch_trkdca1", "All channels, xi p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoMatch_gtrkdca1(new TH2D("hXiRecoMatch_gtrkdca1", "All channels, xi p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoMatch_gtrkdca2(new TH2D("hXiRecoMatch_gtrkdca2", "All channels, xi p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoMatch_lamangle3D(new TH2D("hXiRecoMatch_lamangle3D", "All channels, xi p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoMatch_lamdl3D(new TH2D("hXiRecoMatch_lamdl3D", "All channels, xi p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoMatch_angle3D(new TH2D("hXiRecoMatch_angle3D", "All channels, xi p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hXiRecoMatch_dl3D(new TH2D("hXiRecoMatch_dl3D", "All channels, xi p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));

  std::unique_ptr<TH2D> hOmRecoAll_trkdca1(new TH2D("hOmRecoAll_trkdca1", "All channels, om p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoAll_gtrkdca1(new TH2D("hOmRecoAll_gtrkdca1", "All channels, om p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoAll_gtrkdca2(new TH2D("hOmRecoAll_gtrkdca2", "All channels, om p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoAll_lamangle3D(new TH2D("hOmRecoAll_lamangle3D", "All channels, om p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoAll_lamdl3D(new TH2D("hOmRecoAll_lamdl3D", "All channels, om p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoAll_angle3D(new TH2D("hOmRecoAll_angle3D", "All channels, om p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoAll_dl3D(new TH2D("hOmRecoAll_dl3D", "All channels, om p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoMatch_trkdca1(new TH2D("hOmRecoMatch_trkdca1", "All channels, om p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoMatch_gtrkdca1(new TH2D("hOmRecoMatch_gtrkdca1", "All channels, om p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoMatch_gtrkdca2(new TH2D("hOmRecoMatch_gtrkdca2", "All channels, om p3;dca;p_{T}",
          1000, 0, 100, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoMatch_lamangle3D(new TH2D("hOmRecoMatch_lamangle3D", "All channels, om p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoMatch_lamdl3D(new TH2D("hOmRecoMatch_lamdl3D", "All channels, om p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoMatch_angle3D(new TH2D("hOmRecoMatch_angle3D", "All channels, om p3;angle3D;p_{T}",
          1000, 0.99, 1, 10, 0, 6));
  std::unique_ptr<TH2D> hOmRecoMatch_dl3D(new TH2D("hOmRecoMatch_dl3D", "All channels, om p3;dl3D;p_{T}",
          1000, 0, 200, 10, 0, 6));

  std::unique_ptr<TH2D> hKsGen(new TH2D("hKsGen", "All channels, Ks p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hLamGen(new TH2D("hLamGen", "All channels, Lam p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hXiGen(new TH2D("hXiGen", "All channels, Xi p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hOmGen(new TH2D("hOmGen", "All channels, Om p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet4KsGen(new TH2D("hJet4KsGen", "All channels, Ks p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet8KsGen(new TH2D("hJet8KsGen", "All channels, Ks p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet4LamGen(new TH2D("hJet4LamGen", "All channels, Lam p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet8LamGen(new TH2D("hJet8LamGen", "All channels, Lam p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet4XiGen(new TH2D("hJet4XiGen", "All channels, Xi p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet8XiGen(new TH2D("hJet8XiGen", "All channels, Xi p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet4OmGen(new TH2D("hJet4OmGen", "All channels, Om p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet8OmGen(new TH2D("hJet8OmGen", "All channels, Om p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet4KsGenMatch(new TH2D("hJet4KsGenMatch", "All channels, Ks p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet8KsGenMatch(new TH2D("hJet8KsGenMatch", "All channels, Ks p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet4LamGenMatch(new TH2D("hJet4LamGenMatch", "All channels, Lam p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet8LamGenMatch(new TH2D("hJet8LamGenMatch", "All channels, Lam p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet4XiGenMatch(new TH2D("hJet4XiGenMatch", "All channels, Xi p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet8XiGenMatch(new TH2D("hJet8XiGenMatch", "All channels, Xi p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet4OmGenMatch(new TH2D("hJet4OmGenMatch", "All channels, Om p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));
  std::unique_ptr<TH2D> hJet8OmGenMatch(new TH2D("hJet8OmGenMatch", "All channels, Om p3;pT (GeV);y;",
          400, 0, 80, 10, -2.5, 2.5));

  std::unique_ptr<TH2D> hJetEtaPt(new TH2D("hJetEtaPt", ";pT (GeV);#eta;",
          100, 0, 1000, 10, -2.5, 2.5));

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
  ParticleTreeMC p_ks(&t_ks);

  TChain t_antilam(treeDir_antilam+"/ParticleTree");
  t_antilam.AddFileInfoList(tf.GetList());
  ParticleTreeMC p_antilam(&t_antilam);

  TChain t_antixi(treeDir_antixi+"/ParticleTree");
  t_antixi.AddFileInfoList(tf.GetList());
  ParticleTreeMC p_antixi(&t_antixi);

  TChain t_antiom(treeDir_antiom+"/ParticleTree");
  t_antiom.AddFileInfoList(tf.GetList());
  ParticleTreeMC p_antiom(&t_antiom);

  TChain t_lam(treeDir_lam+"/ParticleTree");
  t_lam.AddFileInfoList(tf.GetList());
//  t_lam.Add(&t_antilam);
  ParticleTreeMC p_lam(&t_lam);

  TChain t_xi(treeDir_xi+"/ParticleTree");
  t_xi.AddFileInfoList(tf.GetList());
//  t_xi.Add(&t_antixi);
  ParticleTreeMC p_xi(&t_xi);

  TChain t_om(treeDir_om+"/ParticleTree");
  t_om.AddFileInfoList(tf.GetList());
//  t_om.Add(&t_antiom);
  ParticleTreeMC p_om(&t_om);

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

    auto gensize_ks = p_ks.gen_mass().size();
    auto gensize_lam = p_lam.gen_mass().size();
    auto gensize_xi = p_xi.gen_mass().size();
    auto gensize_om = p_om.gen_mass().size();
    auto gensize_antilam = p_antilam.gen_mass().size();
    auto gensize_antixi = p_antixi.gen_mass().size();
    auto gensize_antiom = p_antiom.gen_mass().size();

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

    TVector3 jet;
    jet.SetPtEtaPhi(jetpt,jeteta,jetphi);

    if(jetpt>500)
    {
      for (size_t ireco=0; ireco<recosize_ks; ireco++) {
         if( fabs(pdgId_ks[ireco]) != 310 ) continue;
         hKsRecoAll_angle3D->Fill(cos(p_ks.cand_angle3D()[ireco]), p_ks.cand_pT()[ireco]);
         hKsRecoAll_dl3D->Fill(fabs(p_ks.cand_decayLength3D()[ireco]/p_ks.cand_decayLengthError3D()[ireco]),p_ks.cand_pT()[ireco]);

//         if(p_ks.cand_decayLength3D()[ireco]/p_ks.cand_decayLengthError3D()[ireco]<5) continue;
//         if(cos(p_ks.cand_angle3D()[ireco])<0.9999) continue;

         auto dauIdx = dauIdxEvt_ks.at(ireco);

         auto trkIdx0 = p_ks.cand_trkIdx().at(dauIdx.at(0));
         auto trkIdx1 = p_ks.cand_trkIdx().at(dauIdx.at(1));

         hKsRecoAll_trkdca1->Fill(sqrt(pow(p_ks.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_ks.trk_zDCASignificance().at(trkIdx0),2)),p_ks.cand_pT()[ireco]);
         hKsRecoAll_trkdca2->Fill(sqrt(pow(p_ks.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_ks.trk_zDCASignificance().at(trkIdx1),2)),p_ks.cand_pT()[ireco]);
         hKsRecoAll_trkdcaxy1->Fill(p_ks.trk_xyDCASignificance().at(trkIdx0), p_ks.cand_pT()[ireco]);
         hKsRecoAll_trkdcaz1->Fill(p_ks.trk_zDCASignificance().at(trkIdx0), p_ks.cand_pT()[ireco]);
         hKsRecoAll_trkdcaxy2->Fill(p_ks.trk_xyDCASignificance().at(trkIdx1), p_ks.cand_pT()[ireco]);
         hKsRecoAll_trkdcaz2->Fill(p_ks.trk_zDCASignificance().at(trkIdx1), p_ks.cand_pT()[ireco]);

         hKsRecoAll->Fill(p_ks.cand_pT()[ireco], p_ks.cand_y()[ireco]);
         hKsRecoAll_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);

         TVector3 ks;
         ks.SetPtEtaPhi(p_ks.cand_pT()[ireco],p_ks.cand_eta()[ireco],p_ks.cand_phi()[ireco]);
         if(ks.DeltaR(jet)<0.4) 
           hJet4KsRecoAll_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);
         if(ks.DeltaR(jet)<0.8)
           hJet8KsRecoAll_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);

         for (size_t igen=0; igen<gensize_ks; igen++)
         {
           TVector3 ks_mc;
           ks_mc.SetPtEtaPhi(p_ks.gen_pT()[igen],p_ks.gen_eta()[igen],p_ks.gen_phi()[igen]);

           double deltaR = ks.DeltaR(ks_mc);
           if(deltaR<0.01 && fabs(p_ks.cand_p()[ireco]-p_ks.gen_p()[igen])<0.05*p_ks.gen_p()[igen])
           {
             hKsRecoMatch_angle3D->Fill(cos(p_ks.cand_angle3D()[ireco]), p_ks.cand_pT()[ireco]);
             hKsRecoMatch_dl3D->Fill(fabs(p_ks.cand_decayLength3D()[ireco]/p_ks.cand_decayLengthError3D()[ireco]),p_ks.cand_pT()[ireco]);
             hKsRecoMatch_trkdca1->Fill(sqrt(pow(p_ks.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_ks.trk_zDCASignificance().at(trkIdx0),2)),p_ks.cand_pT()[ireco]);
             hKsRecoMatch_trkdca2->Fill(sqrt(pow(p_ks.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_ks.trk_zDCASignificance().at(trkIdx1),2)),p_ks.cand_pT()[ireco]);
             hKsRecoMatch_trkdcaxy1->Fill(p_ks.trk_xyDCASignificance().at(trkIdx0), p_ks.cand_pT()[ireco]);
             hKsRecoMatch_trkdcaz2->Fill(p_ks.trk_zDCASignificance().at(trkIdx0), p_ks.cand_pT()[ireco]);
             hKsRecoMatch_trkdcaxy2->Fill(p_ks.trk_xyDCASignificance().at(trkIdx1), p_ks.cand_pT()[ireco]);
             hKsRecoMatch_trkdcaz2->Fill(p_ks.trk_zDCASignificance().at(trkIdx1), p_ks.cand_pT()[ireco]);

             hKsRecoMatch->Fill(p_ks.cand_pT()[ireco], p_ks.cand_y()[ireco]);
             hKsRecoMatch_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);

             if(ks.DeltaR(jet)<0.4)
             {
               hJet4KsRecoMatch_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);
               hJet4KsGenMatch->Fill(p_ks.gen_pT()[igen], p_ks.gen_y()[igen]);
             }
             if(ks.DeltaR(jet)<0.8)
             {
               hJet8KsRecoMatch_mass->Fill(p_ks.cand_pT()[ireco], p_ks.cand_mass()[ireco]);
               hJet8KsGenMatch->Fill(p_ks.gen_pT()[igen], p_ks.gen_y()[igen]);
             }
           }
         }
      }

      for (size_t ireco=0; ireco<recosize_lam; ireco++) {
         if( fabs(pdgId_lam[ireco]) != 3122 ) continue;

         hLamRecoAll_angle3D->Fill(cos(p_lam.cand_angle3D()[ireco]), p_lam.cand_pT()[ireco]);
         hLamRecoAll_dl3D->Fill(fabs(p_lam.cand_decayLength3D()[ireco]/p_lam.cand_decayLengthError3D()[ireco]),p_lam.cand_pT()[ireco]);

//         if(p_lam.cand_decayLength3D()[ireco]/p_lam.cand_decayLengthError3D()[ireco]<5) continue;
//         if(cos(p_lam.cand_angle3D()[ireco])<0.99995) continue;

         auto dauIdx = dauIdxEvt_lam.at(ireco);

         auto trkIdx0 = p_lam.cand_trkIdx().at(dauIdx.at(0));
         auto trkIdx1 = p_lam.cand_trkIdx().at(dauIdx.at(1));

         hLamRecoAll_trkdca1->Fill(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx0),2)),p_lam.cand_pT()[ireco]);
         hLamRecoAll_trkdca2->Fill(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx1),2)),p_lam.cand_pT()[ireco]);
         hLamRecoAll_trkdcaxy1->Fill(p_lam.trk_xyDCASignificance().at(trkIdx0), p_lam.cand_pT()[ireco]);
         hLamRecoAll_trkdcaz1->Fill(p_lam.trk_zDCASignificance().at(trkIdx0), p_lam.cand_pT()[ireco]);
         hLamRecoAll_trkdcaxy2->Fill(p_lam.trk_xyDCASignificance().at(trkIdx1), p_lam.cand_pT()[ireco]);
         hLamRecoAll_trkdcaz2->Fill(p_lam.trk_zDCASignificance().at(trkIdx1), p_lam.cand_pT()[ireco]);

/*
         if(fabs(p_lam.trk_xyDCASignificance().at(trkIdx0))<2.) continue;
         if(fabs(p_lam.trk_zDCASignificance().at(trkIdx0))<2.) continue;
         if(fabs(p_lam.trk_xyDCASignificance().at(trkIdx1))<2.) continue;
         if(fabs(p_lam.trk_zDCASignificance().at(trkIdx1))<2.) continue;
*/
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
           hJet4LamRecoAll_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
         if(lam.DeltaR(jet)<0.8)
           hJet8LamRecoAll_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);

         for (size_t igen=0; igen<gensize_lam; igen++)
         {
           TVector3 lam_mc;
           lam_mc.SetPtEtaPhi(p_lam.gen_pT()[igen],p_lam.gen_eta()[igen],p_lam.gen_phi()[igen]);

           double deltaR = lam.DeltaR(lam_mc);
           if(deltaR<0.01 && fabs(p_lam.cand_p()[ireco]-p_lam.gen_p()[igen])<0.05*p_lam.gen_p()[igen])
           {
             hLamRecoMatch_angle3D->Fill(cos(p_lam.cand_angle3D()[ireco]),p_lam.cand_pT()[ireco]);
             hLamRecoMatch_dl3D->Fill(fabs(p_lam.cand_decayLength3D()[ireco]/p_lam.cand_decayLengthError3D()[ireco]),p_lam.cand_pT()[ireco]);

             hLamRecoMatch_trkdca1->Fill(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx0),2)),p_lam.cand_pT()[ireco]);
             hLamRecoMatch_trkdca2->Fill(sqrt(pow(p_lam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_lam.trk_zDCASignificance().at(trkIdx1),2)),p_lam.cand_pT()[ireco]);
             hLamRecoMatch_trkdcaxy1->Fill(p_lam.trk_xyDCASignificance().at(trkIdx0), p_lam.cand_pT()[ireco]);
             hLamRecoMatch_trkdcaz1->Fill(p_lam.trk_zDCASignificance().at(trkIdx0), p_lam.cand_pT()[ireco]);
             hLamRecoMatch_trkdcaxy2->Fill(p_lam.trk_xyDCASignificance().at(trkIdx1), p_lam.cand_pT()[ireco]);
             hLamRecoMatch_trkdcaz2->Fill(p_lam.trk_zDCASignificance().at(trkIdx1), p_lam.cand_pT()[ireco]);

             hLamRecoMatch->Fill(p_lam.cand_pT()[ireco], p_lam.cand_y()[ireco]);
             hLamRecoMatch_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);

             if(lam.DeltaR(jet)<0.4)
             {
               hJet4LamRecoMatch_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
               hJet4LamGenMatch->Fill(p_lam.gen_pT()[igen], p_lam.gen_y()[igen]);
             }
             if(lam.DeltaR(jet)<0.8)
             {
               hJet8LamRecoMatch_mass->Fill(p_lam.cand_pT()[ireco], p_lam.cand_mass()[ireco]);
               hJet8LamGenMatch->Fill(p_lam.gen_pT()[igen], p_lam.gen_y()[igen]);
             }
           }
         }
      }

      for (size_t ireco=0; ireco<recosize_antilam; ireco++) {
         if( fabs(pdgId_antilam[ireco]) != 3122 ) continue;

         hLamRecoAll_angle3D->Fill(cos(p_antilam.cand_angle3D()[ireco]),p_antilam.cand_pT()[ireco]);
         hLamRecoAll_dl3D->Fill(fabs(p_antilam.cand_decayLength3D()[ireco]/p_antilam.cand_decayLengthError3D()[ireco]),p_antilam.cand_pT()[ireco]);

//         if(p_antilam.cand_decayLength3D()[ireco]/p_antilam.cand_decayLengthError3D()[ireco]<5) continue;
//         if(cos(p_antilam.cand_angle3D()[ireco])<0.99995) continue;

         auto dauIdx = dauIdxEvt_antilam.at(ireco);

         auto trkIdx0 = p_antilam.cand_trkIdx().at(dauIdx.at(0));
         auto trkIdx1 = p_antilam.cand_trkIdx().at(dauIdx.at(1));

         hLamRecoAll_trkdca1->Fill(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx0),2)),p_antilam.cand_pT()[ireco]);
         hLamRecoAll_trkdca2->Fill(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx1),2)),p_antilam.cand_pT()[ireco]);
         hLamRecoAll_trkdcaxy1->Fill(p_antilam.trk_xyDCASignificance().at(trkIdx0),p_antilam.cand_pT()[ireco]);
         hLamRecoAll_trkdcaz1->Fill(p_antilam.trk_zDCASignificance().at(trkIdx0),p_antilam.cand_pT()[ireco]);
         hLamRecoAll_trkdcaxy2->Fill(p_antilam.trk_xyDCASignificance().at(trkIdx1),p_antilam.cand_pT()[ireco]);
         hLamRecoAll_trkdcaz2->Fill(p_antilam.trk_zDCASignificance().at(trkIdx1),p_antilam.cand_pT()[ireco]);

/*
         if(fabs(p_antilam.trk_xyDCASignificance().at(trkIdx0))<2.) continue;
         if(fabs(p_antilam.trk_zDCASignificance().at(trkIdx0))<2.) continue;
         if(fabs(p_antilam.trk_xyDCASignificance().at(trkIdx1))<2.) continue;
         if(fabs(p_antilam.trk_zDCASignificance().at(trkIdx1))<2.) continue;
*/
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
           hJet4LamRecoAll_mass->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);
         if(antilam.DeltaR(jet)<0.8)
           hJet8LamRecoAll_mass->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);

         for (size_t igen=0; igen<gensize_antilam; igen++)
         {
           TVector3 antilam_mc;
           antilam_mc.SetPtEtaPhi(p_antilam.gen_pT()[igen],p_antilam.gen_eta()[igen],p_antilam.gen_phi()[igen]);

           double deltaR = antilam.DeltaR(antilam_mc);
           if(deltaR<0.01 && fabs(p_antilam.cand_p()[ireco]-p_antilam.gen_p()[igen])<0.05*p_antilam.gen_p()[igen])
           {
             hLamRecoMatch_angle3D->Fill(cos(p_antilam.cand_angle3D()[ireco]),p_antilam.cand_pT()[ireco]);
             hLamRecoMatch_dl3D->Fill(fabs(p_antilam.cand_decayLength3D()[ireco]/p_antilam.cand_decayLengthError3D()[ireco]),p_antilam.cand_pT()[ireco]);

             hLamRecoMatch_trkdca1->Fill(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx0),2)),p_antilam.cand_pT()[ireco]);
             hLamRecoMatch_trkdca2->Fill(sqrt(pow(p_antilam.trk_xyDCASignificance().at(trkIdx1),2)+pow(p_antilam.trk_zDCASignificance().at(trkIdx1),2)),p_antilam.cand_pT()[ireco]);
             hLamRecoMatch_trkdcaxy1->Fill(p_antilam.trk_xyDCASignificance().at(trkIdx0),p_antilam.cand_pT()[ireco]);
             hLamRecoMatch_trkdcaz1->Fill(p_antilam.trk_zDCASignificance().at(trkIdx0),p_antilam.cand_pT()[ireco]);
             hLamRecoMatch_trkdcaxy2->Fill(p_antilam.trk_xyDCASignificance().at(trkIdx1),p_antilam.cand_pT()[ireco]);
             hLamRecoMatch_trkdcaz2->Fill(p_antilam.trk_zDCASignificance().at(trkIdx1),p_antilam.cand_pT()[ireco]);

             hLamRecoMatch->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_y()[ireco]);
             hLamRecoMatch_mass->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);

             if(antilam.DeltaR(jet)<0.4)
             {
               hJet4LamRecoMatch_mass->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);
               hJet4LamGenMatch->Fill(p_antilam.gen_pT()[igen], p_antilam.gen_y()[igen]);
             }
             if(antilam.DeltaR(jet)<0.8)
             {
               hJet8LamRecoMatch_mass->Fill(p_antilam.cand_pT()[ireco], p_antilam.cand_mass()[ireco]);
               hJet8LamGenMatch->Fill(p_antilam.gen_pT()[igen], p_antilam.gen_y()[igen]);
             }
           }
         }
      }
     
      for (size_t ireco=0; ireco<recosize_xi; ireco++) {
         if( fabs(pdgId_xi[ireco]) != 3312 ) continue;

         auto dauIdx = dauIdxEvt_xi.at(ireco);
         if(fabs(p_xi.cand_mass()[dauIdx[1]]-1.116)>0.005) continue;

         hXiRecoAll_angle3D->Fill(cos(p_xi.cand_angle3D()[ireco]),p_xi.cand_pT()[ireco]);
         hXiRecoAll_dl3D->Fill(fabs(p_xi.cand_decayLength3D()[ireco]/p_xi.cand_decayLengthError3D()[ireco]),p_xi.cand_pT()[ireco]);

//         if(p_xi.cand_decayLength3D()[ireco]/p_xi.cand_decayLengthError3D()[ireco]<3) continue;
//         if(cos(p_xi.cand_angle3D()[ireco])<0.9999) continue;

         hXiRecoAll_lamangle3D->Fill(cos(p_xi.cand_angle3D()[dauIdx[1]]),p_xi.cand_pT()[ireco]);
         hXiRecoAll_lamdl3D->Fill(fabs(p_xi.cand_decayLength3D()[dauIdx[1]]/p_xi.cand_decayLengthError3D()[dauIdx[1]]),p_xi.cand_pT()[ireco]);

//         if(p_xi.cand_decayLength3D()[dauIdx[1]]/p_xi.cand_decayLengthError3D()[dauIdx[1]]<12) continue;

         auto trkIdx0 = p_xi.cand_trkIdx().at(dauIdx.at(0));
         hXiRecoAll_trkdca1->Fill(sqrt(pow(p_xi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(trkIdx0),2)),p_xi.cand_pT()[ireco]);
//         if(sqrt(pow(p_xi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(trkIdx0),2))<5) continue;

         auto gdauIdx = dauIdxEvt_xi.at(dauIdx[1]);
         auto gtrkIdx0 = p_xi.cand_trkIdx().at(gdauIdx.at(0));
         auto gtrkIdx1 = p_xi.cand_trkIdx().at(gdauIdx.at(1));

         hXiRecoAll_gtrkdca1->Fill(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx0),2)),p_xi.cand_pT()[ireco]);
         hXiRecoAll_gtrkdca2->Fill(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx1),2)),p_xi.cand_pT()[ireco]);
//        if(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx0),2))<4) continue;
//         if(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx1),2))<3) continue;

         hXiRecoAll->Fill(p_xi.cand_pT()[ireco], p_xi.cand_y()[ireco]);
         hXiRecoAll_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);

         TVector3 xi;
         xi.SetPtEtaPhi(p_xi.cand_pT()[ireco],p_xi.cand_eta()[ireco],p_xi.cand_phi()[ireco]);
         if(xi.DeltaR(jet)<0.4)
           hJet4XiRecoAll_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
         if(xi.DeltaR(jet)<0.8)
           hJet8XiRecoAll_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);

         for (size_t igen=0; igen<gensize_xi; igen++)
         {
           TVector3 xi_mc;
           xi_mc.SetPtEtaPhi(p_xi.gen_pT()[igen],p_xi.gen_eta()[igen],p_xi.gen_phi()[igen]);

           double deltaR = xi.DeltaR(xi_mc);
           if(deltaR<0.01 && fabs(p_xi.cand_p()[ireco]-p_xi.gen_p()[igen])<0.05*p_xi.gen_p()[igen])
           {
             hXiRecoMatch_angle3D->Fill(cos(p_xi.cand_angle3D()[ireco]),p_xi.cand_pT()[ireco]);
             hXiRecoMatch_dl3D->Fill(fabs(p_xi.cand_decayLength3D()[ireco]/p_xi.cand_decayLengthError3D()[ireco]),p_xi.cand_pT()[ireco]);

             hXiRecoMatch_lamangle3D->Fill(cos(p_xi.cand_angle3D()[dauIdx[1]]),p_xi.cand_pT()[ireco]);
             hXiRecoMatch_lamdl3D->Fill(fabs(p_xi.cand_decayLength3D()[dauIdx[1]]/p_xi.cand_decayLengthError3D()[dauIdx[1]]),p_xi.cand_pT()[ireco]);

             hXiRecoMatch_trkdca1->Fill(sqrt(pow(p_xi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(trkIdx0),2)),p_xi.cand_pT()[ireco]);
             hXiRecoMatch_gtrkdca1->Fill(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx0),2)),p_xi.cand_pT()[ireco]);
             hXiRecoMatch_gtrkdca2->Fill(sqrt(pow(p_xi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_xi.trk_zDCASignificance().at(gtrkIdx1),2)),p_xi.cand_pT()[ireco]);

             hXiRecoMatch->Fill(p_xi.cand_pT()[ireco], p_xi.cand_y()[ireco]);
             hXiRecoMatch_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);

             if(xi.DeltaR(jet)<0.4)
             {
               hJet4XiRecoMatch_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
               hJet4XiGenMatch->Fill(p_xi.gen_pT()[igen], p_xi.gen_y()[igen]);
             }
             if(xi.DeltaR(jet)<0.8)
             {
               hJet8XiRecoMatch_mass->Fill(p_xi.cand_pT()[ireco], p_xi.cand_mass()[ireco]);
               hJet8XiGenMatch->Fill(p_xi.gen_pT()[igen], p_xi.gen_y()[igen]);
             }
           }
         }
      }

      for (size_t ireco=0; ireco<recosize_antixi; ireco++) {
         if( fabs(pdgId_antixi[ireco]) != 3312 ) continue;
         
         auto dauIdx = dauIdxEvt_antixi.at(ireco);
         if(fabs(p_antixi.cand_mass()[dauIdx[1]]-1.116)>0.005) continue;

         hXiRecoAll_angle3D->Fill(cos(p_antixi.cand_angle3D()[ireco]),p_antixi.cand_pT()[ireco]);
         hXiRecoAll_dl3D->Fill(fabs(p_antixi.cand_decayLength3D()[ireco]/p_antixi.cand_decayLengthError3D()[ireco]),p_antixi.cand_pT()[ireco]);

//         if(p_antixi.cand_decayLength3D()[ireco]/p_antixi.cand_decayLengthError3D()[ireco]<3) continue;

         hXiRecoAll_lamangle3D->Fill(cos(p_antixi.cand_angle3D()[dauIdx[1]]),p_antixi.cand_pT()[ireco]);
         hXiRecoAll_lamdl3D->Fill(fabs(p_antixi.cand_decayLength3D()[dauIdx[1]]/p_antixi.cand_decayLengthError3D()[dauIdx[1]]),p_antixi.cand_pT()[ireco]);

//         if(p_antixi.cand_decayLength3D()[dauIdx[1]]/p_antixi.cand_decayLengthError3D()[dauIdx[1]]<12) continue;

         auto trkIdx0 = p_antixi.cand_trkIdx().at(dauIdx.at(0));
         hXiRecoAll_trkdca1->Fill(sqrt(pow(p_antixi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(trkIdx0),2)),p_antixi.cand_pT()[ireco]);
//         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(trkIdx0),2))<5) continue;

         auto gdauIdx = dauIdxEvt_antixi.at(dauIdx[1]);
         auto gtrkIdx0 = p_antixi.cand_trkIdx().at(gdauIdx.at(0));
         auto gtrkIdx1 = p_antixi.cand_trkIdx().at(gdauIdx.at(1));

         hXiRecoAll_gtrkdca1->Fill(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx0),2)),p_antixi.cand_pT()[ireco]);
         hXiRecoAll_gtrkdca2->Fill(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx1),2)),p_antixi.cand_pT()[ireco]);
//         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx0),2))<4) continue;
//         if(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx1),2))<3) continue;

         hXiRecoAll->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_y()[ireco]);
         hXiRecoAll_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);

         TVector3 antixi;
         antixi.SetPtEtaPhi(p_antixi.cand_pT()[ireco],p_antixi.cand_eta()[ireco],p_antixi.cand_phi()[ireco]);
         if(antixi.DeltaR(jet)<0.4)
           hJet4XiRecoAll_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);
         if(antixi.DeltaR(jet)<0.8)
           hJet8XiRecoAll_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);

         for (size_t igen=0; igen<gensize_antixi; igen++)
         {
           TVector3 antixi_mc;
           antixi_mc.SetPtEtaPhi(p_antixi.gen_pT()[igen],p_antixi.gen_eta()[igen],p_antixi.gen_phi()[igen]);

           double deltaR = antixi.DeltaR(antixi_mc);
           if(deltaR<0.01 && fabs(p_antixi.cand_p()[ireco]-p_antixi.gen_p()[igen])<0.05*p_antixi.gen_p()[igen])
           {
             hXiRecoMatch_angle3D->Fill(cos(p_antixi.cand_angle3D()[ireco]),p_antixi.cand_pT()[ireco]);
             hXiRecoMatch_dl3D->Fill(fabs(p_antixi.cand_decayLength3D()[ireco]/p_antixi.cand_decayLengthError3D()[ireco]),p_antixi.cand_pT()[ireco]);

             hXiRecoMatch_lamangle3D->Fill(cos(p_antixi.cand_angle3D()[dauIdx[1]]),p_antixi.cand_pT()[ireco]);
             hXiRecoMatch_lamdl3D->Fill(fabs(p_antixi.cand_decayLength3D()[dauIdx[1]]/p_antixi.cand_decayLengthError3D()[dauIdx[1]]),p_antixi.cand_pT()[ireco]);

             hXiRecoMatch_trkdca1->Fill(sqrt(pow(p_antixi.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(trkIdx0),2)),p_antixi.cand_pT()[ireco]);
             hXiRecoMatch_gtrkdca1->Fill(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx0),2)),p_antixi.cand_pT()[ireco]);
             hXiRecoMatch_gtrkdca2->Fill(sqrt(pow(p_antixi.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antixi.trk_zDCASignificance().at(gtrkIdx1),2)),p_antixi.cand_pT()[ireco]);

             hXiRecoMatch->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_y()[ireco]);
             hXiRecoMatch_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);

             if(antixi.DeltaR(jet)<0.4)
             {
               hJet4XiRecoMatch_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);
               hJet4XiGenMatch->Fill(p_antixi.gen_pT()[igen], p_antixi.gen_y()[igen]);
             } 
             if(antixi.DeltaR(jet)<0.8)
             {
               hJet8XiRecoMatch_mass->Fill(p_antixi.cand_pT()[ireco], p_antixi.cand_mass()[ireco]);
               hJet8XiGenMatch->Fill(p_antixi.gen_pT()[igen], p_antixi.gen_y()[igen]);
             } 
           }
         }
      }

      for (size_t ireco=0; ireco<recosize_om; ireco++) {
         if( fabs(pdgId_om[ireco]) != 3334 ) continue;

         auto dauIdx = dauIdxEvt_om.at(ireco);
         if(fabs(p_om.cand_mass()[dauIdx[1]]-1.116)>0.005) continue;

         hOmRecoAll_angle3D->Fill(cos(p_om.cand_angle3D()[ireco]),p_om.cand_pT()[ireco]);
         hOmRecoAll_dl3D->Fill(fabs(p_om.cand_decayLength3D()[ireco]/p_om.cand_decayLengthError3D()[ireco]),p_om.cand_pT()[ireco]);

//         if(p_om.cand_decayLength3D()[ireco]/p_om.cand_decayLengthError3D()[ireco]<2) continue;
//         if(cos(p_om.cand_angle3D()[ireco])<0.9999) continue;

         hOmRecoAll_lamangle3D->Fill(cos(p_om.cand_angle3D()[dauIdx[1]]),p_om.cand_pT()[ireco]);
         hOmRecoAll_lamdl3D->Fill(fabs(p_om.cand_decayLength3D()[dauIdx[1]]/p_om.cand_decayLengthError3D()[dauIdx[1]]),p_om.cand_pT()[ireco]);

//         if(p_om.cand_decayLength3D()[dauIdx[1]]/p_om.cand_decayLengthError3D()[dauIdx[1]]<10) continue;

         auto trkIdx0 = p_om.cand_trkIdx().at(dauIdx.at(0));
         hOmRecoAll_trkdca1->Fill(sqrt(pow(p_om.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_om.trk_zDCASignificance().at(trkIdx0),2)),p_om.cand_pT()[ireco]);
//         if(sqrt(pow(p_om.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_om.trk_zDCASignificance().at(trkIdx0),2))<4) continue;

         auto gdauIdx = dauIdxEvt_om.at(dauIdx[1]);
         auto gtrkIdx0 = p_om.cand_trkIdx().at(gdauIdx.at(0));
         auto gtrkIdx1 = p_om.cand_trkIdx().at(gdauIdx.at(1));
         hOmRecoAll_gtrkdca1->Fill(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx0),2)),p_om.cand_pT()[ireco]);
         hOmRecoAll_gtrkdca2->Fill(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx1),2)),p_om.cand_pT()[ireco]);
//         if(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx0),2))<3) continue;
//         if(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx1),2))<2) continue;

         hOmRecoAll->Fill(p_om.cand_pT()[ireco], p_om.cand_y()[ireco]);
         hOmRecoAll_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);

         TVector3 om;
         om.SetPtEtaPhi(p_om.cand_pT()[ireco],p_om.cand_eta()[ireco],p_om.cand_phi()[ireco]);
         if(om.DeltaR(jet)<0.4)
           hJet4OmRecoAll_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
         if(om.DeltaR(jet)<0.8)
           hJet8OmRecoAll_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);

         for (size_t igen=0; igen<gensize_om; igen++)
         {
           TVector3 om_mc;
           om_mc.SetPtEtaPhi(p_om.gen_pT()[igen],p_om.gen_eta()[igen],p_om.gen_phi()[igen]);

           double deltaR = om.DeltaR(om_mc);
           if(deltaR<0.01 && fabs(p_om.cand_p()[ireco]-p_om.gen_p()[igen])<0.05*p_om.gen_p()[igen])
           {
             hOmRecoMatch_angle3D->Fill(cos(p_om.cand_angle3D()[ireco]),p_om.cand_pT()[ireco]);
             hOmRecoMatch_dl3D->Fill(fabs(p_om.cand_decayLength3D()[ireco]/p_om.cand_decayLengthError3D()[ireco]),p_om.cand_pT()[ireco]);

             hOmRecoMatch_lamangle3D->Fill(cos(p_om.cand_angle3D()[dauIdx[1]]),p_om.cand_pT()[ireco]);
             hOmRecoMatch_lamdl3D->Fill(fabs(p_om.cand_decayLength3D()[dauIdx[1]]/p_om.cand_decayLengthError3D()[dauIdx[1]]),p_om.cand_pT()[ireco]);

             hOmRecoMatch_trkdca1->Fill(sqrt(pow(p_om.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_om.trk_zDCASignificance().at(trkIdx0),2)),p_om.cand_pT()[ireco]);
             hOmRecoMatch_gtrkdca1->Fill(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx0),2)),p_om.cand_pT()[ireco]);
             hOmRecoMatch_gtrkdca2->Fill(sqrt(pow(p_om.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_om.trk_zDCASignificance().at(gtrkIdx1),2)),p_om.cand_pT()[ireco]);

             hOmRecoMatch->Fill(p_om.cand_pT()[ireco], p_om.cand_y()[ireco]);
             hOmRecoMatch_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);

             if(om.DeltaR(jet)<0.4)
             {
               hJet4OmRecoMatch_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
               hJet4OmGenMatch->Fill(p_om.gen_pT()[igen], p_om.gen_y()[igen]);
             }
             if(om.DeltaR(jet)<0.8)
             {
               hJet8OmRecoMatch_mass->Fill(p_om.cand_pT()[ireco], p_om.cand_mass()[ireco]);
               hJet8OmGenMatch->Fill(p_om.gen_pT()[igen], p_om.gen_y()[igen]);
             }
           }
         }
      }

      for (size_t ireco=0; ireco<recosize_antiom; ireco++) {
         if( fabs(pdgId_antiom[ireco]) != 3334 ) continue;

         auto dauIdx = dauIdxEvt_antiom.at(ireco);
         if(fabs(p_antiom.cand_mass()[dauIdx[1]]-1.116)>0.005) continue;

         hOmRecoAll_angle3D->Fill(cos(p_antiom.cand_angle3D()[ireco]),p_antiom.cand_pT()[ireco]);
         hOmRecoAll_dl3D->Fill(fabs(p_antiom.cand_decayLength3D()[ireco]/p_antiom.cand_decayLengthError3D()[ireco]),p_antiom.cand_pT()[ireco]);

//         if(p_antiom.cand_decayLength3D()[ireco]/p_antiom.cand_decayLengthError3D()[ireco]<2) continue;
//         if(cos(p_antiom.cand_angle3D()[ireco])<0.9999) continue;

         hOmRecoAll_lamangle3D->Fill(cos(p_antiom.cand_angle3D()[dauIdx[1]]),p_antiom.cand_pT()[ireco]);
         hOmRecoAll_lamdl3D->Fill(fabs(p_antiom.cand_decayLength3D()[dauIdx[1]]/p_antiom.cand_decayLengthError3D()[dauIdx[1]]),p_antiom.cand_pT()[ireco]);

//         if(p_antiom.cand_decayLength3D()[dauIdx[1]]/p_antiom.cand_decayLengthError3D()[dauIdx[1]]<10) continue;

         auto trkIdx0 = p_antiom.cand_trkIdx().at(dauIdx.at(0));
         hOmRecoAll_trkdca1->Fill(sqrt(pow(p_antiom.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(trkIdx0),2)),p_antiom.cand_pT()[ireco]);
//         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(trkIdx0),2))<4) continue;

         auto gdauIdx = dauIdxEvt_antiom.at(dauIdx[1]);
         auto gtrkIdx0 = p_antiom.cand_trkIdx().at(gdauIdx.at(0));
         auto gtrkIdx1 = p_antiom.cand_trkIdx().at(gdauIdx.at(1));
         hOmRecoAll_gtrkdca1->Fill(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx0),2)),p_antiom.cand_pT()[ireco]);
         hOmRecoAll_gtrkdca2->Fill(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx1),2)),p_antiom.cand_pT()[ireco]);
//         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx0),2))<3) continue;
//         if(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx1),2))<2) continue;

         hOmRecoAll->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_y()[ireco]);
         hOmRecoAll_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);

         TVector3 antiom;
         antiom.SetPtEtaPhi(p_antiom.cand_pT()[ireco],p_antiom.cand_eta()[ireco],p_antiom.cand_phi()[ireco]);
         if(antiom.DeltaR(jet)<0.4)
           hJet4OmRecoAll_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);
         if(antiom.DeltaR(jet)<0.8)
           hJet8OmRecoAll_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);

         for (size_t igen=0; igen<gensize_antiom; igen++)
         {
           TVector3 antiom_mc;
           antiom_mc.SetPtEtaPhi(p_antiom.gen_pT()[igen],p_antiom.gen_eta()[igen],p_antiom.gen_phi()[igen]);
         
           double deltaR = antiom.DeltaR(antiom_mc);
           if(deltaR<0.01 && fabs(p_antiom.cand_p()[ireco]-p_antiom.gen_p()[igen])<0.05*p_antiom.gen_p()[igen])
           {
             hOmRecoMatch_angle3D->Fill(cos(p_antiom.cand_angle3D()[ireco]),p_antiom.cand_pT()[ireco]);
             hOmRecoMatch_dl3D->Fill(fabs(p_antiom.cand_decayLength3D()[ireco]/p_antiom.cand_decayLengthError3D()[ireco]),p_antiom.cand_pT()[ireco]);
             
             hOmRecoMatch_lamangle3D->Fill(cos(p_antiom.cand_angle3D()[dauIdx[1]]),p_antiom.cand_pT()[ireco]);
             hOmRecoMatch_lamdl3D->Fill(fabs(p_antiom.cand_decayLength3D()[dauIdx[1]]/p_antiom.cand_decayLengthError3D()[dauIdx[1]]),p_antiom.cand_pT()[ireco]);
             
             hOmRecoMatch_trkdca1->Fill(sqrt(pow(p_antiom.trk_xyDCASignificance().at(trkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(trkIdx0),2)),p_antiom.cand_pT()[ireco]);
             hOmRecoMatch_gtrkdca1->Fill(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx0),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx0),2)),p_antiom.cand_pT()[ireco]);
             hOmRecoMatch_gtrkdca2->Fill(sqrt(pow(p_antiom.trk_xyDCASignificance().at(gtrkIdx1),2)+pow(p_antiom.trk_zDCASignificance().at(gtrkIdx1),2)),p_antiom.cand_pT()[ireco]);

             hOmRecoMatch->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_y()[ireco]);
             hOmRecoMatch_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);
         
             if(antiom.DeltaR(jet)<0.4)
             {
               hJet4OmRecoMatch_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);
               hJet4OmGenMatch->Fill(p_antiom.gen_pT()[igen], p_antiom.gen_y()[igen]);
             }
             if(antiom.DeltaR(jet)<0.8)
             {
               hJet8OmRecoMatch_mass->Fill(p_antiom.cand_pT()[ireco], p_antiom.cand_mass()[ireco]);
               hJet8OmGenMatch->Fill(p_antiom.gen_pT()[igen], p_antiom.gen_y()[igen]);
             } 
           }
         }
      }

      for (size_t igen=0; igen<gensize_ks; igen++)
      { 
        TVector3 ks_mc;
        ks_mc.SetPtEtaPhi(p_ks.gen_pT()[igen],p_ks.gen_eta()[igen],p_ks.gen_phi()[igen]);
       
        hKsGen->Fill(p_ks.gen_pT()[igen], p_ks.gen_y()[igen]);
   
        if(ks_mc.DeltaR(jet)<0.4)
          hJet4KsGen->Fill(p_ks.gen_pT()[igen], p_ks.gen_y()[igen]);
        if(ks_mc.DeltaR(jet)<0.8)
          hJet8KsGen->Fill(p_ks.gen_pT()[igen], p_ks.gen_y()[igen]);
      }

      for (size_t igen=0; igen<gensize_lam; igen++)
      {    
        TVector3 lam_mc;
        lam_mc.SetPtEtaPhi(p_lam.gen_pT()[igen],p_lam.gen_eta()[igen],p_lam.gen_phi()[igen]);
             
        hLamGen->Fill(p_lam.gen_pT()[igen], p_lam.gen_y()[igen]);

        if(lam_mc.DeltaR(jet)<0.4)
          hJet4LamGen->Fill(p_lam.gen_pT()[igen], p_lam.gen_y()[igen]);
        if(lam_mc.DeltaR(jet)<0.8)
          hJet8LamGen->Fill(p_lam.gen_pT()[igen], p_lam.gen_y()[igen]);
      }

      for (size_t igen=0; igen<gensize_xi; igen++)
      {
        TVector3 xi_mc;
        xi_mc.SetPtEtaPhi(p_xi.gen_pT()[igen],p_xi.gen_eta()[igen],p_xi.gen_phi()[igen]);

        hXiGen->Fill(p_xi.gen_pT()[igen], p_xi.gen_y()[igen]);

        if(xi_mc.DeltaR(jet)<0.4)
          hJet4XiGen->Fill(p_xi.gen_pT()[igen], p_xi.gen_y()[igen]);
        if(xi_mc.DeltaR(jet)<0.8)
          hJet8XiGen->Fill(p_xi.gen_pT()[igen], p_xi.gen_y()[igen]);
      }

      for (size_t igen=0; igen<gensize_om; igen++)
      {
        TVector3 om_mc;
        om_mc.SetPtEtaPhi(p_om.gen_pT()[igen],p_om.gen_eta()[igen],p_om.gen_phi()[igen]);

        hOmGen->Fill(p_om.gen_pT()[igen], p_om.gen_y()[igen]);

        if(om_mc.DeltaR(jet)<0.4)
          hJet4OmGen->Fill(p_om.gen_pT()[igen], p_om.gen_y()[igen]);
        if(om_mc.DeltaR(jet)<0.8)
          hJet8OmGen->Fill(p_om.gen_pT()[igen], p_om.gen_y()[igen]);
      }
    }
  }

  TFile* fout = new TFile(Form("outputs/jetv0MC_output_%s.root",inputList.Data()),"recreate");

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

  hKsRecoMatch->Write();
  hKsRecoMatch_mass->Write();
  hJet4KsRecoMatch_mass->Write();
  hJet8KsRecoMatch_mass->Write();

  hLamRecoMatch->Write();
  hLamRecoMatch_mass->Write();
  hJet4LamRecoMatch_mass->Write();
  hJet8LamRecoMatch_mass->Write();

  hXiRecoMatch->Write();
  hXiRecoMatch_mass->Write();
  hJet4XiRecoMatch_mass->Write();
  hJet8XiRecoMatch_mass->Write();

  hOmRecoMatch->Write();
  hOmRecoMatch_mass->Write();
  hJet4OmRecoMatch_mass->Write();
  hJet8OmRecoMatch_mass->Write();

  hKsRecoAll_trkdca1->Write();
  hKsRecoAll_trkdca2->Write();
  hKsRecoAll_trkdcaxy1->Write();
  hKsRecoAll_trkdcaxy2->Write();
  hKsRecoAll_trkdcaz1->Write();
  hKsRecoAll_trkdcaz2->Write();
  hKsRecoAll_angle3D->Write();
  hKsRecoAll_dl3D->Write();
  hKsRecoMatch_trkdca1->Write();
  hKsRecoMatch_trkdca2->Write();
  hKsRecoMatch_trkdcaxy1->Write();
  hKsRecoMatch_trkdcaxy2->Write();
  hKsRecoMatch_trkdcaz1->Write();
  hKsRecoMatch_trkdcaz2->Write();
  hKsRecoMatch_angle3D->Write();
  hKsRecoMatch_dl3D->Write();

  hLamRecoAll_trkdca1->Write();
  hLamRecoAll_trkdca2->Write();
  hLamRecoAll_trkdcaxy1->Write();
  hLamRecoAll_trkdcaxy2->Write();
  hLamRecoAll_trkdcaz1->Write();
  hLamRecoAll_trkdcaz2->Write();
  hLamRecoAll_angle3D->Write();
  hLamRecoAll_dl3D->Write();
  hLamRecoMatch_trkdca1->Write();
  hLamRecoMatch_trkdca2->Write();
  hLamRecoMatch_trkdcaxy1->Write();
  hLamRecoMatch_trkdcaxy2->Write();
  hLamRecoMatch_trkdcaz1->Write();
  hLamRecoMatch_trkdcaz2->Write();
  hLamRecoMatch_angle3D->Write();
  hLamRecoMatch_dl3D->Write();

  hXiRecoAll_trkdca1->Write();
  hXiRecoAll_gtrkdca1->Write();
  hXiRecoAll_gtrkdca2->Write();
  hXiRecoAll_lamangle3D->Write();
  hXiRecoAll_lamdl3D->Write();
  hXiRecoAll_angle3D->Write();
  hXiRecoAll_dl3D->Write();
  hXiRecoMatch_trkdca1->Write();
  hXiRecoMatch_gtrkdca1->Write();
  hXiRecoMatch_gtrkdca2->Write();
  hXiRecoMatch_lamangle3D->Write();
  hXiRecoMatch_lamdl3D->Write();
  hXiRecoMatch_angle3D->Write();
  hXiRecoMatch_dl3D->Write();

  hOmRecoAll_trkdca1->Write();
  hOmRecoAll_gtrkdca1->Write();
  hOmRecoAll_gtrkdca2->Write();
  hOmRecoAll_lamangle3D->Write();
  hOmRecoAll_lamdl3D->Write();
  hOmRecoAll_angle3D->Write();
  hOmRecoAll_dl3D->Write();
  hOmRecoMatch_trkdca1->Write();
  hOmRecoMatch_gtrkdca1->Write();
  hOmRecoMatch_gtrkdca2->Write();
  hOmRecoMatch_lamangle3D->Write();
  hOmRecoMatch_lamdl3D->Write();
  hOmRecoMatch_angle3D->Write();
  hOmRecoMatch_dl3D->Write();

  hKsGen->Write();
  hLamGen->Write();
  hXiGen->Write();
  hOmGen->Write();
  hJet4KsGen->Write();
  hJet8KsGen->Write();
  hJet4LamGen->Write();
  hJet8LamGen->Write();
  hJet4XiGen->Write();
  hJet8XiGen->Write();
  hJet4OmGen->Write();
  hJet8OmGen->Write();
  hJet4KsGenMatch->Write();
  hJet8KsGenMatch->Write();
  hJet4LamGenMatch->Write();
  hJet8LamGenMatch->Write();
  hJet4XiGenMatch->Write();
  hJet8XiGenMatch->Write();
  hJet4OmGenMatch->Write();
  hJet8OmGenMatch->Write();

  hJetEtaPt->Write();
  
  fout->Close();
}
