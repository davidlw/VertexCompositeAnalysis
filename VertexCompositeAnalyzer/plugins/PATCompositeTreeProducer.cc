// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// constants, enums and typedefs
//

#define PI 3.1416
#define MAXCAN 50000
#define MAXDAU 3
#define MAXGDAU 2
#define MAXTRG 1024
#define MAXSEL 100

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;

//
// class decleration
//

class PATCompositeTreeProducer : public edm::EDAnalyzer {
public:
  explicit PATCompositeTreeProducer(const edm::ParameterSet&);
  ~PATCompositeTreeProducer();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&);
  virtual void fillGEN(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void initHistogram();
  virtual void initTree();
  reco::GenParticleRef findLastPar(const reco::GenParticleRef&);
  reco::GenParticleRef findMother(const reco::GenParticleRef&);
  bool hasDaughters(const reco::GenParticleRef&, const std::vector<int>&);

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree* PATCompositeNtuple;
  TH2F*  hMassVsMVA[6][10];
  TH2F*  hpTVsMVA[6][10];
  TH2F*  hetaVsMVA[6][10];
  TH2F*  hyVsMVA[6][10];
  TH2F*  hVtxProbVsMVA[6][10];
  TH2F*  h3DCosPointingAngleVsMVA[6][10];
  TH2F*  h3DPointingAngleVsMVA[6][10];
  TH2F*  h2DCosPointingAngleVsMVA[6][10];
  TH2F*  h2DPointingAngleVsMVA[6][10];
  TH2F*  h3DDecayLengthSignificanceVsMVA[6][10];
  TH2F*  h3DDecayLengthVsMVA[6][10];
  TH2F*  h2DDecayLengthSignificanceVsMVA[6][10];
  TH2F*  h2DDecayLengthVsMVA[6][10];
  TH2F*  h3DDCAVsMVA[6][10];
  TH2F*  h2DDCAVsMVA[6][10];
  TH2F*  hzDCASignificanceDaugtherVsMVA[MAXDAU][6][10];
  TH2F*  hxyDCASignificanceDaugtherVsMVA[MAXDAU][6][10];
  TH2F*  hNHitDVsMVA[MAXDAU][6][10];
  TH2F*  hpTDVsMVA[MAXDAU][6][10];
  TH2F*  hpTerrDVsMVA[MAXDAU][6][10];
  TH2F*  hEtaDVsMVA[MAXDAU][6][10];
  TH2F*  hdedxHarmonic2DVsMVA[MAXDAU][6][10];
  TH2F*  hdedxHarmonic2DVsP[MAXDAU][6][10];

  bool   saveTree_;
  bool   saveHistogram_;
  bool   saveAllHistogram_;
  double massHistPeak_;
  double massHistWidth_;
  int    massHistBins_;

  //options
  bool doRecoNtuple_;
  bool doGenNtuple_;   
  bool doGenMatching_;
  bool doGenMatchingTOF_;
  bool decayInGen_;
  bool twoLayerDecay_;
  bool threeProngDecay_;
  bool doMuon_;
  bool doMuonFull_;
  const std::vector<int> PID_dau_;
  const ushort NDAU_;
  const ushort NGDAU_ = MAXGDAU;

  //cut variables
  double multMax_;
  double multMin_;
  double deltaR_; //deltaR for Gen matching

  std::vector<double> pTBins_;
  std::vector<double> yBins_;

  //tree branches
  //event info
  uint  runNb;
  uint  eventNb;
  uint  lsNb;
  short trigPrescale[MAXTRG];
  short centrality;
  int   Ntrkoffline;
  int   Npixel;
  short nPV;
  uint candSize;
  bool  trigHLT[MAXTRG];
  bool  evtSel[MAXSEL];
  float HFsumETPlus;
  float HFsumETMinus;
  float ZDCPlus;
  float ZDCMinus;
  float bestvx;
  float bestvy;
  float bestvz;
  float ephfpAngle[3];
  float ephfmAngle[3];
  float eptrackmidAngle[3];
  float ephfpQ[3];
  float ephfmQ[3];
  float eptrackmidQ[3];
  float ephfpSumW;
  float ephfmSumW;
  float eptrackmidSumW;

  //Composite candidate info
  float mva[MAXCAN];
  float pt[MAXCAN];
  float eta[MAXCAN];
  float phi[MAXCAN];
  float flavor[MAXCAN];
  float y[MAXCAN];
  float mass[MAXCAN];
  float VtxProb[MAXCAN];
  float dlos[MAXCAN];
  float dl[MAXCAN];
  float dlerror[MAXCAN];
  float agl[MAXCAN];
  float vtxChi2[MAXCAN];
  float ndf[MAXCAN];
  float agl_abs[MAXCAN];
  float agl2D[MAXCAN];
  float agl2D_abs[MAXCAN];
  float dlos2D[MAXCAN];
  float dl2D[MAXCAN];
  bool isSwap[MAXCAN];
  bool matchGEN[MAXCAN];
  int idmom_reco[MAXCAN];
    
  //dau candidate info
  float grand_mass[MAXCAN];
  float grand_VtxProb[MAXCAN];
  float grand_dlos[MAXCAN];
  float grand_dl[MAXCAN];
  float grand_dlerror[MAXCAN];
  float grand_agl[MAXCAN];
  float grand_vtxChi2[MAXCAN];
  float grand_ndf[MAXCAN];
  float grand_agl_abs[MAXCAN];
  float grand_agl2D[MAXCAN];
  float grand_agl2D_abs[MAXCAN];
  float grand_dlos2D[MAXCAN];

  //dau info
  float dzos[MAXDAU][MAXCAN];
  float dxyos[MAXDAU][MAXCAN];
  float nhit[MAXDAU][MAXCAN];
  bool  trkquality[MAXDAU][MAXCAN];
  float ptDau[MAXDAU][MAXCAN];
  float ptErr[MAXDAU][MAXCAN];
  float pDau[MAXDAU][MAXCAN];
  float etaDau[MAXDAU][MAXCAN];
  float phiDau[MAXDAU][MAXCAN];
  short chargeDau[MAXDAU][MAXCAN];
  int   pid[MAXDAU][MAXCAN];
  float tof[MAXDAU][MAXCAN];
  float H2dedx[MAXDAU][MAXCAN];
  float T4dedx[MAXDAU][MAXCAN];
  float trkChi[MAXDAU][MAXCAN];
   
  //grand-dau info
  float grand_dzos[MAXGDAU][MAXCAN];
  float grand_dxyos[MAXGDAU][MAXCAN];
  float grand_nhit[MAXGDAU][MAXCAN];
  bool  grand_trkquality[MAXGDAU][MAXCAN];
  float grand_pt[MAXGDAU][MAXCAN];
  float grand_ptErr[MAXGDAU][MAXCAN];
  float grand_p[MAXGDAU][MAXCAN];
  float grand_eta[MAXGDAU][MAXCAN];
  short grand_charge[MAXGDAU][MAXCAN];
  float grand_H2dedx[MAXGDAU][MAXCAN];
  float grand_T4dedx[MAXGDAU][MAXCAN];
  float grand_trkChi[MAXGDAU][MAXCAN];
    
  //dau muon info
  bool  onestmuon[MAXDAU][MAXCAN];
  bool  pfmuon[MAXDAU][MAXCAN];
  bool  glbmuon[MAXDAU][MAXCAN];
  bool  trkmuon[MAXDAU][MAXCAN];
  bool  tightmuon[MAXDAU][MAXCAN];
  bool  softmuon[MAXDAU][MAXCAN];
  bool  hybridmuon[MAXDAU][MAXCAN];
  bool  hpmuon[MAXDAU][MAXCAN];
  std::vector<std::vector<UChar_t> > trgmuon[MAXDAU];
  short nmatchedst[MAXDAU][MAXCAN];
  short ntrackerlayer[MAXDAU][MAXCAN];
  short npixellayer[MAXDAU][MAXCAN];
  short npixelhit[MAXDAU][MAXCAN];
  short nmuonhit[MAXDAU][MAXCAN];
  float glbtrkchi[MAXDAU][MAXCAN];
  float muonbestdxy[MAXDAU][MAXCAN];
  float muonbestdz[MAXDAU][MAXCAN];
  float muondxy[MAXDAU][MAXCAN];
  float muondz[MAXDAU][MAXCAN];
  short nmatchedch[MAXDAU][MAXCAN];
  float matchedenergy[MAXDAU][MAXCAN];
  float dx_seg[MAXDAU][MAXCAN];
  float dy_seg[MAXDAU][MAXCAN];
  float dxSig_seg[MAXDAU][MAXCAN];
  float dySig_seg[MAXDAU][MAXCAN];
  float ddxdz_seg[MAXDAU][MAXCAN];
  float ddydz_seg[MAXDAU][MAXCAN];
  float ddxdzSig_seg[MAXDAU][MAXCAN];
  float ddydzSig_seg[MAXDAU][MAXCAN];

  // gen info
  std::vector<reco::GenParticleRef> genVec_;
  float weight_gen;
  uint candSize_gen;
  float pt_gen[MAXCAN];
  float eta_gen[MAXCAN];
  float y_gen[MAXCAN];
  short status_gen[MAXCAN];
  int pid_gen[MAXCAN];
  int idmom_gen[MAXCAN];
  short idxrec_gen[MAXCAN];
  int idDau_gen[MAXDAU][MAXCAN];
  short chargeDau_gen[MAXDAU][MAXCAN];
  float ptDau_gen[MAXDAU][MAXCAN];
  float etaDau_gen[MAXDAU][MAXCAN];
  float phiDau_gen[MAXDAU][MAXCAN];

  bool useAnyMVA_;
  bool isSkimMVA_;
  bool isCentrality_;
  bool isEventPlane_;
  bool useDeDxData_;

  //token
  edm::EDGetTokenT<reco::BeamSpot> tok_offlineBS_;
  edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> patCompositeCandidateCollection_Token_;
  edm::EDGetTokenT<MVACollection> MVAValues_Token_;

  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
  edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
  edm::EDGetTokenT<GenEventInfoProduct> tok_genInfo_;

  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

  edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;

  //trigger
  const std::vector<std::string> triggerNames_;
  const std::vector<std::string> filterNames_;
  edm::EDGetTokenT<edm::TriggerResults> tok_triggerResults_;
  const ushort NTRG_;

  //event selection
  const std::vector<std::string> eventFilters_;
  edm::EDGetTokenT<edm::TriggerResults> tok_filterResults_;
  const ushort NSEL_;
  const std::string selectEvents_;

  //prescale provider
  HLTPrescaleProvider hltPrescaleProvider_;
};

//
// static data member definitions
//

//
// constructors and destructor
//

PATCompositeTreeProducer::PATCompositeTreeProducer(const edm::ParameterSet& iConfig) :
  PID_dau_(iConfig.getUntrackedParameter<std::vector<int> >("PID_dau")),
  NDAU_(PID_dau_.size()>MAXDAU ? MAXDAU : PID_dau_.size()),
  patCompositeCandidateCollection_Token_(consumes<pat::CompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"))),
  triggerNames_(iConfig.getUntrackedParameter<std::vector<std::string> >("triggerPathNames")),
  filterNames_(iConfig.getUntrackedParameter<std::vector<std::string> >("triggerFilterNames")),
  NTRG_(triggerNames_.size()>MAXTRG ? MAXTRG : triggerNames_.size()),
  eventFilters_(iConfig.getUntrackedParameter<std::vector<std::string> >("eventFilterNames")),
  NSEL_(eventFilters_.size()>MAXSEL ? MAXSEL : eventFilters_.size()),
  selectEvents_(iConfig.getUntrackedParameter<std::string>("selectEvents")),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
  //options
  doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
  doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
  twoLayerDecay_ = iConfig.getUntrackedParameter<bool>("twoLayerDecay");
  threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
  doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
  doGenMatchingTOF_ = iConfig.getUntrackedParameter<bool>("doGenMatchingTOF");
  decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
  doMuon_ = iConfig.getUntrackedParameter<bool>("doMuon");
  doMuonFull_ = iConfig.getUntrackedParameter<bool>("doMuonFull");

  saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
  saveHistogram_ = iConfig.getUntrackedParameter<bool>("saveHistogram");
  saveAllHistogram_ = iConfig.getUntrackedParameter<bool>("saveAllHistogram");
  massHistPeak_ = iConfig.getUntrackedParameter<double>("massHistPeak");
  massHistWidth_ = iConfig.getUntrackedParameter<double>("massHistWidth");
  massHistBins_ = iConfig.getUntrackedParameter<int>("massHistBins");

  //cut variables
  multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
  multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
  deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);

  pTBins_ = iConfig.getUntrackedParameter<std::vector<double> >("pTBins");
  yBins_  = iConfig.getUntrackedParameter<std::vector<double> >("yBins");

  //input tokens
  tok_offlineBS_ = consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"));
  tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
  MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
  tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));
  tok_genInfo_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

  useDeDxData_ = (iConfig.exists("useDeDxData") ? iConfig.getParameter<bool>("useDeDxData") : false);
  if(useDeDxData_)
  {
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
  }

  isCentrality_ = (iConfig.exists("isCentrality") ? iConfig.getParameter<bool>("isCentrality") : false);
  if(isCentrality_)
  {
    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

  isEventPlane_ = (iConfig.exists("isEventPlane") ? iConfig.getParameter<bool>("isEventPlane") : false);
  if(isEventPlane_) tok_eventplaneSrc_ = consumes<reco::EvtPlaneCollection>(iConfig.getParameter<edm::InputTag>("eventplaneSrc"));

  useAnyMVA_ = (iConfig.exists("useAnyMVA") ? iConfig.getParameter<bool>("useAnyMVA") : false);
  if(useAnyMVA_) MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
  isSkimMVA_ = iConfig.getUntrackedParameter<bool>("isSkimMVA");
  
  tok_triggerResults_ = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultCollection"));
  tok_filterResults_ = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("FilterResultCollection"));
}


PATCompositeTreeProducer::~PATCompositeTreeProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PATCompositeTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //check event
  if(selectEvents_!="")
  {
    edm::Handle<edm::TriggerResults> filterResults;
    iEvent.getByToken(tok_filterResults_, filterResults);
    const auto& filterNames = iEvent.triggerNames(*filterResults);
    const auto& index = filterNames.triggerIndex(selectEvents_);
    if(index<filterNames.size() && filterResults->wasrun(index) && !filterResults->accept(index)) return;
  }
  genVec_.clear();
  if(doRecoNtuple_) fillRECO(iEvent,iSetup);
  if(doGenNtuple_) fillGEN(iEvent,iSetup);
  if(saveTree_) PATCompositeNtuple->Fill();
}


void
PATCompositeTreeProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get collection
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(tok_offlineBS_, beamspot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);
  if(!vertices.isValid()) throw cms::Exception("PATCompositeAnalyzer") << "Primary vertices  collection not found!" << std::endl;

  edm::Handle<pat::CompositeCandidateCollection> v0candidates;
  iEvent.getByToken(patCompositeCandidateCollection_Token_, v0candidates);
  if(!v0candidates.isValid()) throw cms::Exception("PATCompositeAnalyzer") << "V0 candidate collection not found!" << std::endl;

  edm::Handle<MVACollection> mvavalues;
  if(useAnyMVA_)
  {
    iEvent.getByToken(MVAValues_Token_, mvavalues);
    if(!mvavalues.isValid()) throw cms::Exception("PATCompositeAnalyzer") << "MVA collection not found!" << std::endl;
    assert( (*mvavalues).size() == v0candidates->size() );
  }

  edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1 , dEdxHandle2;
  if(useDeDxData_)
  {
    iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
    iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
  }

  runNb = iEvent.id().run();
  eventNb = iEvent.id().event();
  lsNb = iEvent.luminosityBlock();

  //Trigger Information
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(tok_triggerResults_, triggerResults);
  if(triggerNames_.size()>0)
  {
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
    for(ushort iTr=0; iTr<NTRG_; iTr++)
    {
      //Initiliaze the arrays
      trigHLT[iTr] = false;
      trigPrescale[iTr] = -9;
      //Find the trigger index
      const auto& trigName = triggerNames_.at(iTr);
      std::vector<ushort> trgIdxFound;
      for(ushort trgIdx=0; trgIdx<triggerNames.size(); trgIdx++)
      {
        if(triggerNames.triggerName(trgIdx).find(trigName)!=std::string::npos && triggerResults->wasrun(trgIdx)) { trgIdxFound.push_back(trgIdx); }
      }
      short triggerIndex = -1;
      if(trgIdxFound.size()>1)
      {
        for(const auto& trgIdx : trgIdxFound) { if(triggerResults->accept(trgIdx)) { triggerIndex = trgIdx; break; } }
        if(triggerIndex<0) triggerIndex = trgIdxFound[0];
      }
      else if(trgIdxFound.size()==1) triggerIndex = trgIdxFound[0];
      else continue;
      //Check if trigger fired
      bool isTriggerFired = false;
      if(triggerResults->accept(triggerIndex)) isTriggerFired = true;
      //Get the trigger prescale
      int prescaleValue = -1;
      if(hltPrescaleProvider_.hltConfigProvider().inited() && hltPrescaleProvider_.prescaleSet(iEvent,iSetup)>=0)
      {
        const auto& presInfo = hltPrescaleProvider_.prescaleValuesInDetail(iEvent, iSetup, triggerNames.triggerName(triggerIndex));
        const auto& hltPres = presInfo.second;
        const short& l1Pres = ((presInfo.first.size()==1) ? presInfo.first.at(0).second : ((presInfo.first.size()>1) ? 1 : -1));
        prescaleValue = hltPres*l1Pres;
      }
      trigPrescale[iTr] = prescaleValue;
      if(isTriggerFired) trigHLT[iTr] = true;
    }
  }

  //Event selection information
  edm::Handle<edm::TriggerResults> filterResults;
  iEvent.getByToken(tok_filterResults_, filterResults);
  if(eventFilters_.size()>0)
  {
    const edm::TriggerNames& filterNames = iEvent.triggerNames(*filterResults);
    for(ushort iFr=0; iFr<eventFilters_.size(); ++iFr)
    {
      evtSel[iFr] = false;
      const auto& index = filterNames.triggerIndex(eventFilters_.at(iFr));
      if(index < filterNames.size()) evtSel[iFr] = (filterResults->wasrun(index) && filterResults->accept(index));
    }
  }

  centrality = -1;
  if(isCentrality_)
  {
    edm::Handle<reco::Centrality> cent;
    iEvent.getByToken(tok_centSrc_, cent);
    HFsumETPlus = (cent.isValid() ? cent->EtHFtowerSumPlus() : -1.);
    HFsumETMinus = (cent.isValid() ? cent->EtHFtowerSumMinus() : -1.);
    Npixel = (cent.isValid() ? cent->multiplicityPixel() : -1);
    ZDCPlus = (cent.isValid() ? cent->zdcSumPlus() : -1.);
    ZDCMinus = (cent.isValid() ? cent->zdcSumMinus() : -1.);
    Ntrkoffline = (cent.isValid() ? cent->Ntracks() : -1);
      
    edm::Handle<int> cbin;
    iEvent.getByToken(tok_centBinLabel_, cbin);
    centrality = (cbin.isValid() ? *cbin : -1);
  }

  if(isEventPlane_)
  {
    edm::Handle<reco::EvtPlaneCollection> eventplanes;
    iEvent.getByToken(tok_eventplaneSrc_, eventplanes);

    ephfpAngle[0] = (eventplanes.isValid() ? (*eventplanes)[0].angle(2) : -99.);
    ephfpAngle[1] = (eventplanes.isValid() ? (*eventplanes)[6].angle(2) : -99.);
    ephfpAngle[2] = (eventplanes.isValid() ? (*eventplanes)[13].angle(2) : -99.);

    ephfmAngle[0] = (eventplanes.isValid() ? (*eventplanes)[1].angle(2) : -99.);
    ephfmAngle[1] = (eventplanes.isValid() ? (*eventplanes)[7].angle(2) : -99.);
    ephfmAngle[2] = (eventplanes.isValid() ? (*eventplanes)[14].angle(2) : -99.);
    
    eptrackmidAngle[0] = -99.9;
    eptrackmidAngle[1] = (eventplanes.isValid() ? (*eventplanes)[9].angle(2) : -99.);
    eptrackmidAngle[2] = (eventplanes.isValid() ? (*eventplanes)[16].angle(2) : -99.);

    ephfpQ[0] = (eventplanes.isValid() ? (*eventplanes)[0].q(2) : -99.);
    ephfpQ[1] = (eventplanes.isValid() ? (*eventplanes)[6].q(2) : -99.);
    ephfpQ[2] = (eventplanes.isValid() ? (*eventplanes)[13].q(2) : -99.);

    ephfmQ[0] = (eventplanes.isValid() ? (*eventplanes)[1].q(2) : -99.);
    ephfmQ[1] = (eventplanes.isValid() ? (*eventplanes)[7].q(2) : -99.);
    ephfmQ[2] = (eventplanes.isValid() ? (*eventplanes)[14].q(2) : -99.);

    eptrackmidQ[0] = -99.9;
    eptrackmidQ[1] = (eventplanes.isValid() ? (*eventplanes)[9].q(2) : -99.);
    eptrackmidQ[2] = (eventplanes.isValid() ? (*eventplanes)[16].q(2) : -99.);
    
    ephfpSumW = (eventplanes.isValid() ? (*eventplanes)[6].sumw() : -99.);
    ephfmSumW = (eventplanes.isValid() ? (*eventplanes)[7].sumw() : -99.);
    eptrackmidSumW = (eventplanes.isValid() ? (*eventplanes)[9].sumw() : -99.);
  }

  nPV = vertices->size();
  //best vertex
  const auto& vtxPrimary = (vertices->size()>0 ? (*vertices)[0] : reco::Vertex());
  const bool& isPV = (!vtxPrimary.isFake() && vtxPrimary.tracksSize()>=2);
  const auto& bs = (!isPV ? reco::Vertex(beamspot->position(), beamspot->covariance3D()) : reco::Vertex());
  const reco::Vertex& vtx = (isPV ? vtxPrimary : bs);
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  const math::XYZPoint bestvtx(bestvx, bestvy, bestvz);
  const double& bestvzError = vtx.zError(), bestvxError = vtx.xError(), bestvyError = vtx.yError();

  //RECO Candidate info
  candSize = v0candidates->size();
  if(candSize>MAXCAN) throw cms::Exception("PATCompositeAnalyzer") << "Number of candidates (" << candSize << ") exceeds limit!" << std::endl; 
  for(uint it=0; it<candSize; ++it)
  { 
    const auto& trk = (*v0candidates)[it];

    pt[it] = trk.pt();
    eta[it] = trk.eta();
    phi[it] = trk.phi();
    mass[it] = trk.mass();
    y[it] = trk.rapidity();
    flavor[it] = (trk.pdgId()!=0 ? trk.pdgId()/fabs(trk.pdgId()) : 0.);

    mva[it] = (useAnyMVA_ ? (*mvavalues)[it] : 0.0);

    //vtxChi2
    vtxChi2[it] = trk.userFloat("vertexChi2");
    ndf[it] = trk.userFloat("vertexNdof");
    VtxProb[it] = TMath::Prob(vtxChi2[it],ndf[it]);

    const double& secvz = trk.vz(), secvx = trk.vx(), secvy = trk.vy();
    const double& px = trk.px(), py = trk.py(), pz = trk.pz();

    //PAngle
    const TVector3 ptosvec(secvx-bestvx, secvy-bestvy, secvz-bestvz);
    const TVector3 secvec(px, py, pz);
    const TVector3 ptosvec2D(secvx-bestvx, secvy-bestvy, 0);
    const TVector3 secvec2D(px,py,0);

    agl[it] = std::cos(secvec.Angle(ptosvec));
    agl_abs[it] = secvec.Angle(ptosvec);
    agl2D[it] = std::cos(secvec2D.Angle(ptosvec2D));
    agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
        
    //Decay length 3D
    const SMatrixSym3D& trkCovMat = *trk.userData<reco::Vertex::CovarianceMatrix>("vertexCovariance");
    const SMatrixSym3D& totalCov = vtx.covariance() + trkCovMat;
    const SVector3 distanceVector(secvx-bestvx, secvy-bestvy, secvz-bestvz);

    dl[it] = ROOT::Math::Mag(distanceVector);
    dlerror[it] = std::sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl[it];
    dlos[it] = dl[it]/dlerror[it];

    //Decay length 2D
    const SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1), vtx.covariance(1,1), 0, 0, 0);
    const SVector6 v2(trkCovMat(0,0), trkCovMat(0,1), trkCovMat(1,1), 0, 0, 0);
    const SMatrixSym3D& totalCov2D = SMatrixSym3D(v1) + SMatrixSym3D(v2);
    const SVector3 distanceVector2D(secvx-bestvx, secvy-bestvy, 0);

    dl2D[it] = ROOT::Math::Mag(distanceVector2D);
    const double& dl2Derror = std::sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D[it];
    dlos2D[it] = dl2D[it]/dl2Derror;

    const ushort& nDau = trk.numberOfDaughters();
    if(nDau!=NDAU_) throw cms::Exception("PATCompositeAnalyzer") << "Expected " << NDAU_ << " daughters but V0 candidate has " << nDau << " daughters!" << std::endl;

    //Gen match
    if(doGenMatching_)
    {
      isSwap[it] = false;
      bool foundMom = true;
      reco::GenParticleRef genMom;
      for(ushort iDau=0; iDau<nDau; iDau++)
      {
        const auto& recDau = dynamic_cast<const pat::Muon*>(trk.daughter(iDau));
        const auto& genDau = (recDau ? recDau->genParticleRef() : reco::GenParticleRef());
        const auto& mom = findMother(genDau);
        if(!recDau || genDau.isNull() || mom.isNull()) { foundMom = false; break; }
        if(iDau==0) genMom = mom;
        else if(genMom!=mom) { foundMom = false; break; }
        if(fabs(recDau->mass() - genDau->mass())>0.01) { isSwap[it] = true; }
      }
      matchGEN[it] = foundMom;
      isSwap[it] = (foundMom ? isSwap[it] : false);
      idmom_reco[it] = (foundMom ? genMom->pdgId() : -77);
      genVec_.push_back(foundMom ? genMom : reco::GenParticleRef());
    }

    for(ushort iDau=0; iDau<nDau; iDau++)
    {
      const auto& dau = *(trk.daughter(iDau));

      ptDau[iDau][it] = dau.pt();
      pDau[iDau][it] = dau.p();
      etaDau[iDau][it] = dau.eta();
      phiDau[iDau][it] = dau.phi();
      chargeDau[iDau][it] = dau.charge();

      pid[iDau][it] = -77;
      if(doGenMatchingTOF_)
      {
        const auto& recDau = dynamic_cast<const pat::Muon*>(trk.daughter(iDau));
        const auto& genDau = (recDau ? recDau->genParticleRef() : reco::GenParticleRef());
        if(genDau.isNonnull()) { pid[iDau][it] = genDau->pdgId(); }
      }

      //trk info
      if(!twoLayerDecay_)
      {
        const auto& dtrk = dau.get<reco::TrackRef>();

        //trk quality
        trkquality[iDau][it] = (dtrk.isNonnull() ? dtrk->quality(reco::TrackBase::highPurity) : false);

        //trk dEdx
        H2dedx[iDau][it] = -999.9;
        if(dtrk.isNonnull() && dEdxHandle1.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle1.product();
          H2dedx[iDau][it] = dEdxTrack[dtrk].dEdx();
        }
        T4dedx[iDau][it] = -999.9;
        if(dtrk.isNonnull() && dEdxHandle2.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle2.product();
          T4dedx[iDau][it] = dEdxTrack[dtrk].dEdx();
        }

        //track Chi2
        trkChi[iDau][it] = (dtrk.isNonnull() ? dtrk->normalizedChi2() : 99.);

        //track pT error
        ptErr[iDau][it] = (dtrk.isNonnull() ? dtrk->ptError() : -1.);

        //trkNHits
        nhit[iDau][it] = (dtrk.isNonnull() ? dtrk->numberOfValidHits() : -1);

        //DCA
        dzos[iDau][it] = 99.;
        dxyos[iDau][it] = 99.;
        if (dtrk.isNonnull())
        {
          const double& dzbest = dtrk->dz(bestvtx);
          const double& dxybest = dtrk->dxy(bestvtx);
          const double& dzerror = std::sqrt(dtrk->dzError()*dtrk->dzError() + bestvzError*bestvzError);
          const double& dxyerror = std::sqrt(dtrk->d0Error()*dtrk->d0Error() + bestvxError*bestvyError);
          dzos[iDau][it] = dzbest/dzerror;
          dxyos[iDau][it] = dxybest/dxyerror;
        }
      }
 
      if(doMuon_)
      {
        const auto& muon = (dau.isMuon() ? *dynamic_cast<const pat::Muon*>(trk.daughter(iDau)) : pat::Muon());

        // Tight ID Muon POG Run 2
        glbmuon[iDau][it] = (dau.isMuon() ? muon.isGlobalMuon() : false);
        pfmuon[iDau][it]  = (dau.isMuon() ? muon.isPFMuon() : false);
        glbtrkchi[iDau][it] = (muon.globalTrack().isNonnull() ? muon.globalTrack()->normalizedChi2() : 99.);
        nmuonhit[iDau][it] = (muon.globalTrack().isNonnull() ? muon.globalTrack()->hitPattern().numberOfValidMuonHits() : -1);
        nmatchedst[iDau][it] = (dau.isMuon() ? muon.numberOfMatchedStations() : -1);
        npixelhit[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().numberOfValidPixelHits() : -1);
        ntrackerlayer[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -1);
        muonbestdxy[iDau][it] = (muon.muonBestTrack().isNonnull() ? muon.muonBestTrack()->dxy(bestvtx) : 99.);
        muonbestdz[iDau][it] = (muon.muonBestTrack().isNonnull() ? muon.muonBestTrack()->dz(bestvtx) : 99.);
        tightmuon[iDau][it] = (
                               glbmuon[iDau][it] &&
                               pfmuon[iDau][it] &&
                               (glbtrkchi[iDau][it] < 10.) &&
                               (nmuonhit[iDau][it] > 0) &&
                               (nmatchedst[iDau][it] > 1) &&
                               (npixelhit[iDau][it] > 0) &&
                               (ntrackerlayer[iDau][it] > 5) &&
                               (fabs(muonbestdxy[iDau][it]) < 0.2) &&
                               (fabs(muonbestdz[iDau][it]) < 0.5)
                               );

        // Soft ID Muon POG Run 2
        onestmuon[iDau][it] = (dau.isMuon() ? muon::isGoodMuon(muon, muon::SelectionType::TMOneStationTight) : false);
        npixellayer[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() : -1);
        hpmuon[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->quality(reco::TrackBase::highPurity) : false);
        muondxy[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->dxy(bestvtx) : 99.);
        muondz[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->dz(bestvtx) : 99.);
        softmuon[iDau][it] = (
                              onestmuon[iDau][it] &&
                              (ntrackerlayer[iDau][it] > 5) &&
                              (npixellayer[iDau][it] > 0) &&
                              hpmuon[iDau][it] &&
                              (fabs(muondxy[iDau][it]) < 0.3) &&
                              (fabs(muondz[iDau][it]) < 20.)
                              );

        // Hybrid Soft ID HIN PAG Run 2 PbPb
        trkmuon[iDau][it] = (dau.isMuon() ? muon.isTrackerMuon() : false);
        hybridmuon[iDau][it] = (
                                glbmuon[iDau][it] &&
                                (ntrackerlayer[iDau][it] > 5) &&
                                (npixellayer[iDau][it] > 0) &&
                                (fabs(muondxy[iDau][it]) < 0.3) &&
                                (fabs(muondz[iDau][it]) < 20.)
                                );

        // Muon Trigger Matching
        if(it==0)
        {
          trgmuon[iDau].clear();
          trgmuon[iDau] = std::vector<std::vector<UChar_t>>(filterNames_.size(), std::vector<UChar_t>(candSize, 0));
        }
        if(dau.isMuon())
        {
          for(ushort iTr=0; iTr<filterNames_.size(); iTr++)
          {
            const auto& muHLTMatchesFilter = muon.triggerObjectMatchesByFilter(filterNames_.at(iTr));
            if(muHLTMatchesFilter.size()>0) trgmuon[iDau][iTr][it] = 1;
          }
        }

        if(doMuonFull_)
        {
          nmatchedch[iDau][it] = (dau.isMuon() ? muon.numberOfMatches() : -1);
          matchedenergy[iDau][it] = (dau.isMuon() ? muon.calEnergy().hadMax : -99.);

          dx_seg[iDau][it] = 999.9;
          dy_seg[iDau][it] = 999.9;
          dxSig_seg[iDau][it] = 999.9;
          dySig_seg[iDau][it] = 999.9;
          ddxdz_seg[iDau][it] = 999.9;
          ddydz_seg[iDau][it] = 999.9;
          ddxdzSig_seg[iDau][it] = 999.9;
          ddydzSig_seg[iDau][it] = 999.9;
          const std::vector<reco::MuonChamberMatch>& muchmatches = muon.matches();
          for(ushort ich=0; ich<muchmatches.size(); ich++)
          {
            const double& x_exp = muchmatches[ich].x;
            const double& y_exp = muchmatches[ich].y;
            const double& xerr_exp = muchmatches[ich].xErr;
            const double& yerr_exp = muchmatches[ich].yErr;
            const double& dxdz_exp = muchmatches[ich].dXdZ;
            const double& dydz_exp = muchmatches[ich].dYdZ;
            const double& dxdzerr_exp = muchmatches[ich].dXdZErr;
            const double& dydzerr_exp = muchmatches[ich].dYdZErr;

            const std::vector<reco::MuonSegmentMatch>& musegmatches = muchmatches[ich].segmentMatches;
            for(ushort jseg=0; jseg<musegmatches.size(); jseg++)
            {
              const double& x_seg = musegmatches[jseg].x;
              const double& y_seg = musegmatches[jseg].y;
              const double& xerr_seg = musegmatches[jseg].xErr;
              const double& yerr_seg = musegmatches[jseg].yErr;
              const double& dxdz_seg = musegmatches[jseg].dXdZ;
              const double& dydz_seg = musegmatches[jseg].dYdZ;
              const double& dxdzerr_seg = musegmatches[jseg].dXdZErr;
              const double& dydzerr_seg = musegmatches[jseg].dYdZErr;

              const double& dseg = std::sqrt((x_seg-x_exp)*(x_seg-x_exp) + (y_seg-y_exp)*(y_seg-y_exp));
              const double& dxerr_seg = std::sqrt(xerr_seg*xerr_seg + xerr_exp*xerr_exp);
              const double& dyerr_seg = std::sqrt(yerr_seg*yerr_seg + yerr_exp*yerr_exp);
              const double& ddxdzerr_seg = std::sqrt(dxdzerr_seg*dxdzerr_seg + dxdzerr_exp*dxdzerr_exp);
              const double& ddydzerr_seg = std::sqrt(dydzerr_seg*dydzerr_seg + dydzerr_exp*dydzerr_exp);

              if(dseg < std::sqrt(dx_seg[iDau][it]*dx_seg[iDau][it] + dy_seg[iDau][it]*dy_seg[iDau][it]))
              {
                dx_seg[iDau][it] = x_seg - x_exp;
                dy_seg[iDau][it] = y_seg - y_exp;
                dxSig_seg[iDau][it] = dx_seg[iDau][it] / dxerr_seg;
                dySig_seg[iDau][it] = dy_seg[iDau][it] / dyerr_seg;
                ddxdz_seg[iDau][it] = dxdz_seg - dxdz_exp;
                ddydz_seg[iDau][it] = dydz_seg - dydz_exp;
                ddxdzSig_seg[iDau][it] = ddxdz_seg[iDau][it] / ddxdzerr_seg;
                ddydzSig_seg[iDau][it] = ddydz_seg[iDau][it] / ddydzerr_seg;
              }
            }
          }
        }
      }
    }
 
    if(twoLayerDecay_)
    {
      const auto& d = *(trk.daughter(0));
      grand_mass[it] = d.mass();
      for(ushort iGDau=0; iGDau<NGDAU_; iGDau++)
      {
        if(!d.daughter(iGDau)) continue;
        const auto& gd = *(d.daughter(iGDau));
        const auto& gdau = gd.get<reco::TrackRef>();

        //trk quality
        grand_trkquality[iGDau][it] = (gdau.isNonnull() ? gdau->quality(reco::TrackBase::highPurity) : false);

        //trk dEdx
        grand_H2dedx[iGDau][it] = -999.9;
        if(gdau.isNonnull() && dEdxHandle1.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle1.product();
          grand_H2dedx[iGDau][it] = dEdxTrack[gdau].dEdx();
        }   
        grand_T4dedx[iGDau][it] = -999.9;
        if(gdau.isNonnull() && dEdxHandle2.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle2.product();
          grand_T4dedx[iGDau][it] = dEdxTrack[gdau].dEdx();
        }

        //track pt
        grand_pt[iGDau][it] = gd.pt();

        //track momentum
        grand_p[iGDau][it] = gd.p();

        //track eta
        grand_eta[iGDau][it] = gd.eta();

        //track charge
        grand_charge[iGDau][it] = gd.charge();

        //track Chi2
        grand_trkChi[iGDau][it] = (gdau.isNonnull() ? gdau->normalizedChi2() : 99.);

        //track pT error
        grand_ptErr[iGDau][it] = (gdau.isNonnull() ? gdau->ptError() : -1.);

        //trkNHits
        grand_nhit[iGDau][it] = (gdau.isNonnull() ? gdau->numberOfValidHits() : -1);

        //DCA
        grand_dzos[iGDau][it] = 99.;
        grand_dxyos[iGDau][it] = 99.;
        if(gdau.isNonnull())
        {
          const double& gdzbest = gdau->dz(bestvtx);
          const double& gdxybest = gdau->dxy(bestvtx);
          const double& gdzerror = std::sqrt(gdau->dzError()*gdau->dzError() + bestvzError*bestvzError);
          const double& gdxyerror = std::sqrt(gdau->d0Error()*gdau->d0Error() + bestvxError*bestvyError);
          grand_dzos[iGDau][it] = gdzbest/gdzerror;
          grand_dxyos[iGDau][it] = gdxybest/gdxyerror;
        }
      }
   
      //vtxChi2
      grand_vtxChi2[it] = d.vertexChi2();
      grand_ndf[it] = d.vertexNdof();
      grand_VtxProb[it] = TMath::Prob(grand_vtxChi2[it], grand_ndf[it]);

      //PAngle
      const double& secvz = d.vz(), secvx = d.vx(), secvy = d.vy();
      const TVector3 ptosvec(secvx-bestvx, secvy-bestvy, secvz-bestvz);
      const TVector3 secvec(d.px(), d.py(), d.pz());            
      const TVector3 ptosvec2D(secvx-bestvx, secvy-bestvy, 0);
      const TVector3 secvec2D(d.px(), d.py(), 0);

      grand_agl[it] = std::cos(secvec.Angle(ptosvec));
      grand_agl_abs[it] = secvec.Angle(ptosvec);
      grand_agl2D[it] = std::cos(secvec2D.Angle(ptosvec2D));
      grand_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);

      //Decay length 3D
      const SMatrixSym3D& totalCov = vtx.covariance() + d.vertexCovariance();
      const SVector3 distanceVector(secvx-bestvx, secvy-bestvy, secvz-bestvz);

      grand_dl[it] = ROOT::Math::Mag(distanceVector);
      grand_dlerror[it] = std::sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/grand_dl[it];
      grand_dlos[it] = grand_dl[it]/grand_dlerror[it];

      //Decay length 2D
      const SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1), vtx.covariance(1,1), 0, 0, 0);
      const SVector6 v2(d.vertexCovariance(0,0), d.vertexCovariance(0,1), d.vertexCovariance(1,1), 0, 0, 0);
      const SMatrixSym3D totalCov2D = SMatrixSym3D(v1) + SMatrixSym3D(v2);
      const SVector3 distanceVector2D(secvx-bestvx, secvy-bestvy, 0);

      const double& gdl2D = ROOT::Math::Mag(distanceVector2D);
      const double& gdl2Derror = std::sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/gdl2D;
      grand_dlos2D[it] = gdl2D/gdl2Derror;
    }

    if(saveHistogram_)
    {
      for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
      {
        for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
        {
          if(pt[it]<pTBins_[ipt+1] && pt[it]>pTBins_[ipt] && y[it]<yBins_[iy+1] && y[it]>yBins_[iy])
          {
            hMassVsMVA[iy][ipt]->Fill(mva[it], mass[it]);

            if(saveAllHistogram_)
            {
              hpTVsMVA[iy][ipt]->Fill(mva[it], pt[it]);
              hetaVsMVA[iy][ipt]->Fill(mva[it], eta[it]);
              hyVsMVA[iy][ipt]->Fill(mva[it], y[it]);
              hVtxProbVsMVA[iy][ipt]->Fill(mva[it], VtxProb[it]);
              h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl[it]);
              h3DPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl_abs[it]);
              h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl2D[it]);
              h2DPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl2D_abs[it]);
              h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it], dlos[it]);
              h3DDecayLengthVsMVA[iy][ipt]->Fill(mva[it], dl[it]);
              h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it], dlos2D[it]);
              h2DDecayLengthVsMVA[iy][ipt]->Fill(mva[it], dl2D[it]);
              for (ushort iDau=0; iDau<NDAU_; iDau++)
              {
                hzDCASignificanceDaugtherVsMVA[iDau][iy][ipt]->Fill(mva[it], dzos[iDau][it]);
                hxyDCASignificanceDaugtherVsMVA[iDau][iy][ipt]->Fill(mva[it], dxyos[iDau][it]);
                hNHitDVsMVA[iDau][iy][ipt]->Fill(mva[it], nhit[iDau][it]);
                hpTDVsMVA[iDau][iy][ipt]->Fill(mva[it], ptDau[iDau][it]);
                hpTerrDVsMVA[iDau][iy][ipt]->Fill(mva[it], ptErr[iDau][it]/ptDau[iDau][it]);
                hEtaDVsMVA[iDau][iy][ipt]->Fill(mva[it], etaDau[iDau][it]);
                hdedxHarmonic2DVsMVA[iDau][iy][ipt]->Fill(mva[it], H2dedx[iDau][it]);
                hdedxHarmonic2DVsP[iDau][iy][ipt]->Fill(pDau[iDau][it], H2dedx[iDau][it]);
              }
            }
          }
        }
      }
    }
  }
}


void
PATCompositeTreeProducer::fillGEN(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<GenEventInfoProduct> geninfo;
  iEvent.getByToken(tok_genInfo_, geninfo);
  weight_gen = (geninfo.isValid() ? geninfo->weight() : -1.0);

  edm::Handle<reco::GenParticleCollection> genpars;
  iEvent.getByToken(tok_genParticle_, genpars);
  if(!genpars.isValid()) throw cms::Exception("PATCompositeAnalyzer") << "Gen matching cannot be done without Gen collection!" << std::endl;

  candSize_gen = 0;
  for(uint idx=0; idx<genpars->size(); ++idx)
  {
    const auto& trk = reco::GenParticleRef(genpars, idx);

    if (trk.isNull()) continue; //check gen particle ref
    if(!hasDaughters(trk, PID_dau_)) continue; //check if has the daughters

    pt_gen[candSize_gen] = trk->pt();
    eta_gen[candSize_gen] = trk->eta();
    y_gen[candSize_gen] = trk->rapidity();
    pid_gen[candSize_gen] = trk->pdgId();
    status_gen[candSize_gen] = trk->status();

    const auto& recIt = std::find(genVec_.begin(), genVec_.end(), trk);
    idxrec_gen[candSize_gen] = (recIt!=genVec_.end() ? std::distance(genVec_.begin(), recIt) : -1);

    const auto& mom = findMother(trk);
    idmom_gen[candSize_gen] = (mom.isNonnull() ? mom->pdgId() : -77);

    if(decayInGen_)
    {
      for(ushort iDau=0; iDau<NDAU_; iDau++)
      {
        const auto& Dd = findLastPar(trk->daughterRef(iDau));
        idDau_gen[iDau][candSize_gen] = (Dd.isNonnull() ? Dd->pdgId() : 99999);
        chargeDau_gen[iDau][candSize_gen] = (Dd.isNonnull() ? Dd->charge() : 9);
        ptDau_gen[iDau][candSize_gen] = (Dd.isNonnull() ? Dd->pt() : -1.);
        etaDau_gen[iDau][candSize_gen] = (Dd.isNonnull() ? Dd->eta() : 9.);
        phiDau_gen[iDau][candSize_gen] = (Dd.isNonnull() ? Dd->phi() : 9.);
      }
    }

    candSize_gen++;
  }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
PATCompositeTreeProducer::beginJob()
{
  TH1D::SetDefaultSumw2();

  // Check inputs
  if((!threeProngDecay_ && NDAU_!=2) || (threeProngDecay_ && NDAU_!=3))
  {
    throw cms::Exception("PATCompositeAnalyzer") << "Want threeProngDecay but PID daughter vector size is: " << NDAU_ << " !" << std::endl;
  }
  if(!doRecoNtuple_ && !doGenNtuple_) throw cms::Exception("PATCompositeAnalyzer") << "No output for either RECO or GEN!! Fix config!!" << std::endl;
  if(twoLayerDecay_ && doMuon_) throw cms::Exception("PATCompositeAnalyzer") << "Muons cannot be coming from two layer decay!! Fix config!!" << std::endl;

  if(saveHistogram_) initHistogram();
  if(saveTree_) initTree();
}


void
PATCompositeTreeProducer::initHistogram()
{
  for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
  {
    for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
    {
      hMassVsMVA[iy][ipt] = fs->make<TH2F>(Form("hMassVsMVA_y%d_pt%d",iy,ipt),";mva;mass(GeV)",100,-1.,1.,massHistBins_,massHistPeak_-massHistWidth_,massHistPeak_+massHistWidth_);
      if(saveAllHistogram_)
      {
        hpTVsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTVsMVA_y%d_pt%d",iy,ipt),";mva;pT;",100,-1,1,100,0,10);
        hetaVsMVA[iy][ipt] = fs->make<TH2F>(Form("hetaVsMVA_y%d_pt%d",iy,ipt),";mva;eta;",100,-1.,1.,40,-4,4);
        hyVsMVA[iy][ipt] = fs->make<TH2F>(Form("hyVsMVA_y%d_pt%d",iy,ipt),";mva;y;",100,-1.,1.,40,-4,4);
        hVtxProbVsMVA[iy][ipt] = fs->make<TH2F>(Form("hVtxProbVsMVA_y%d_pt%d",iy,ipt),";mva;VtxProb;",100,-1.,1.,100,0,1);
        h3DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DCosPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;3DCosPointingAngle;",100,-1.,1.,100,-1,1);
        h3DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;3DPointingAngle;",100,-1.,1.,50,-3.14,3.14);
        h2DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DCosPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;2DCosPointingAngle;",100,-1.,1.,100,-1,1);
        h2DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;2DPointingAngle;",100,-1.,1.,50,-3.14,3.14);
        h3DDecayLengthSignificanceVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDecayLengthSignificanceVsMVA_y%d_pt%d",iy,ipt),";mva;3DDecayLengthSignificance;",100,-1.,1.,300,0,30);
        h2DDecayLengthSignificanceVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDecayLengthSignificanceVsMVA_y%d_pt%d",iy,ipt),";mva;2DDecayLengthSignificance;",100,-1.,1.,300,0,30);
        h3DDecayLengthVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDecayLengthVsMVA_y%d_pt%d",iy,ipt),";mva;3DDecayLength;",100,-1.,1.,300,0,30);
        h2DDecayLengthVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDecayLengthVsMVA_y%d_pt%d",iy,ipt),";mva;2DDecayLength;",100,-1.,1.,300,0,30);
        for(ushort d=1; d<=NDAU_; d++)
        {
          hzDCASignificanceDaugtherVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hzDCASignificanceDaugther%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;zDCASignificanceDaugther%d;",d),100,-1.,1.,100,-10,10);
          hxyDCASignificanceDaugtherVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificanceDaugther%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;xyDCASignificanceDaugther%d;",d),100,-1.,1.,100,-10,10);
          hNHitDVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hNHitD%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;NHitD%d;",d),100,-1.,1.,100,0,100);
          hpTDVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hpTD%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;pTD%d;",d),100,-1.,1.,100,0,10);
          hpTerrDVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hpTerrD%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;pTerrD%d;",d),100,-1.,1.,50,0,0.5);
          hEtaDVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hEtaD%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;EtaD%d;",d),100,-1.,1.,40,-4,4);
          if(useDeDxData_)
          {
            hdedxHarmonic2DVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;dedxHarmonic2D%d;",d),100,-1.,1.,100,0,10);
            hdedxHarmonic2DVsP[d-1][iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D%dVsP_y%d_pt%d",d,iy,ipt),Form(";p (GeV);dedxHarmonic2D%d",d),100,0,10,100,0,10);
          }
        }
      }
    }
  }
}


void 
PATCompositeTreeProducer::initTree()
{ 
  PATCompositeNtuple = fs->make< TTree>("VertexCompositeNtuple","VertexCompositeNtuple");

  if(doRecoNtuple_)
  {
    // Event info
    
    PATCompositeNtuple->Branch("RunNb",&runNb,"RunNb/i");
    PATCompositeNtuple->Branch("LSNb",&lsNb,"LSNb/i");
    PATCompositeNtuple->Branch("EventNb",&eventNb,"EventNb/i");
    PATCompositeNtuple->Branch("nPV",&nPV,"nPV/S");
    PATCompositeNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
    PATCompositeNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
    PATCompositeNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");
    PATCompositeNtuple->Branch("candSize",&candSize,"candSize/i");
    if(isCentrality_) 
    {
      PATCompositeNtuple->Branch("centrality",&centrality,"centrality/S");
      PATCompositeNtuple->Branch("Npixel",&Npixel,"Npixel/I");
      PATCompositeNtuple->Branch("HFsumETPlus",&HFsumETPlus,"HFsumETPlus/F");
      PATCompositeNtuple->Branch("HFsumETMinus",&HFsumETMinus,"HFsumETMinus/F");
      PATCompositeNtuple->Branch("ZDCPlus",&ZDCPlus,"ZDCPlus/F");
      PATCompositeNtuple->Branch("ZDCMinus",&ZDCMinus,"ZDCMinus/F");
      PATCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
    }
    if(isEventPlane_) {
      PATCompositeNtuple->Branch("ephfpAngle",ephfpAngle,"ephfpAngle[3]/F");
      PATCompositeNtuple->Branch("ephfmAngle",ephfmAngle,"ephfmAngle[3]/F");
      PATCompositeNtuple->Branch("eptrackmidAngle",eptrackmidAngle,"eptrackmidAngle[3]/F");
      PATCompositeNtuple->Branch("ephfpQ",ephfpQ,"ephfpQ[3]/F");
      PATCompositeNtuple->Branch("ephfmQ",ephfmQ,"ephfmQ[3]/F");
      PATCompositeNtuple->Branch("eptrackmidQ",eptrackmidQ,"eptrackmidQ[3]/F");
      PATCompositeNtuple->Branch("ephfpSumW",&ephfpSumW,"ephfpSumW/F");
      PATCompositeNtuple->Branch("ephfmSumW",&ephfmSumW,"ephfmSumW/F");
      PATCompositeNtuple->Branch("eptrackmidSumW",&eptrackmidSumW,"eptrackmidSumW/F");
    }
    PATCompositeNtuple->Branch("trigPrescale",trigPrescale,Form("trigPrescale[%d]/S",NTRG_));
    PATCompositeNtuple->Branch("trigHLT",trigHLT,Form("trigHLT[%d]/O",NTRG_));
    PATCompositeNtuple->Branch("evtSel",evtSel,Form("evtSel[%d]/O",NSEL_));

    // particle info
    PATCompositeNtuple->Branch("pT",pt,"pT[candSize]/F");
    PATCompositeNtuple->Branch("eta",eta,"eta[candSize]/F");
    PATCompositeNtuple->Branch("phi",phi,"phi[candSize]/F");
    PATCompositeNtuple->Branch("mass",mass,"mass[candSize]/F");
    PATCompositeNtuple->Branch("y",y,"y[candSize]/F");
    if(useAnyMVA_) PATCompositeNtuple->Branch("mva",mva,"mva[candSize]/F");

    if(!isSkimMVA_)
    {
      //Composite candidate info RECO
      PATCompositeNtuple->Branch("flavor",flavor,"flavor[candSize]/F");
      PATCompositeNtuple->Branch("VtxProb",VtxProb,"VtxProb[candSize]/F");
      PATCompositeNtuple->Branch("3DCosPointingAngle",agl,"3DCosPointingAngle[candSize]/F");
      PATCompositeNtuple->Branch("3DPointingAngle",agl_abs,"3DPointingAngle[candSize]/F");
      PATCompositeNtuple->Branch("2DCosPointingAngle",agl2D,"2DCosPointingAngle[candSize]/F");
      PATCompositeNtuple->Branch("2DPointingAngle",agl2D_abs,"2DPointingAngle[candSize]/F");
      PATCompositeNtuple->Branch("3DDecayLengthSignificance",dlos,"3DDecayLengthSignificance[candSize]/F");
      PATCompositeNtuple->Branch("3DDecayLength",dl,"3DDecayLength[candSize]/F");
      PATCompositeNtuple->Branch("3DDecayLengthError",dlerror,"3DDecayLengthError[candSize]/F");
      PATCompositeNtuple->Branch("2DDecayLengthSignificance",dlos2D,"2DDecayLengthSignificance[candSize]/F");
      PATCompositeNtuple->Branch("2DDecayLength",dl2D,"2DDecayLength[candSize]/F");

      if(doGenMatching_)
      {
        PATCompositeNtuple->Branch("isSwap",isSwap,"isSwap[candSize]/O");
        PATCompositeNtuple->Branch("idmom_reco",idmom_reco,"idmom_reco[candSize]/I");
        PATCompositeNtuple->Branch("matchGEN",matchGEN,"matchGEN[candSize]/O");
      }
 
      if(doGenMatchingTOF_)
      {
        for(ushort iDau=1; iDau<=NDAU_; iDau++)
        {
          PATCompositeNtuple->Branch(Form("PIDD%d",iDau),pid[iDau-1],Form("PIDD%d[candSize]/I",iDau));
        }
      }

      //daugther & grand daugther info
      if(twoLayerDecay_)
      {
        PATCompositeNtuple->Branch("massDaugther1",grand_mass,"massDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("VtxProbDaugther1",grand_VtxProb,"VtxProbDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DCosPointingAngleDaugther1",grand_agl,"3DCosPointingAngleDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DPointingAngleDaugther1",grand_agl_abs,"3DPointingAngleDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("2DCosPointingAngleDaugther1",grand_agl2D,"2DCosPointingAngleDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("2DPointingAngleDaugther1",grand_agl2D_abs,"2DPointingAngleDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DDecayLengthSignificanceDaugther1",grand_dlos,"3DDecayLengthSignificanceDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DDecayLengthDaugther1",grand_dl,"3DDecayLengthDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DDecayLengthErrorDaugther1",grand_dlerror,"3DDecayLengthErrorDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("2DDecayLengthSignificanceDaugther1",grand_dlos2D,"2DDecayLengthSignificanceDaugther1[candSize]/F");
        for(ushort iGDau=1; iGDau<=NGDAU_; iGDau++)
        {
          PATCompositeNtuple->Branch(Form("zDCASignificanceGrandDaugther%d",iGDau),grand_dzos[iGDau-1],Form("zDCASignificanceGrandDaugther%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("xyDCASignificanceGrandDaugther%d",iGDau),grand_dxyos[iGDau-1],Form("xyDCASignificanceGrandDaugther%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("NHitGrandD%d",iGDau),grand_nhit[iGDau-1],Form("NHitGrandD%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("HighPurityGrandDaugther%d",iGDau),grand_trkquality[iGDau-1],Form("HighPurityGrandDaugther%d[candSize]/O",iGDau));
          PATCompositeNtuple->Branch(Form("pTGrandD%d",iGDau),grand_pt[iGDau-1],Form("pTGrandD%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("pTerrGrandD%d",iGDau),grand_ptErr[iGDau-1],Form("pTerrGrandD%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("EtaGrandD%d",iGDau),grand_eta[iGDau-1],Form("EtaGrandD%d[candSize]/F",iGDau));
          if(useDeDxData_)
          {
            PATCompositeNtuple->Branch(Form("dedxPixelHarmonic2GrandD%d",iGDau),grand_T4dedx[iGDau-1],Form("dedxPixelHarmonic2GrandD%d[candSize]/F",iGDau));
            PATCompositeNtuple->Branch(Form("dedxHarmonic2GrandD%d",iGDau),grand_H2dedx[iGDau-1],Form("dedxHarmonic2GrandD%d[candSize]/F",iGDau));
          }
        }
      }
      for(ushort iDau=1; iDau<=NDAU_; iDau++)
      {
        PATCompositeNtuple->Branch(Form("zDCASignificanceDaugther%d",iDau),dzos[iDau-1],Form("zDCASignificanceDaugther%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("xyDCASignificanceDaugther%d",iDau),dxyos[iDau-1],Form("xyDCASignificanceDaugther%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("NHitD%d",iDau),nhit[iDau-1],Form("NHitD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("HighPurityDaugther%d",iDau),trkquality[iDau-1],Form("HighPurityDaugther%d[candSize]/O",iDau));
        PATCompositeNtuple->Branch(Form("pTD%d",iDau),ptDau[iDau-1],Form("pTD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("pTerrD%d",iDau),ptErr[iDau-1],Form("pTerrD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("EtaD%d",iDau),etaDau[iDau-1],Form("EtaD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("PhiD%d",iDau),phiDau[iDau-1],Form("PhiD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("chargeD%d",iDau),chargeDau[iDau-1],Form("chargeD%d[candSize]/S",iDau));
        if(useDeDxData_)
        {
          PATCompositeNtuple->Branch(Form("dedxPixelHarmonic2D%d",iDau),T4dedx[iDau-1],Form("dedxPixelHarmonic2D%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("dedxHarmonic2D%d",iDau),H2dedx[iDau-1],Form("dedxHarmonic2D%d[candSize]/F",iDau));
        }
      }
 
      if(doMuon_)
      {
        for(ushort iDau=1; iDau<=NDAU_; iDau++)
        {
          if(fabs(PID_dau_[iDau-1])!=13) continue;
          PATCompositeNtuple->Branch(Form("OneStMuon%d",iDau),onestmuon[iDau-1],Form("OneStMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("PFMuon%d",iDau),pfmuon[iDau-1],Form("PFMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("GlbMuon%d",iDau),glbmuon[iDau-1],Form("GlbMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("trkMuon%d",iDau),trkmuon[iDau-1],Form("trkMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("tightMuon%d",iDau),tightmuon[iDau-1],Form("tightMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("softMuon%d",iDau),softmuon[iDau-1],Form("softMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("hybridMuon%d",iDau),hybridmuon[iDau-1],Form("hybridMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("HPMuon%d",iDau),hpmuon[iDau-1],Form("hybridMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("trigMuon%d",iDau),&(trgmuon[iDau-1]));
          PATCompositeNtuple->Branch(Form("nMatchedStationD%d",iDau),nmatchedst[iDau-1],Form("nMatchedStationD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("nTrackerLayerD%d",iDau),ntrackerlayer[iDau-1],Form("nTrackerLayerD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("nPixelLayerD%d",iDau),npixellayer[iDau-1],Form("nPixelLayerD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("nPixelHitD%d",iDau),npixelhit[iDau-1],Form("nPixelHitD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("nMuonHitD%d",iDau),nmuonhit[iDau-1],Form("nMuonHitD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("GlbTrkChiD%d",iDau),glbtrkchi[iDau-1],Form("GlbTrkChiD%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("muondXYD%d",iDau),muonbestdxy[iDau-1],Form("muondXYD%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("muondZD%d",iDau),muonbestdz[iDau-1],Form("muondZD%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("dXYD%d",iDau),muondxy[iDau-1],Form("dXYD%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("dZD%d",iDau),muondz[iDau-1],Form("dZD%d[candSize]/F",iDau));
          if(doMuonFull_)
          {
            PATCompositeNtuple->Branch(Form("nMatchedChamberD%d",iDau),nmatchedch[iDau-1],Form("nMatchedChamberD%d[candSize]/S",iDau));
            PATCompositeNtuple->Branch(Form("EnergyDepositionD%d",iDau),matchedenergy[iDau-1],Form("EnergyDepositionD%d[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("dx%d_seg",iDau),        dx_seg[iDau-1], Form("dx%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("dy%d_seg",iDau),        dy_seg[iDau-1], Form("dy%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("dxSig%d_seg",iDau),     dxSig_seg[iDau-1], Form("dxSig%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("dySig%d_seg",iDau),     dySig_seg[iDau-1], Form("dySig%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("ddxdz%d_seg",iDau),     ddxdz_seg[iDau-1], Form("ddxdz%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("ddydz%d_seg",iDau),     ddydz_seg[iDau-1], Form("ddydz%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("ddxdzSig%d_seg",iDau),  ddxdzSig_seg[iDau-1], Form("ddxdzSig%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("ddydzSig%d_seg",iDau),  ddydzSig_seg[iDau-1], Form("ddydzSig%d_seg[candSize]/F",iDau));
          }
        }
      }
    }
  } // doRecoNtuple_

  if(doGenNtuple_)
  {
    PATCompositeNtuple->Branch("weight_gen",&weight_gen,"weight_gen/F");
    PATCompositeNtuple->Branch("candSize_gen",&candSize_gen,"candSize_gen/i");
    PATCompositeNtuple->Branch("pT_gen",pt_gen,"pT_gen[candSize_gen]/F");
    PATCompositeNtuple->Branch("eta_gen",eta_gen,"eta_gen[candSize_gen]/F");
    PATCompositeNtuple->Branch("y_gen",y_gen,"y_gen[candSize_gen]/F");
    PATCompositeNtuple->Branch("status_gen",status_gen,"status_gen[candSize_gen]/S");
    PATCompositeNtuple->Branch("PID_gen",pid_gen,"PID_gen[candSize_gen]/I");
    PATCompositeNtuple->Branch("MotherID_gen",idmom_gen,"MotherID_gen[candSize_gen]/I");
    PATCompositeNtuple->Branch("RecIdx_gen",idxrec_gen,"RecIdx_gen[candSize_gen]/S");

    if(decayInGen_)
    {
      for(ushort iDau=1; iDau<=NDAU_; iDau++)
      {
        PATCompositeNtuple->Branch(Form("PIDD%d_gen",iDau),idDau_gen[iDau-1],Form("PIDD%d_gen[candSize_gen]/I",iDau));
        PATCompositeNtuple->Branch(Form("chargeD%d_gen",iDau),chargeDau_gen[iDau-1],Form("chargeD%d_gen[candSize_gen]/S",iDau));
        PATCompositeNtuple->Branch(Form("pTD%d_gen",iDau),ptDau_gen[iDau-1],Form("pTD%d_gen[candSize_gen]/F",iDau));
        PATCompositeNtuple->Branch(Form("EtaD%d_gen",iDau),etaDau_gen[iDau-1],Form("EtaD%d_gen[candSize_gen]/F",iDau));
        PATCompositeNtuple->Branch(Form("PhiD%d_gen",iDau),phiDau_gen[iDau-1],Form("PhiD%d_gen[candSize_gen]/F",iDau));
      }
    }
  }
}


//--------------------------------------------------------------------------------------------------
void 
PATCompositeTreeProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = true;
  EDConsumerBase::Labels triggerResultsLabel;
  EDConsumerBase::labelsForToken(tok_triggerResults_, triggerResultsLabel);
  hltPrescaleProvider_.init(iRun, iSetup, triggerResultsLabel.process, changed);
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
PATCompositeTreeProducer::endJob()
{
}


reco::GenParticleRef
PATCompositeTreeProducer::findLastPar(const reco::GenParticleRef& genParRef)
{
  if(genParRef.isNull()) return genParRef;
  reco::GenParticleRef genLastParRef = genParRef;
  const int& pdg_OLD = genParRef->pdgId();
  while(genLastParRef->numberOfDaughters()>0 && genLastParRef->daughterRef(0)->pdgId()==pdg_OLD)
  {
    genLastParRef = genLastParRef->daughterRef(0);
  }
  return genLastParRef;
}


reco::GenParticleRef
PATCompositeTreeProducer::findMother(const reco::GenParticleRef& genParRef)
{
  if(genParRef.isNull()) return genParRef;
  reco::GenParticleRef genMomRef = genParRef;
  int pdg = genParRef->pdgId(); const int pdg_OLD = pdg;
  while(pdg==pdg_OLD && genMomRef->numberOfMothers()>0)
  {
    genMomRef = genMomRef->motherRef(0);
    pdg = genMomRef->pdgId();
  }
  if(pdg==pdg_OLD) genMomRef = reco::GenParticleRef();
  return genMomRef;
}


bool
PATCompositeTreeProducer::hasDaughters(const reco::GenParticleRef& genParRef, const std::vector<int>& PID_dau)
{
  bool hasDau = false;
  if(genParRef->numberOfDaughters()==PID_dau.size())
  {
    std::vector<int> PIDvec = PID_dau;
    for(ushort iDau=0; iDau<genParRef->numberOfDaughters(); iDau++)
    {
      const auto& dau = *(genParRef->daughter(iDau));
      bool found = false;
      for(ushort iPID=0; iPID<PIDvec.size(); iPID++)
      {
        if(fabs(dau.pdgId())==PIDvec[iPID])
        {
          PIDvec.erase(PIDvec.begin()+iPID);
          found=true;
          break;
        }
      }
      if(!found) break;
    }
    hasDau = (PIDvec.size()==0);
  }
  return hasDau;
}


//define this as a plug-in
DEFINE_FWK_MODULE(PATCompositeTreeProducer);
