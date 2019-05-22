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

class PATCompositeNtupleProducer : public edm::EDAnalyzer {
public:
  explicit PATCompositeNtupleProducer(const edm::ParameterSet&);
  ~PATCompositeNtupleProducer();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void initHistogram();
  virtual void initTree();
  reco::GenParticleRef findMother(const reco::GenParticleRef&);

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
  bool doGenMatching_;
  bool doGenMatchingTOF_;
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
  float ephfpQ[3];
  float ephfmQ[3];
  float ephfpSumW;
  float ephfmSumW;

  //Composite candidate info
  float mva;
  float pt;
  float eta;
  float phi;
  float flavor;
  float y;
  float mass;
  float VtxProb;
  float dlos;
  float dl;
  float dlerror;
  float agl;
  float vtxChi2;
  float ndf;
  float agl_abs;
  float agl2D;
  float agl2D_abs;
  float dlos2D;
  float dl2D;
  bool isSwap;
  bool matchGEN;
  int idmom_reco;
    
  //dau candidate info
  float grand_mass;
  float grand_VtxProb;
  float grand_dlos;
  float grand_dl;
  float grand_dlerror;
  float grand_agl;
  float grand_vtxChi2;
  float grand_ndf;
  float grand_agl_abs;
  float grand_agl2D;
  float grand_agl2D_abs;
  float grand_dlos2D;

  //dau info
  float dzos[MAXDAU];
  float dxyos[MAXDAU];
  float nhit[MAXDAU];
  bool  trkquality[MAXDAU];
  float ptDau[MAXDAU];
  float ptErr[MAXDAU];
  float pDau[MAXDAU];
  float etaDau[MAXDAU];
  float phiDau[MAXDAU];
  short chargeDau[MAXDAU];
  int   pid[MAXDAU];
  float tof[MAXDAU];
  float H2dedx[MAXDAU];
  float T4dedx[MAXDAU];
  float trkChi[MAXDAU];
   
  //grand-dau info
  float grand_dzos[MAXGDAU];
  float grand_dxyos[MAXGDAU];
  float grand_nhit[MAXGDAU];
  bool  grand_trkquality[MAXGDAU];
  float grand_pt[MAXGDAU];
  float grand_ptErr[MAXGDAU];
  float grand_p[MAXGDAU];
  float grand_eta[MAXGDAU];
  short grand_charge[MAXGDAU];
  float grand_H2dedx[MAXGDAU];
  float grand_T4dedx[MAXGDAU];
  float grand_trkChi[MAXGDAU];
    
  //dau muon info
  bool  onestmuon[MAXDAU];
  bool  pfmuon[MAXDAU];
  bool  glbmuon[MAXDAU];
  bool  trkmuon[MAXDAU];
  bool  tightmuon[MAXDAU];
  bool  softmuon[MAXDAU];
  bool  hybridmuon[MAXDAU];
  bool  hpmuon[MAXDAU];
  bool  trgmuon[MAXDAU][MAXTRG];
  short nmatchedst[MAXDAU];
  short ntrackerlayer[MAXDAU];
  short npixellayer[MAXDAU];
  short npixelhit[MAXDAU];
  short nmuonhit[MAXDAU];
  float glbtrkchi[MAXDAU];
  float muonbestdxy[MAXDAU];
  float muonbestdz[MAXDAU];
  float muondxy[MAXDAU];
  float muondz[MAXDAU];
  short nmatchedch[MAXDAU];
  float matchedenergy[MAXDAU];
  float dx_seg[MAXDAU];
  float dy_seg[MAXDAU];
  float dxSig_seg[MAXDAU];
  float dySig_seg[MAXDAU];
  float ddxdz_seg[MAXDAU];
  float ddydz_seg[MAXDAU];
  float ddxdzSig_seg[MAXDAU];
  float ddydzSig_seg[MAXDAU];

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

PATCompositeNtupleProducer::PATCompositeNtupleProducer(const edm::ParameterSet& iConfig) :
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
  twoLayerDecay_ = iConfig.getUntrackedParameter<bool>("twoLayerDecay");
  threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
  doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
  doGenMatchingTOF_ = iConfig.getUntrackedParameter<bool>("doGenMatchingTOF");
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


PATCompositeNtupleProducer::~PATCompositeNtupleProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PATCompositeNtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  if(doRecoNtuple_) fillRECO(iEvent,iSetup);
}


void
PATCompositeNtupleProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    ephfpQ[0] = (eventplanes.isValid() ? (*eventplanes)[0].q(2) : -99.);
    ephfpQ[1] = (eventplanes.isValid() ? (*eventplanes)[6].q(2) : -99.);
    ephfpQ[2] = (eventplanes.isValid() ? (*eventplanes)[13].q(2) : -99.);

    ephfmQ[0] = (eventplanes.isValid() ? (*eventplanes)[1].q(2) : -99.);
    ephfmQ[1] = (eventplanes.isValid() ? (*eventplanes)[7].q(2) : -99.);
    ephfmQ[2] = (eventplanes.isValid() ? (*eventplanes)[14].q(2) : -99.);

    ephfpSumW = (eventplanes.isValid() ? (*eventplanes)[6].sumw() : -99.);
    ephfmSumW = (eventplanes.isValid() ? (*eventplanes)[7].sumw() : -99.);
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
  for(uint it=0; it<candSize; ++it)
  { 
    const auto& trk = (*v0candidates)[it];

    pt = trk.pt();
    eta = trk.eta();
    phi = trk.phi();
    mass = trk.mass();
    y = trk.rapidity();
    flavor = (trk.pdgId()!=0 ? trk.pdgId()/fabs(trk.pdgId()) : 0.);

    mva = (useAnyMVA_ ? (*mvavalues)[it] : 0.0);

    //vtxChi2
    vtxChi2 = trk.userFloat("vertexChi2");
    ndf = trk.userFloat("vertexNdof");
    VtxProb = TMath::Prob(vtxChi2,ndf);

    const double& secvz = trk.vz(), secvx = trk.vx(), secvy = trk.vy();
    const double& px = trk.px(), py = trk.py(), pz = trk.pz();

    //PAngle
    const TVector3 ptosvec(secvx-bestvx, secvy-bestvy, secvz-bestvz);
    const TVector3 secvec(px, py, pz);
    const TVector3 ptosvec2D(secvx-bestvx, secvy-bestvy, 0);
    const TVector3 secvec2D(px,py,0);

    agl = std::cos(secvec.Angle(ptosvec));
    agl_abs = secvec.Angle(ptosvec);
    agl2D = std::cos(secvec2D.Angle(ptosvec2D));
    agl2D_abs = secvec2D.Angle(ptosvec2D);
        
    //Decay length 3D
    const SMatrixSym3D& trkCovMat = *trk.userData<reco::Vertex::CovarianceMatrix>("vertexCovariance");
    const SMatrixSym3D& totalCov = vtx.covariance() + trkCovMat;
    const SVector3 distanceVector(secvx-bestvx, secvy-bestvy, secvz-bestvz);

    dl = ROOT::Math::Mag(distanceVector);
    dlerror = std::sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
    dlos = dl/dlerror;

    //Decay length 2D
    const SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1), vtx.covariance(1,1), 0, 0, 0);
    const SVector6 v2(trkCovMat(0,0), trkCovMat(0,1), trkCovMat(1,1), 0, 0, 0);
    const SMatrixSym3D& totalCov2D = SMatrixSym3D(v1) + SMatrixSym3D(v2);
    const SVector3 distanceVector2D(secvx-bestvx, secvy-bestvy, 0);

    dl2D = ROOT::Math::Mag(distanceVector2D);
    const double& dl2Derror = std::sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D;
    dlos2D = dl2D/dl2Derror;

    const ushort& nDau = trk.numberOfDaughters();
    if(nDau!=NDAU_) throw cms::Exception("PATCompositeAnalyzer") << "Expected " << NDAU_ << " daughters but V0 candidate has " << nDau << " daughters!" << std::endl;

    //Gen match
    if(doGenMatching_)
    {
      isSwap = false;
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
        if(fabs(recDau->mass() - genDau->mass())>0.01) { isSwap = true; }
      }
      matchGEN = foundMom;
      isSwap = (foundMom ? isSwap : false);
      idmom_reco = (foundMom ? genMom->pdgId() : -77);
    }

    for(ushort iDau=0; iDau<nDau; iDau++)
    {
      const auto& dau = *(trk.daughter(iDau));

      ptDau[iDau] = dau.pt();
      pDau[iDau] = dau.p();
      etaDau[iDau] = dau.eta();
      phiDau[iDau] = dau.phi();
      chargeDau[iDau] = dau.charge();

      pid[iDau] = -77;
      if(doGenMatchingTOF_)
      {
        const auto& recDau = dynamic_cast<const pat::Muon*>(trk.daughter(iDau));
        const auto& genDau = (recDau ? recDau->genParticleRef() : reco::GenParticleRef());
        if(genDau.isNonnull()) { pid[iDau] = genDau->pdgId(); }
      }

      //trk info
      if(!twoLayerDecay_)
      {
        const auto& dtrk = dau.get<reco::TrackRef>();

        //trk quality
        trkquality[iDau] = (dtrk.isNonnull() ? dtrk->quality(reco::TrackBase::highPurity) : false);

        //trk dEdx
        H2dedx[iDau] = -999.9;
        if(dtrk.isNonnull() && dEdxHandle1.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle1.product();
          H2dedx[iDau] = dEdxTrack[dtrk].dEdx();
        }
        T4dedx[iDau] = -999.9;
        if(dtrk.isNonnull() && dEdxHandle2.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle2.product();
          T4dedx[iDau] = dEdxTrack[dtrk].dEdx();
        }

        //track Chi2
        trkChi[iDau] = (dtrk.isNonnull() ? dtrk->normalizedChi2() : 99.);

        //track pT error
        ptErr[iDau] = (dtrk.isNonnull() ? dtrk->ptError() : -1.);

        //trkNHits
        nhit[iDau] = (dtrk.isNonnull() ? dtrk->numberOfValidHits() : -1);

        //DCA
        dzos[iDau] = 99.;
        dxyos[iDau] = 99.;
        if (dtrk.isNonnull())
        {
          const double& dzbest = dtrk->dz(bestvtx);
          const double& dxybest = dtrk->dxy(bestvtx);
          const double& dzerror = std::sqrt(dtrk->dzError()*dtrk->dzError() + bestvzError*bestvzError);
          const double& dxyerror = std::sqrt(dtrk->d0Error()*dtrk->d0Error() + bestvxError*bestvyError);
          dzos[iDau] = dzbest/dzerror;
          dxyos[iDau] = dxybest/dxyerror;
        }
      }
 
      if(doMuon_)
      {
        const auto& muon = (dau.isMuon() ? *dynamic_cast<const pat::Muon*>(trk.daughter(iDau)) : pat::Muon());

        // Tight ID Muon POG Run 2
        glbmuon[iDau] = (dau.isMuon() ? muon.isGlobalMuon() : false);
        pfmuon[iDau]  = (dau.isMuon() ? muon.isPFMuon() : false);
        glbtrkchi[iDau] = (muon.globalTrack().isNonnull() ? muon.globalTrack()->normalizedChi2() : 99.);
        nmuonhit[iDau] = (muon.globalTrack().isNonnull() ? muon.globalTrack()->hitPattern().numberOfValidMuonHits() : -1);
        nmatchedst[iDau] = (dau.isMuon() ? muon.numberOfMatchedStations() : -1);
        npixelhit[iDau] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().numberOfValidPixelHits() : -1);
        ntrackerlayer[iDau] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -1);
        muonbestdxy[iDau] = (muon.muonBestTrack().isNonnull() ? muon.muonBestTrack()->dxy(bestvtx) : 99.);
        muonbestdz[iDau] = (muon.muonBestTrack().isNonnull() ? muon.muonBestTrack()->dz(bestvtx) : 99.);
        tightmuon[iDau] = (
                               glbmuon[iDau] &&
                               pfmuon[iDau] &&
                               (glbtrkchi[iDau] < 10.) &&
                               (nmuonhit[iDau] > 0) &&
                               (nmatchedst[iDau] > 1) &&
                               (npixelhit[iDau] > 0) &&
                               (ntrackerlayer[iDau] > 5) &&
                               (fabs(muonbestdxy[iDau]) < 0.2) &&
                               (fabs(muonbestdz[iDau]) < 0.5)
                               );

        // Soft ID Muon POG Run 2
        onestmuon[iDau] = (dau.isMuon() ? muon::isGoodMuon(muon, muon::SelectionType::TMOneStationTight) : false);
        npixellayer[iDau] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() : -1);
        hpmuon[iDau] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->quality(reco::TrackBase::highPurity) : false);
        muondxy[iDau] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->dxy(bestvtx) : 99.);
        muondz[iDau] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->dz(bestvtx) : 99.);
        softmuon[iDau] = (
                              onestmuon[iDau] &&
                              (ntrackerlayer[iDau] > 5) &&
                              (npixellayer[iDau] > 0) &&
                              hpmuon[iDau] &&
                              (fabs(muondxy[iDau]) < 0.3) &&
                              (fabs(muondz[iDau]) < 20.)
                              );

        // Hybrid Soft ID HIN PAG Run 2 PbPb
        trkmuon[iDau] = (dau.isMuon() ? muon.isTrackerMuon() : false);
        hybridmuon[iDau] = (
                                glbmuon[iDau] &&
                                (ntrackerlayer[iDau] > 5) &&
                                (npixellayer[iDau] > 0) &&
                                (fabs(muondxy[iDau]) < 0.3) &&
                                (fabs(muondz[iDau]) < 20.)
                                );

        // Muon Trigger Matching
        for(ushort iTr=0; iTr<filterNames_.size(); iTr++)
        {
          trgmuon[iDau][iTr] = false;
          if(dau.isMuon())
          {
            const auto& muHLTMatchesFilter = muon.triggerObjectMatchesByFilter(filterNames_.at(iTr));
            if(muHLTMatchesFilter.size()>0) trgmuon[iDau][iTr] = true;
          }
        }

        if(doMuonFull_)
        {
          nmatchedch[iDau] = (dau.isMuon() ? muon.numberOfMatches() : -1);
          matchedenergy[iDau] = (dau.isMuon() ? muon.calEnergy().hadMax : -99.);

          dx_seg[iDau] = 999.9;
          dy_seg[iDau] = 999.9;
          dxSig_seg[iDau] = 999.9;
          dySig_seg[iDau] = 999.9;
          ddxdz_seg[iDau] = 999.9;
          ddydz_seg[iDau] = 999.9;
          ddxdzSig_seg[iDau] = 999.9;
          ddydzSig_seg[iDau] = 999.9;
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

              if(dseg < std::sqrt(dx_seg[iDau]*dx_seg[iDau] + dy_seg[iDau]*dy_seg[iDau]))
              {
                dx_seg[iDau] = x_seg - x_exp;
                dy_seg[iDau] = y_seg - y_exp;
                dxSig_seg[iDau] = dx_seg[iDau] / dxerr_seg;
                dySig_seg[iDau] = dy_seg[iDau] / dyerr_seg;
                ddxdz_seg[iDau] = dxdz_seg - dxdz_exp;
                ddydz_seg[iDau] = dydz_seg - dydz_exp;
                ddxdzSig_seg[iDau] = ddxdz_seg[iDau] / ddxdzerr_seg;
                ddydzSig_seg[iDau] = ddydz_seg[iDau] / ddydzerr_seg;
              }
            }
          }
        }
      }
    }
 
    if(twoLayerDecay_)
    {
      const auto& d = *(trk.daughter(0));
      grand_mass = d.mass();
      for(ushort iGDau=0; iGDau<NGDAU_; iGDau++)
      {
        if(!d.daughter(iGDau)) continue;
        const auto& gd = *(d.daughter(iGDau));
        const auto& gdau = gd.get<reco::TrackRef>();

        //trk quality
        grand_trkquality[iGDau] = (gdau.isNonnull() ? gdau->quality(reco::TrackBase::highPurity) : false);

        //trk dEdx
        grand_H2dedx[iGDau] = -999.9;
        if(gdau.isNonnull() && dEdxHandle1.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle1.product();
          grand_H2dedx[iGDau] = dEdxTrack[gdau].dEdx();
        }   
        grand_T4dedx[iGDau] = -999.9;
        if(gdau.isNonnull() && dEdxHandle2.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle2.product();
          grand_T4dedx[iGDau] = dEdxTrack[gdau].dEdx();
        }

        //track pt
        grand_pt[iGDau] = gd.pt();

        //track momentum
        grand_p[iGDau] = gd.p();

        //track eta
        grand_eta[iGDau] = gd.eta();

        //track charge
        grand_charge[iGDau] = gd.charge();

        //track Chi2
        grand_trkChi[iGDau] = (gdau.isNonnull() ? gdau->normalizedChi2() : 99.);

        //track pT error
        grand_ptErr[iGDau] = (gdau.isNonnull() ? gdau->ptError() : -1.);

        //trkNHits
        grand_nhit[iGDau] = (gdau.isNonnull() ? gdau->numberOfValidHits() : -1);

        //DCA
        grand_dzos[iGDau] = 99.;
        grand_dxyos[iGDau] = 99.;
        if(gdau.isNonnull())
        {
          const double& gdzbest = gdau->dz(bestvtx);
          const double& gdxybest = gdau->dxy(bestvtx);
          const double& gdzerror = std::sqrt(gdau->dzError()*gdau->dzError() + bestvzError*bestvzError);
          const double& gdxyerror = std::sqrt(gdau->d0Error()*gdau->d0Error() + bestvxError*bestvyError);
          grand_dzos[iGDau] = gdzbest/gdzerror;
          grand_dxyos[iGDau] = gdxybest/gdxyerror;
        }
      }
   
      //vtxChi2
      grand_vtxChi2 = d.vertexChi2();
      grand_ndf = d.vertexNdof();
      grand_VtxProb = TMath::Prob(grand_vtxChi2, grand_ndf);

      //PAngle
      const double& secvz = d.vz(), secvx = d.vx(), secvy = d.vy();
      const TVector3 ptosvec(secvx-bestvx, secvy-bestvy, secvz-bestvz);
      const TVector3 secvec(d.px(), d.py(), d.pz());            
      const TVector3 ptosvec2D(secvx-bestvx, secvy-bestvy, 0);
      const TVector3 secvec2D(d.px(), d.py(), 0);

      grand_agl = std::cos(secvec.Angle(ptosvec));
      grand_agl_abs = secvec.Angle(ptosvec);
      grand_agl2D = std::cos(secvec2D.Angle(ptosvec2D));
      grand_agl2D_abs = secvec2D.Angle(ptosvec2D);

      //Decay length 3D
      const SMatrixSym3D& totalCov = vtx.covariance() + d.vertexCovariance();
      const SVector3 distanceVector(secvx-bestvx, secvy-bestvy, secvz-bestvz);

      grand_dl = ROOT::Math::Mag(distanceVector);
      grand_dlerror = std::sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/grand_dl;
      grand_dlos = grand_dl/grand_dlerror;

      //Decay length 2D
      const SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1), vtx.covariance(1,1), 0, 0, 0);
      const SVector6 v2(d.vertexCovariance(0,0), d.vertexCovariance(0,1), d.vertexCovariance(1,1), 0, 0, 0);
      const SMatrixSym3D totalCov2D = SMatrixSym3D(v1) + SMatrixSym3D(v2);
      const SVector3 distanceVector2D(secvx-bestvx, secvy-bestvy, 0);

      const double& gdl2D = ROOT::Math::Mag(distanceVector2D);
      const double& gdl2Derror = std::sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/gdl2D;
      grand_dlos2D = gdl2D/gdl2Derror;
    }

    if(saveTree_) PATCompositeNtuple->Fill();

    if(saveHistogram_)
    {
      for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
      {
        for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
        {
          if(pt<pTBins_[ipt+1] && pt>pTBins_[ipt] && y<yBins_[iy+1] && y>yBins_[iy])
          {
            hMassVsMVA[iy][ipt]->Fill(mva, mass);

            if(saveAllHistogram_)
            {
              hpTVsMVA[iy][ipt]->Fill(mva, pt);
              hetaVsMVA[iy][ipt]->Fill(mva, eta);
              hyVsMVA[iy][ipt]->Fill(mva, y);
              hVtxProbVsMVA[iy][ipt]->Fill(mva, VtxProb);
              h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva, agl);
              h3DPointingAngleVsMVA[iy][ipt]->Fill(mva, agl_abs);
              h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva, agl2D);
              h2DPointingAngleVsMVA[iy][ipt]->Fill(mva, agl2D_abs);
              h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva, dlos);
              h3DDecayLengthVsMVA[iy][ipt]->Fill(mva, dl);
              h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva, dlos2D);
              h2DDecayLengthVsMVA[iy][ipt]->Fill(mva, dl2D);
              for (ushort iDau=0; iDau<NDAU_; iDau++)
              {
                hzDCASignificanceDaugtherVsMVA[iDau][iy][ipt]->Fill(mva, dzos[iDau]);
                hxyDCASignificanceDaugtherVsMVA[iDau][iy][ipt]->Fill(mva, dxyos[iDau]);
                hNHitDVsMVA[iDau][iy][ipt]->Fill(mva, nhit[iDau]);
                hpTDVsMVA[iDau][iy][ipt]->Fill(mva, ptDau[iDau]);
                hpTerrDVsMVA[iDau][iy][ipt]->Fill(mva, ptErr[iDau]/ptDau[iDau]);
                hEtaDVsMVA[iDau][iy][ipt]->Fill(mva, etaDau[iDau]);
                hdedxHarmonic2DVsMVA[iDau][iy][ipt]->Fill(mva, H2dedx[iDau]);
                hdedxHarmonic2DVsP[iDau][iy][ipt]->Fill(pDau[iDau], H2dedx[iDau]);
              }
            }
          }
        }
      }
    }
  }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
PATCompositeNtupleProducer::beginJob()
{
  TH1D::SetDefaultSumw2();

  // Check inputs
  if((!threeProngDecay_ && NDAU_!=2) || (threeProngDecay_ && NDAU_!=3))
  {
    throw cms::Exception("PATCompositeAnalyzer") << "Want threeProngDecay but PID daughter vector size is: " << NDAU_ << " !" << std::endl;
  }
  if(!doRecoNtuple_) throw cms::Exception("PATCompositeAnalyzer") << "No output for RECO!! Fix config!!" << std::endl;
  if(twoLayerDecay_ && doMuon_) throw cms::Exception("PATCompositeAnalyzer") << "Muons cannot be coming from two layer decay!! Fix config!!" << std::endl;

  if(saveHistogram_) initHistogram();
  if(saveTree_) initTree();
}


void
PATCompositeNtupleProducer::initHistogram()
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
PATCompositeNtupleProducer::initTree()
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
      for(ushort iEP=1; iEP<=3; iEP++)
      {
        PATCompositeNtuple->Branch(Form("ephfpAngle%d",iEP),&ephfpAngle[iEP-1],Form("ephfpAngle%d/F",iEP));
        PATCompositeNtuple->Branch(Form("ephfmAngle%d",iEP),&ephfmAngle[iEP-1],Form("ephfmAngle%d/F",iEP));
        PATCompositeNtuple->Branch(Form("ephfpQ%d",iEP),&ephfpQ[iEP-1],Form("ephfpQ%d/F",iEP));
        PATCompositeNtuple->Branch(Form("ephfmQ%d",iEP),&ephfmQ[iEP-1],Form("ephfmQ%d/F",iEP));
      }
      PATCompositeNtuple->Branch("ephfpSumW",&ephfpSumW,"ephfpSumW/F");
      PATCompositeNtuple->Branch("ephfmSumW",&ephfmSumW,"ephfmSumW/F");
    }
    for(ushort iTRG=0; iTRG<NTRG_; iTRG++)
    {
      PATCompositeNtuple->Branch(Form("trigPrescale%d",iTRG),&trigPrescale[iTRG],Form("trigPrescale%d/S",iTRG));
      PATCompositeNtuple->Branch(Form("trigHLT%d",iTRG),&trigHLT[iTRG],Form("trigHLT%d/O",iTRG));
    }
    for(ushort iSEL=0; iSEL<NSEL_; iSEL++)
    {
      PATCompositeNtuple->Branch(Form("evtSel%d",iSEL),&evtSel[iSEL],Form("evtSel%d/O",iSEL));
    }

    // particle info
    PATCompositeNtuple->Branch("pT",&pt,"pT/F");
    PATCompositeNtuple->Branch("eta",&eta,"eta/F");
    PATCompositeNtuple->Branch("phi",&phi,"phi/F");
    PATCompositeNtuple->Branch("mass",&mass,"mass/F");
    PATCompositeNtuple->Branch("y",&y,"y/F");
    if(useAnyMVA_) PATCompositeNtuple->Branch("mva",&mva,"mva/F");

    if(!isSkimMVA_)
    {
      //Composite candidate info RECO
      PATCompositeNtuple->Branch("flavor",&flavor,"flavor/F");
      PATCompositeNtuple->Branch("VtxProb",&VtxProb,"VtxProb/F");
      PATCompositeNtuple->Branch("3DCosPointingAngle",&agl,"3DCosPointingAngle/F");
      PATCompositeNtuple->Branch("3DPointingAngle",&agl_abs,"3DPointingAngle/F");
      PATCompositeNtuple->Branch("2DCosPointingAngle",&agl2D,"2DCosPointingAngle/F");
      PATCompositeNtuple->Branch("2DPointingAngle",&agl2D_abs,"2DPointingAngle/F");
      PATCompositeNtuple->Branch("3DDecayLengthSignificance",&dlos,"3DDecayLengthSignificance/F");
      PATCompositeNtuple->Branch("3DDecayLength",&dl,"3DDecayLength/F");
      PATCompositeNtuple->Branch("3DDecayLengthError",&dlerror,"3DDecayLengthError/F");
      PATCompositeNtuple->Branch("2DDecayLengthSignificance",&dlos2D,"2DDecayLengthSignificance/F");
      PATCompositeNtuple->Branch("2DDecayLength",&dl2D,"2DDecayLength/F");

      if(doGenMatching_)
      {
        PATCompositeNtuple->Branch("isSwap",&isSwap,"isSwap/O");
        PATCompositeNtuple->Branch("idmom_reco",&idmom_reco,"idmom_reco/I");
        PATCompositeNtuple->Branch("matchGEN",&matchGEN,"matchGEN/O");
      }
 
      if(doGenMatchingTOF_)
      {
        for(ushort iDau=1; iDau<=NDAU_; iDau++)
        {
          PATCompositeNtuple->Branch(Form("PIDD%d",iDau),&pid[iDau-1],Form("PIDD%d/I",iDau));
        }
      }

      //daugther & grand daugther info
      if(twoLayerDecay_)
      {
        PATCompositeNtuple->Branch("massDaugther1",&grand_mass,"massDaugther1/F");
        PATCompositeNtuple->Branch("VtxProbDaugther1",&grand_VtxProb,"VtxProbDaugther1/F");
        PATCompositeNtuple->Branch("3DCosPointingAngleDaugther1",&grand_agl,"3DCosPointingAngleDaugther1/F");
        PATCompositeNtuple->Branch("3DPointingAngleDaugther1",&grand_agl_abs,"3DPointingAngleDaugther1/F");
        PATCompositeNtuple->Branch("2DCosPointingAngleDaugther1",&grand_agl2D,"2DCosPointingAngleDaugther1/F");
        PATCompositeNtuple->Branch("2DPointingAngleDaugther1",&grand_agl2D_abs,"2DPointingAngleDaugther1/F");
        PATCompositeNtuple->Branch("3DDecayLengthSignificanceDaugther1",&grand_dlos,"3DDecayLengthSignificanceDaugther1/F");
        PATCompositeNtuple->Branch("3DDecayLengthDaugther1",&grand_dl,"3DDecayLengthDaugther1/F");
        PATCompositeNtuple->Branch("3DDecayLengthErrorDaugther1",&grand_dlerror,"3DDecayLengthErrorDaugther1/F");
        PATCompositeNtuple->Branch("2DDecayLengthSignificanceDaugther1",&grand_dlos2D,"2DDecayLengthSignificanceDaugther1/F");
        for(ushort iGDau=1; iGDau<=NGDAU_; iGDau++)
        {
          PATCompositeNtuple->Branch(Form("zDCASignificanceGrandDaugther%d",iGDau),&grand_dzos[iGDau-1],Form("zDCASignificanceGrandDaugther%d/F",iGDau));
          PATCompositeNtuple->Branch(Form("xyDCASignificanceGrandDaugther%d",iGDau),&grand_dxyos[iGDau-1],Form("xyDCASignificanceGrandDaugther%d/F",iGDau));
          PATCompositeNtuple->Branch(Form("NHitGrandD%d",iGDau),&grand_nhit[iGDau-1],Form("NHitGrandD%d/F",iGDau));
          PATCompositeNtuple->Branch(Form("HighPurityGrandDaugther%d",iGDau),&grand_trkquality[iGDau-1],Form("HighPurityGrandDaugther%d/O",iGDau));
          PATCompositeNtuple->Branch(Form("pTGrandD%d",iGDau),&grand_pt[iGDau-1],Form("pTGrandD%d/F",iGDau));
          PATCompositeNtuple->Branch(Form("pTerrGrandD%d",iGDau),&grand_ptErr[iGDau-1],Form("pTerrGrandD%d/F",iGDau));
          PATCompositeNtuple->Branch(Form("EtaGrandD%d",iGDau),&grand_eta[iGDau-1],Form("EtaGrandD%d/F",iGDau));
          if(useDeDxData_)
          {
            PATCompositeNtuple->Branch(Form("dedxPixelHarmonic2GrandD%d",iGDau),&grand_T4dedx[iGDau-1],Form("dedxPixelHarmonic2GrandD%d/F",iGDau));
            PATCompositeNtuple->Branch(Form("dedxHarmonic2GrandD%d",iGDau),&grand_H2dedx[iGDau-1],Form("dedxHarmonic2GrandD%d/F",iGDau));
          }
        }
      }
      for(ushort iDau=1; iDau<=NDAU_; iDau++)
      {
        PATCompositeNtuple->Branch(Form("zDCASignificanceDaugther%d",iDau),&dzos[iDau-1],Form("zDCASignificanceDaugther%d/F",iDau));
        PATCompositeNtuple->Branch(Form("xyDCASignificanceDaugther%d",iDau),&dxyos[iDau-1],Form("xyDCASignificanceDaugther%d/F",iDau));
        PATCompositeNtuple->Branch(Form("NHitD%d",iDau),&nhit[iDau-1],Form("NHitD%d/F",iDau));
        PATCompositeNtuple->Branch(Form("HighPurityDaugther%d",iDau),&trkquality[iDau-1],Form("HighPurityDaugther%d/O",iDau));
        PATCompositeNtuple->Branch(Form("pTD%d",iDau),&ptDau[iDau-1],Form("pTD%d/F",iDau));
        PATCompositeNtuple->Branch(Form("pTerrD%d",iDau),&ptErr[iDau-1],Form("pTerrD%d/F",iDau));
        PATCompositeNtuple->Branch(Form("EtaD%d",iDau),&etaDau[iDau-1],Form("EtaD%d/F",iDau));
        PATCompositeNtuple->Branch(Form("PhiD%d",iDau),&phiDau[iDau-1],Form("PhiD%d/F",iDau));
        PATCompositeNtuple->Branch(Form("chargeD%d",iDau),&chargeDau[iDau-1],Form("chargeD%d/S",iDau));
        if(useDeDxData_)
        {
          PATCompositeNtuple->Branch(Form("dedxPixelHarmonic2D%d",iDau),&T4dedx[iDau-1],Form("dedxPixelHarmonic2D%d/F",iDau));
          PATCompositeNtuple->Branch(Form("dedxHarmonic2D%d",iDau),&H2dedx[iDau-1],Form("dedxHarmonic2D%d/F",iDau));
        }
      }
 
      if(doMuon_)
      {
        for(ushort iDau=1; iDau<=NDAU_; iDau++)
        {
          if(fabs(PID_dau_[iDau-1])!=13) continue;
          PATCompositeNtuple->Branch(Form("OneStMuon%d",iDau),&onestmuon[iDau-1],Form("OneStMuon%d/O",iDau));
          PATCompositeNtuple->Branch(Form("PFMuon%d",iDau),&pfmuon[iDau-1],Form("PFMuon%d/O",iDau));
          PATCompositeNtuple->Branch(Form("GlbMuon%d",iDau),&glbmuon[iDau-1],Form("GlbMuon%d/O",iDau));
          PATCompositeNtuple->Branch(Form("trkMuon%d",iDau),&trkmuon[iDau-1],Form("trkMuon%d/O",iDau));
          PATCompositeNtuple->Branch(Form("tightMuon%d",iDau),&tightmuon[iDau-1],Form("tightMuon%d/O",iDau));
          PATCompositeNtuple->Branch(Form("softMuon%d",iDau),&softmuon[iDau-1],Form("softMuon%d/O",iDau));
          PATCompositeNtuple->Branch(Form("hybridMuon%d",iDau),&hybridmuon[iDau-1],Form("hybridMuon%d/O",iDau));
          PATCompositeNtuple->Branch(Form("HPMuon%d",iDau),&hpmuon[iDau-1],Form("hybridMuon%d/O",iDau));
          for(ushort iTRG=0; iTRG<NTRG_; iTRG++)
          {
            PATCompositeNtuple->Branch(Form("trigMuon%d_%d",iDau,iTRG),&trgmuon[iDau-1][iTRG],Form("trigMuon%d_%d/O",iDau,iTRG));
          }
          PATCompositeNtuple->Branch(Form("nMatchedStationD%d",iDau),&nmatchedst[iDau-1],Form("nMatchedStationD%d/S",iDau));
          PATCompositeNtuple->Branch(Form("nTrackerLayerD%d",iDau),&ntrackerlayer[iDau-1],Form("nTrackerLayerD%d/S",iDau));
          PATCompositeNtuple->Branch(Form("nPixelLayerD%d",iDau),&npixellayer[iDau-1],Form("nPixelLayerD%d/S",iDau));
          PATCompositeNtuple->Branch(Form("nPixelHitD%d",iDau),&npixelhit[iDau-1],Form("nPixelHitD%d/S",iDau));
          PATCompositeNtuple->Branch(Form("nMuonHitD%d",iDau),&nmuonhit[iDau-1],Form("nMuonHitD%d/S",iDau));
          PATCompositeNtuple->Branch(Form("GlbTrkChiD%d",iDau),&glbtrkchi[iDau-1],Form("GlbTrkChiD%d/F",iDau));
          PATCompositeNtuple->Branch(Form("muondXYD%d",iDau),&muonbestdxy[iDau-1],Form("muondXYD%d/F",iDau));
          PATCompositeNtuple->Branch(Form("muondZD%d",iDau),&muonbestdz[iDau-1],Form("muondZD%d/F",iDau));
          PATCompositeNtuple->Branch(Form("dXYD%d",iDau),&muondxy[iDau-1],Form("dXYD%d/F",iDau));
          PATCompositeNtuple->Branch(Form("dZD%d",iDau),&muondz[iDau-1],Form("dZD%d/F",iDau));
          if(doMuonFull_)
          {
            PATCompositeNtuple->Branch(Form("nMatchedChamberD%d",iDau),&nmatchedch[iDau-1],Form("nMatchedChamberD%d/S",iDau));
            PATCompositeNtuple->Branch(Form("EnergyDepositionD%d",iDau),&matchedenergy[iDau-1],Form("EnergyDepositionD%d/F",iDau));
            PATCompositeNtuple->Branch(Form("dx%d_seg",iDau),        &dx_seg[iDau-1], Form("dx%d_seg/F",iDau));
            PATCompositeNtuple->Branch(Form("dy%d_seg",iDau),        &dy_seg[iDau-1], Form("dy%d_seg/F",iDau));
            PATCompositeNtuple->Branch(Form("dxSig%d_seg",iDau),     &dxSig_seg[iDau-1], Form("dxSig%d_seg/F",iDau));
            PATCompositeNtuple->Branch(Form("dySig%d_seg",iDau),     &dySig_seg[iDau-1], Form("dySig%d_seg/F",iDau));
            PATCompositeNtuple->Branch(Form("ddxdz%d_seg",iDau),     &ddxdz_seg[iDau-1], Form("ddxdz%d_seg/F",iDau));
            PATCompositeNtuple->Branch(Form("ddydz%d_seg",iDau),     &ddydz_seg[iDau-1], Form("ddydz%d_seg/F",iDau));
            PATCompositeNtuple->Branch(Form("ddxdzSig%d_seg",iDau),  &ddxdzSig_seg[iDau-1], Form("ddxdzSig%d_seg/F",iDau));
            PATCompositeNtuple->Branch(Form("ddydzSig%d_seg",iDau),  &ddydzSig_seg[iDau-1], Form("ddydzSig%d_seg/F",iDau));
          }
        }
      }
    }
  } // doRecoNtuple_
}


//--------------------------------------------------------------------------------------------------
void 
PATCompositeNtupleProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = true;
  EDConsumerBase::Labels triggerResultsLabel;
  EDConsumerBase::labelsForToken(tok_triggerResults_, triggerResultsLabel);
  hltPrescaleProvider_.init(iRun, iSetup, triggerResultsLabel.process, changed);
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
PATCompositeNtupleProducer::endJob()
{
}


reco::GenParticleRef
PATCompositeNtupleProducer::findMother(const reco::GenParticleRef& genParRef)
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


//define this as a plug-in
DEFINE_FWK_MODULE(PATCompositeNtupleProducer);
