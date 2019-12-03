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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
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
#define MAXTRG 1024
#define MAXSEL 100

//
// class decleration
//

class EventInfoTreeProducer : public edm::EDAnalyzer {
public:
  explicit EventInfoTreeProducer(const edm::ParameterSet&);
  ~EventInfoTreeProducer();

private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void initTree();

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree* EventInfoNtuple;

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

  bool isCentrality_;
  bool isEventPlane_;

  //token
  edm::EDGetTokenT<reco::BeamSpot> tok_offlineBS_;
  edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;

  //trigger
  const std::vector<std::string> triggerNames_;
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

EventInfoTreeProducer::EventInfoTreeProducer(const edm::ParameterSet& iConfig) :
  triggerNames_(iConfig.getUntrackedParameter<std::vector<std::string> >("triggerPathNames")),
  NTRG_(triggerNames_.size()>MAXTRG ? MAXTRG : triggerNames_.size()),
  eventFilters_(iConfig.getUntrackedParameter<std::vector<std::string> >("eventFilterNames")),
  NSEL_(eventFilters_.size()>MAXSEL ? MAXSEL : eventFilters_.size()),
  selectEvents_(iConfig.getUntrackedParameter<std::string>("selectEvents")),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
  //input tokens
  tok_offlineBS_ = consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"));
  tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));

  isCentrality_ = (iConfig.exists("isCentrality") ? iConfig.getParameter<bool>("isCentrality") : false);
  if(isCentrality_)
  {
    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

  isEventPlane_ = (iConfig.exists("isEventPlane") ? iConfig.getParameter<bool>("isEventPlane") : false);
  if(isEventPlane_) tok_eventplaneSrc_ = consumes<reco::EvtPlaneCollection>(iConfig.getParameter<edm::InputTag>("eventplaneSrc"));

  tok_triggerResults_ = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultCollection"));
  tok_filterResults_ = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("FilterResultCollection"));
}


EventInfoTreeProducer::~EventInfoTreeProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EventInfoTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  fillRECO(iEvent, iSetup);
  EventInfoNtuple->Fill();
}


void
EventInfoTreeProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get collection
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(tok_offlineBS_, beamspot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);
  if(!vertices.isValid()) throw cms::Exception("EventInfoTreeProducer") << "Primary vertices  collection not found!" << std::endl;

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
}

// ------------ method called once each job just before starting event
//loop  ------------
void
EventInfoTreeProducer::beginJob()
{
  TH1D::SetDefaultSumw2();

  initTree();
}

void 
EventInfoTreeProducer::initTree()
{ 
  EventInfoNtuple = fs->make< TTree>("EventInfoNtuple","EventInfoNtuple");

  // Event info
  EventInfoNtuple->Branch("RunNb",&runNb,"RunNb/i");
  EventInfoNtuple->Branch("LSNb",&lsNb,"LSNb/i");
  EventInfoNtuple->Branch("EventNb",&eventNb,"EventNb/i");
  EventInfoNtuple->Branch("nPV",&nPV,"nPV/S");
  EventInfoNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
  EventInfoNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
  EventInfoNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");
  if(isCentrality_) 
  {
    EventInfoNtuple->Branch("centrality",&centrality,"centrality/S");
    EventInfoNtuple->Branch("Npixel",&Npixel,"Npixel/I");
    EventInfoNtuple->Branch("HFsumETPlus",&HFsumETPlus,"HFsumETPlus/F");
    EventInfoNtuple->Branch("HFsumETMinus",&HFsumETMinus,"HFsumETMinus/F");
    EventInfoNtuple->Branch("ZDCPlus",&ZDCPlus,"ZDCPlus/F");
    EventInfoNtuple->Branch("ZDCMinus",&ZDCMinus,"ZDCMinus/F");
    EventInfoNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
  }
  if(isEventPlane_) {
    EventInfoNtuple->Branch("ephfpAngle",ephfpAngle,"ephfpAngle[3]/F");
    EventInfoNtuple->Branch("ephfmAngle",ephfmAngle,"ephfmAngle[3]/F");
    EventInfoNtuple->Branch("ephfpQ",ephfpQ,"ephfpQ[3]/F");
    EventInfoNtuple->Branch("ephfmQ",ephfmQ,"ephfmQ[3]/F");
    EventInfoNtuple->Branch("ephfpSumW",&ephfpSumW,"ephfpSumW/F");
    EventInfoNtuple->Branch("ephfmSumW",&ephfmSumW,"ephfmSumW/F");
  }
  EventInfoNtuple->Branch("trigPrescale",trigPrescale,Form("trigPrescale[%d]/S",NTRG_));
  EventInfoNtuple->Branch("trigHLT",trigHLT,Form("trigHLT[%d]/O",NTRG_));
  EventInfoNtuple->Branch("evtSel",evtSel,Form("evtSel[%d]/O",NSEL_));
}


//--------------------------------------------------------------------------------------------------
void 
EventInfoTreeProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = true;
  EDConsumerBase::Labels triggerResultsLabel;
  EDConsumerBase::labelsForToken(tok_triggerResults_, triggerResultsLabel);
  hltPrescaleProvider_.init(iRun, iSetup, triggerResultsLabel.process, changed);
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
EventInfoTreeProducer::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventInfoTreeProducer);
