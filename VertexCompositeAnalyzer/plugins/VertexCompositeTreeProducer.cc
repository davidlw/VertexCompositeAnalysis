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
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TMath.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define PI 3.1416
#define MAXCAN 50000

using namespace std;

class VertexCompositeTreeProducer : public edm::EDAnalyzer {
public:
  explicit VertexCompositeTreeProducer(const edm::ParameterSet&);
  ~VertexCompositeTreeProducer();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&) ;
  virtual void fillGEN(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  virtual void initHistogram();
  virtual void initTree();

  int muAssocToTrack( const reco::TrackRef& trackref, const edm::Handle<reco::MuonCollection>& muonh) const;

  // ----------member data ---------------------------
    
    edm::Service<TFileService> fs;

    TTree* VertexCompositeNtuple;
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
    TH2F*  hzDCASignificanceDaugther1VsMVA[6][10];
    TH2F*  hxyDCASignificanceDaugther1VsMVA[6][10];
    TH2F*  hNHitD1VsMVA[6][10];
    TH2F*  hpTD1VsMVA[6][10];
    TH2F*  hpTerrD1VsMVA[6][10];
    TH2F*  hEtaD1VsMVA[6][10];
    TH2F*  hdedxHarmonic2D1VsMVA[6][10];
    TH2F*  hdedxHarmonic2D1VsP[6][10];
    TH2F*  hzDCASignificanceDaugther2VsMVA[6][10];
    TH2F*  hxyDCASignificanceDaugther2VsMVA[6][10];
    TH2F*  hNHitD2VsMVA[6][10];
    TH2F*  hpTD2VsMVA[6][10];
    TH2F*  hpTerrD2VsMVA[6][10];
    TH2F*  hEtaD2VsMVA[6][10];
    TH2F*  hdedxHarmonic2D2VsMVA[6][10];
    TH2F*  hdedxHarmonic2D2VsP[6][10];
    TH2F*  hzDCASignificanceDaugther3VsMVA[6][10];
    TH2F*  hxyDCASignificanceDaugther3VsMVA[6][10];
    TH2F*  hNHitD3VsMVA[6][10];
    TH2F*  hpTD3VsMVA[6][10];
    TH2F*  hpTerrD3VsMVA[6][10];
    TH2F*  hEtaD3VsMVA[6][10];
    TH2F*  hdedxHarmonic2D3VsMVA[6][10];
    TH2F*  hdedxHarmonic2D3VsP[6][10];
    
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
    bool hasSwap_;
    bool decayInGen_;
    bool twoLayerDecay_;
    bool threeProngDecay_;
    bool doMuon_;
    bool doMuonFull_;
    int PID_;
    int PID_dau1_;
    int PID_dau2_;
    int PID_dau3_;
    
    //cut variables
    double multMax_;
    double multMin_;
    double deltaR_; //deltaR for Gen matching

    vector<double> pTBins_;
    vector<double> yBins_;

    //tree branches
    //event info
    int centrality;
    int Ntrkoffline;
    int Npixel;
    float HFsumET;
    float bestvx;
    float bestvy;
    float bestvz;
    int candSize;
    
    //Composite candidate info
    float mva[MAXCAN];
    float pt[MAXCAN];
    float eta[MAXCAN];
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
    float dzos1[MAXCAN];
    float dzos2[MAXCAN];
    float dzos3[MAXCAN];
    float dxyos1[MAXCAN];
    float dxyos2[MAXCAN];
    float dxyos3[MAXCAN];
    float nhit1[MAXCAN];
    float nhit2[MAXCAN];
    float nhit3[MAXCAN];
    bool trkquality1[MAXCAN];
    bool trkquality2[MAXCAN];
    bool trkquality3[MAXCAN];
    float pt1[MAXCAN];
    float pt2[MAXCAN];
    float pt3[MAXCAN];
    float ptErr1[MAXCAN];
    float ptErr2[MAXCAN];
    float ptErr3[MAXCAN];
    float p1[MAXCAN];
    float p2[MAXCAN];
    float p3[MAXCAN];
    float eta1[MAXCAN];
    float eta2[MAXCAN];
    float eta3[MAXCAN];
    float phi1[MAXCAN];
    float phi2[MAXCAN];
    float phi3[MAXCAN];
    int charge1[MAXCAN];
    int charge2[MAXCAN];
    int charge3[MAXCAN];
    int pid1[MAXCAN];
    int pid2[MAXCAN];
    int pid3[MAXCAN];
    float tof1[MAXCAN];
    float tof2[MAXCAN];
    float tof3[MAXCAN];
    float H2dedx1[MAXCAN];
    float H2dedx2[MAXCAN];
    float H2dedx3[MAXCAN];
    float T4dedx1[MAXCAN];
    float T4dedx2[MAXCAN];
    float T4dedx3[MAXCAN];
    float trkChi1[MAXCAN];
    float trkChi2[MAXCAN];
    float trkChi3[MAXCAN];
   
    //grand-dau info
    float grand_dzos1[MAXCAN];
    float grand_dzos2[MAXCAN];
    float grand_dxyos1[MAXCAN];
    float grand_dxyos2[MAXCAN];
    float grand_nhit1[MAXCAN];
    float grand_nhit2[MAXCAN];
    bool grand_trkquality1[MAXCAN];
    bool grand_trkquality2[MAXCAN];
    float grand_pt1[MAXCAN];
    float grand_pt2[MAXCAN];
    float grand_ptErr1[MAXCAN];
    float grand_ptErr2[MAXCAN];
    float grand_p1[MAXCAN];
    float grand_p2[MAXCAN];
    float grand_eta1[MAXCAN];
    float grand_eta2[MAXCAN];
    int grand_charge1[MAXCAN];
    int grand_charge2[MAXCAN];
    float grand_H2dedx1[MAXCAN];
    float grand_H2dedx2[MAXCAN];
    float grand_T4dedx1[MAXCAN];
    float grand_T4dedx2[MAXCAN];
    float grand_trkChi1[MAXCAN];
    float grand_trkChi2[MAXCAN];
    
    //dau muon info
    bool  onestmuon1[MAXCAN];
    bool  onestmuon2[MAXCAN];
    bool  pfmuon1[MAXCAN];
    bool  pfmuon2[MAXCAN];
    bool  glbmuon1[MAXCAN];
    bool  glbmuon2[MAXCAN];
    bool  trkmuon1[MAXCAN];
    bool  trkmuon2[MAXCAN];
    bool  calomuon1[MAXCAN];
    bool  calomuon2[MAXCAN];
    float nmatchedst1[MAXCAN];
    float nmatchedch1[MAXCAN];
    float ntrackerlayer1[MAXCAN];
    float npixellayer1[MAXCAN];
    float matchedenergy1[MAXCAN];
    float nmatchedst2[MAXCAN];
    float nmatchedch2[MAXCAN];
    float ntrackerlayer2[MAXCAN];
    float npixellayer2[MAXCAN];
    float matchedenergy2[MAXCAN];
    float dx1_seg_[MAXCAN];
    float dy1_seg_[MAXCAN];
    float dxSig1_seg_[MAXCAN];
    float dySig1_seg_[MAXCAN];
    float ddxdz1_seg_[MAXCAN];
    float ddydz1_seg_[MAXCAN];
    float ddxdzSig1_seg_[MAXCAN];
    float ddydzSig1_seg_[MAXCAN];
    float dx2_seg_[MAXCAN];
    float dy2_seg_[MAXCAN];
    float dxSig2_seg_[MAXCAN];
    float dySig2_seg_[MAXCAN];
    float ddxdz2_seg_[MAXCAN];
    float ddydz2_seg_[MAXCAN];
    float ddxdzSig2_seg_[MAXCAN];
    float ddydzSig2_seg_[MAXCAN];

    // gen info    
    int candSize_gen;
    float pt_gen[MAXCAN];
    float eta_gen[MAXCAN];
    int status_gen[MAXCAN];
    int idmom[MAXCAN];
    float y_gen[MAXCAN];
    int iddau1[MAXCAN];
    int iddau2[MAXCAN];
    int iddau3[MAXCAN];

    //vector for gen match
    vector< vector<double> > *pVect;
    vector<double> *Dvector1;
    vector<double> *Dvector2;
    vector<double> *Dvector3;
    vector<int> *pVectIDmom;
    
    bool useAnyMVA_;
    bool isSkimMVA_;
    bool isCentrality_;

    edm::Handle<int> cbin_;

    //tokens
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    //edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
    edm::EDGetTokenT<pat::CompositeCandidateCollection> patCompositeCandidateCollection_Token_;
    edm::EDGetTokenT<MVACollection> MVAValues_Token_;

    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
    edm::EDGetTokenT<reco::MuonCollection> tok_muon_;

    edm::EDGetTokenT<int> tok_centBinLabel_;
    edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

VertexCompositeTreeProducer::VertexCompositeTreeProducer(const edm::ParameterSet& iConfig)
{
    //options
    doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
    doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
    twoLayerDecay_ = iConfig.getUntrackedParameter<bool>("twoLayerDecay");
    threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
    doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
    doGenMatchingTOF_ = iConfig.getUntrackedParameter<bool>("doGenMatchingTOF");
    hasSwap_ = iConfig.getUntrackedParameter<bool>("hasSwap");
    decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
    doMuon_ = iConfig.getUntrackedParameter<bool>("doMuon");
    doMuonFull_ = iConfig.getUntrackedParameter<bool>("doMuonFull");
    PID_ = iConfig.getUntrackedParameter<int>("PID");
    PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
    PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
    if(threeProngDecay_) PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");
    
    saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
    saveHistogram_ = iConfig.getUntrackedParameter<bool>("saveHistogram");
    saveAllHistogram_ = iConfig.getUntrackedParameter<bool>("saveAllHistogram");
    massHistPeak_ = iConfig.getUntrackedParameter<double>("massHistPeak");
    massHistWidth_ = iConfig.getUntrackedParameter<double>("massHistWidth");
    massHistBins_ = iConfig.getUntrackedParameter<int>("massHistBins");

    useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
    isSkimMVA_ = iConfig.getUntrackedParameter<bool>("isSkimMVA"); 

    //cut variables
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
    deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);

    pTBins_ = iConfig.getUntrackedParameter< std::vector<double> >("pTBins");
    yBins_  = iConfig.getUntrackedParameter< std::vector<double> >("yBins");

    //input tokens
    tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
    //recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"));
    patCompositeCandidateCollection_Token_ = consumes<pat::CompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"));
    MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
    tok_muon_ = consumes<reco::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection"));
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
    tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

    isCentrality_ = false;
    if(iConfig.exists("isCentrality")) isCentrality_ = iConfig.getParameter<bool>("isCentrality");
    if(isCentrality_)
    {
      tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
      tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
    }

    if(useAnyMVA_ && iConfig.exists("MVACollection"))
      MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
}


VertexCompositeTreeProducer::~VertexCompositeTreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VertexCompositeTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup&
iSetup)
{
    using std::vector;
    using namespace edm;
    using namespace reco;

    if(doGenNtuple_) fillGEN(iEvent,iSetup);
    if(doRecoNtuple_) fillRECO(iEvent,iSetup);

    if(saveTree_) VertexCompositeNtuple->Fill();
}

void
VertexCompositeTreeProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //get collections
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);

    //edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
    edm::Handle<pat::CompositeCandidateCollection> v0candidates;
    //iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,v0candidates);
    //const reco::VertexCompositeCandidateCollection * v0candidates_ = v0candidates.product();
    iEvent.getByToken(patCompositeCandidateCollection_Token_,v0candidates);
    const pat::CompositeCandidateCollection * v0candidates_ = v0candidates.product();
    
    edm::Handle<MVACollection> mvavalues;
    if(useAnyMVA_)
    {
      iEvent.getByToken(MVAValues_Token_,mvavalues);
      assert( (*mvavalues).size() == v0candidates->size() );
    }

    edm::Handle<reco::GenParticleCollection> genpars;
    if(doGenMatching_ || doGenMatchingTOF_) iEvent.getByToken(tok_genParticle_,genpars);

    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
    iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
    
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
    iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
    
    centrality=-1;
    if(isCentrality_)
    {
      edm::Handle<reco::Centrality> cent;
      iEvent.getByToken(tok_centSrc_, cent);

      iEvent.getByToken(tok_centBinLabel_,cbin_);
      centrality = *cbin_;  

      HFsumET = cent->EtHFtowerSum();
      Npixel = cent->multiplicityPixel();
//      int ntrk = cent->Ntracks();
    }
    //best vertex
    bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    //Ntrkoffline
    Ntrkoffline = 0;
    if(multMax_!=-1 && multMin_!=-1)
    {
      for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();
        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        Ntrkoffline++;
      }
    }

    //Gen info for matching
    if(doGenMatching_)
    {
        pVect = new vector< vector<double>>;
        pVectIDmom = new vector<int>;
        
//        edm::Handle<reco::GenParticleCollection> genpars;
//        iEvent.getByToken(tok_genParticle_,genpars);
        
        if(!genpars.isValid())
        {
            cout<<"Gen matching cannot be done without Gen collection!!"<<endl;
            return;
        }

        for(unsigned it=0; it<genpars->size(); ++it){
            
            const reco::GenParticle & trk = (*genpars)[it];
            
            int id = trk.pdgId();
            if(fabs(id)!=PID_) continue; //check is target
            if(decayInGen_ && trk.numberOfDaughters()!=2 && !threeProngDecay_) continue; //check 2-pron decay if target decays in Gen
            if(decayInGen_ && trk.numberOfDaughters()!=3 && threeProngDecay_) continue; //check 2-pron decay if target decays in Gen
            
            int idmom_tmp = -77;
            
            if(trk.numberOfMothers()!=0)
            {
                const reco::Candidate * mom = trk.mother();
                idmom_tmp = mom->pdgId();
            }
            
            const reco::Candidate * Dd1 = trk.daughter(0);
            const reco::Candidate * Dd2 = trk.daughter(1);
            const reco::Candidate * Dd3 = 0;            

            if(!threeProngDecay_ && !(fabs(Dd1->pdgId())==PID_dau1_ && fabs(Dd2->pdgId())==PID_dau2_) && !(fabs(Dd2->pdgId())==PID_dau1_ && fabs(Dd1->pdgId())==PID_dau2_)) continue; //check daughter id                

            if(threeProngDecay_)
            {
              Dd3 = trk.daughter(2);
              if(!(fabs(Dd1->pdgId())==PID_dau1_ && fabs(Dd2->pdgId())==PID_dau2_ && fabs(Dd3->pdgId())==PID_dau3_)
              && !(fabs(Dd1->pdgId())==PID_dau1_ && fabs(Dd2->pdgId())==PID_dau3_ && fabs(Dd3->pdgId())==PID_dau2_)
              && !(fabs(Dd1->pdgId())==PID_dau2_ && fabs(Dd2->pdgId())==PID_dau1_ && fabs(Dd3->pdgId())==PID_dau3_)
              && !(fabs(Dd1->pdgId())==PID_dau2_ && fabs(Dd2->pdgId())==PID_dau3_ && fabs(Dd3->pdgId())==PID_dau1_)
              && !(fabs(Dd1->pdgId())==PID_dau3_ && fabs(Dd2->pdgId())==PID_dau1_ && fabs(Dd3->pdgId())==PID_dau2_)
              && !(fabs(Dd1->pdgId())==PID_dau3_ && fabs(Dd2->pdgId())==PID_dau2_ && fabs(Dd3->pdgId())==PID_dau1_) ) continue;
            }

            Dvector1 = new vector<double>;
            Dvector2 = new vector<double>;
            
            Dvector1->push_back(Dd1->pt());
            Dvector1->push_back(Dd1->eta());
            Dvector1->push_back(Dd1->phi());
            Dvector1->push_back(Dd1->charge());
            Dvector1->push_back(Dd1->mass());
            
            Dvector2->push_back(Dd2->pt());
            Dvector2->push_back(Dd2->eta());
            Dvector2->push_back(Dd2->phi());
            Dvector2->push_back(Dd2->charge());
            Dvector2->push_back(Dd2->mass());
            
            pVect->push_back(*Dvector1);
            pVect->push_back(*Dvector2);
            
            pVectIDmom->push_back(idmom_tmp);
            
            delete Dvector1;
            delete Dvector2;

            if(threeProngDecay_)
            {
              Dvector3 = new vector<double>;

              Dvector3->push_back(Dd3->pt());
              Dvector3->push_back(Dd3->eta());
              Dvector3->push_back(Dd3->phi());
              Dvector3->push_back(Dd3->charge());
              Dvector3->push_back(Dd3->mass());

              pVect->push_back(*Dvector3);
              delete Dvector3;
            }
        }
    }

    //RECO Candidate info
    candSize = v0candidates_->size();
    for(unsigned it=0; it<v0candidates_->size(); ++it){
        
        //const reco::VertexCompositeCandidate & trk = (*v0candidates_)[it];
        const pat::CompositeCandidate & trk = (*v0candidates_)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();

        eta[it] = trk.eta();
        y[it] = trk.rapidity();
        pt[it] = trk.pt();
        flavor[it] = trk.pdgId()/abs(trk.pdgId());

        mva[it] = 0.0;
        if(useAnyMVA_) mva[it] = (*mvavalues)[it];

        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        mass[it] = trk.mass();
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        const reco::Candidate * d3 = 0;        
        if(threeProngDecay_) d3 = trk.daughter(2);

        //Gen match
        if(doGenMatching_)
        {
            matchGEN[it] = false;
            int nGenDau = (int)pVect->size();
            isSwap[it] = false;
            idmom_reco[it] = -77;
            
            for(int i=0;i<nGenDau;i++)
            {
                vector<double> Dvector1 = (*pVect)[i]; //get GEN daugther vector
                if(d1->charge()!=Dvector1.at(3)) continue; //check match charge
                double deltaR = sqrt(pow(d1->eta()-Dvector1.at(1),2)+pow(d1->phi()-Dvector1.at(2),2));
                
                if(deltaR > deltaR_) continue; //check deltaR matching
                if(fabs((d1->pt()-Dvector1.at(0))/d1->pt()) > 0.5) continue; //check deltaPt matching
                double d1massGEN = Dvector1.at(4);
                double d1mass = d1->mass();
                double d2massGEN=0, d2mass=0;
                double d3massGEN=0, d3mass=0;

                if(nGenDau==2)
                {
                  if(i%2==0)
                  {
                    vector<double> Dvector2 = (*pVect)[i+1]; //get GEN daugther vector for track2
                    if(d2->charge()!=Dvector2.at(3)) continue; //check match charge
                    double deltaR = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                    
                    if(deltaR > deltaR_) continue; //check deltaR matching
                    if(fabs((d2->pt()-Dvector2.at(0))/d2->pt()) > 0.5) continue; //check deltaPt matching
                    d2massGEN = Dvector2.at(4);
                    d2mass = d2->mass();
                    
                    matchGEN[it] = true; //matched gen
                  }

                  if(i%2==1)
                  {
                    vector<double> Dvector2 = (*pVect)[i-1]; //get GEN daugther vector for track2
                    if(d2->charge()!=Dvector2.at(3)) continue; //check match charge
                    double deltaR = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                    
                    if(deltaR > deltaR_) continue; //check deltaR matching
                    if(fabs((d2->pt()-Dvector2.at(0))/d2->pt()) > 0.5) continue; //check deltaPt matching
                    d2massGEN = Dvector2.at(4);
                    d2mass = d2->mass();
                    
                    matchGEN[it] = true; //matched gen
                  }

                  //check swap
                  if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01) isSwap[it] = true;

                  //check prompt & record mom id
                  idmom_reco[it] = pVectIDmom->at(i/2);
                }

                if(nGenDau==3)
                {
                  if(i%3==0)
                  {
                    vector<double> Dvector2 = (*pVect)[i+1]; //get GEN daugther vector for track2
                    vector<double> Dvector3 = (*pVect)[i+2]; //get GEN daugther vector for track3

                    if(!(d2->charge()==Dvector2.at(3) && d3->charge()==Dvector3.at(3)) 
                    && !(d3->charge()==Dvector2.at(3) && d2->charge()==Dvector3.at(3))) continue; //check match charge

                    double deltaR22 = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                    double deltaR33 = sqrt(pow(d3->eta()-Dvector3.at(1),2)+pow(d3->phi()-Dvector3.at(2),2));
                    double deltaR23 = sqrt(pow(d2->eta()-Dvector3.at(1),2)+pow(d2->phi()-Dvector3.at(2),2));
                    double deltaR32 = sqrt(pow(d3->eta()-Dvector2.at(1),2)+pow(d3->phi()-Dvector2.at(2),2));

                    if(!(deltaR22 < deltaR_ && deltaR33 < deltaR_) && !(deltaR23 < deltaR_ && deltaR32 < deltaR_) ) continue;

                    double deltaPt22 = fabs((d2->pt()-Dvector2.at(0))/d2->pt());
                    double deltaPt33 = fabs((d3->pt()-Dvector3.at(0))/d3->pt());
                    double deltaPt23 = fabs((d2->pt()-Dvector3.at(0))/d2->pt());
                    double deltaPt32 = fabs((d3->pt()-Dvector2.at(0))/d3->pt());

                    if( !(deltaPt22 < 0.5 && deltaPt33 < 0.5) && !(deltaPt23 < 0.5 && deltaPt32 < 0.5) ) continue; //check deltaPt matching

                    d2massGEN = Dvector2.at(4);
                    d2mass = d2->mass();
                    d3massGEN = Dvector3.at(4);
                    d3mass = d3->mass();

                    matchGEN[it] = true; //matched gen
                  }

                  if(i%3==1)
                  {
                    vector<double> Dvector2 = (*pVect)[i-1]; //get GEN daugther vector for track2
                    vector<double> Dvector3 = (*pVect)[i+1]; //get GEN daugther vector for track3

                    if(!(d2->charge()==Dvector2.at(3) && d3->charge()==Dvector3.at(3))
                    && !(d3->charge()==Dvector2.at(3) && d2->charge()==Dvector3.at(3))) continue; //check match charge

                    double deltaR22 = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                    double deltaR33 = sqrt(pow(d3->eta()-Dvector3.at(1),2)+pow(d3->phi()-Dvector3.at(2),2));
                    double deltaR23 = sqrt(pow(d2->eta()-Dvector3.at(1),2)+pow(d2->phi()-Dvector3.at(2),2));
                    double deltaR32 = sqrt(pow(d3->eta()-Dvector2.at(1),2)+pow(d3->phi()-Dvector2.at(2),2));

                    if(!(deltaR22 < deltaR_ && deltaR33 < deltaR_) && !(deltaR23 < deltaR_ && deltaR32 < deltaR_) ) continue;

                    double deltaPt22 = fabs((d2->pt()-Dvector2.at(0))/d2->pt());
                    double deltaPt33 = fabs((d3->pt()-Dvector3.at(0))/d3->pt());
                    double deltaPt23 = fabs((d2->pt()-Dvector3.at(0))/d2->pt());
                    double deltaPt32 = fabs((d3->pt()-Dvector2.at(0))/d3->pt());

                    if( !(deltaPt22 < 0.5 && deltaPt33 < 0.5) && !(deltaPt23 < 0.5 && deltaPt32 < 0.5) ) continue; //check deltaPt matching

                    d2massGEN = Dvector2.at(4);
                    d2mass = d2->mass();
                    d3massGEN = Dvector3.at(4);
                    d3mass = d3->mass();

                    matchGEN[it] = true; //matched gen
                  }

                  if(i%3==2)
                  { 
                    vector<double> Dvector2 = (*pVect)[i-2]; //get GEN daugther vector for track2
                    vector<double> Dvector3 = (*pVect)[i-1]; //get GEN daugther vector for track3

                    if(!(d2->charge()==Dvector2.at(3) && d3->charge()==Dvector3.at(3))
                    && !(d3->charge()==Dvector2.at(3) && d2->charge()==Dvector3.at(3))) continue; //check match charge

                    double deltaR22 = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                    double deltaR33 = sqrt(pow(d3->eta()-Dvector3.at(1),2)+pow(d3->phi()-Dvector3.at(2),2));
                    double deltaR23 = sqrt(pow(d2->eta()-Dvector3.at(1),2)+pow(d2->phi()-Dvector3.at(2),2));
                    double deltaR32 = sqrt(pow(d3->eta()-Dvector2.at(1),2)+pow(d3->phi()-Dvector2.at(2),2));

                    if(!(deltaR22 < deltaR_ && deltaR33 < deltaR_) && !(deltaR23 < deltaR_ && deltaR32 < deltaR_) ) continue;

                    double deltaPt22 = fabs((d2->pt()-Dvector2.at(0))/d2->pt());
                    double deltaPt33 = fabs((d3->pt()-Dvector3.at(0))/d3->pt());
                    double deltaPt23 = fabs((d2->pt()-Dvector3.at(0))/d2->pt());
                    double deltaPt32 = fabs((d3->pt()-Dvector2.at(0))/d3->pt());

                    if( !(deltaPt22 < 0.5 && deltaPt33 < 0.5) && !(deltaPt23 < 0.5 && deltaPt32 < 0.5) ) continue; //check deltaPt matching

                    d2massGEN = Dvector2.at(4);
                    d2mass = d2->mass();
                    d3massGEN = Dvector3.at(4);
                    d3mass = d3->mass();

                    matchGEN[it] = true; //matched gen
                  }

                  //check swap
                  if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01 || abs(d3massGEN - d3mass)>0.01) isSwap[it] = true;
                  
                  //check prompt & record mom id
                  idmom_reco[it] = pVectIDmom->at(i/3);
                }

            }
        }
        
        double pxd1 = d1->px();
        double pyd1 = d1->py();
        double pzd1 = d1->pz();
        double pxd2 = d2->px();
        double pyd2 = d2->py();
        double pzd2 = d2->pz();
        
        TVector3 dauvec1(pxd1,pyd1,pzd1);
        TVector3 dauvec2(pxd2,pyd2,pzd2);
        
        //pt
        pt1[it] = d1->pt();
        pt2[it] = d2->pt();
        
        //momentum
        p1[it] = d1->p();
        p2[it] = d2->p();
        
        //eta
        eta1[it] = d1->eta();
        eta2[it] = d2->eta();
        
        //phi
        phi1[it] = d1->phi();
        phi2[it] = d2->phi();
        
        //charge
        charge1[it] = d1->charge();
        charge2[it] = d2->charge();
        
        double pxd3 = -999.9;
        double pyd3 = -999.9;
        double pzd3 = -999.9;
        if(threeProngDecay_ && d3)
        {
          pxd3 = d3->px();
          pyd3 = d3->py();
          pzd3 = d3->pz();
          pt3[it] = d3->pt();
          p3[it] = d3->p();
          eta3[it] = d3->eta();
          phi3[it] = d3->phi();
          charge3[it] = d3->charge();
        }
        TVector3 dauvec3(pxd3,pyd3,pzd3);

        pid1[it] = -99999;
        pid2[it] = -99999;
        if(doGenMatchingTOF_)
        {
          for(unsigned it=0; it<genpars->size(); ++it){

              const reco::GenParticle & trk = (*genpars)[it];

              if(trk.pt()<0.001) continue;

              int id = trk.pdgId();
              TVector3 trkvect(trk.px(),trk.py(),trk.pz());

              if(fabs(id)!=PID_ && trk.charge())
              {
                // matching daughter 1
                double deltaR = trkvect.DeltaR(dauvec1);
                if(deltaR < deltaR_ && fabs((trk.pt()-pt1[it])/pt1[it]) < 0.5 && trk.charge()==charge1[it] && pid1[it]==-99999)
                {
                  pid1[it] = id;
                } 

                // matching daughter 2
                deltaR = trkvect.DeltaR(dauvec2);
                if(deltaR < deltaR_ && fabs((trk.pt()-pt2[it])/pt2[it]) < 0.5 && trk.charge()==charge2[it] && pid2[it]==-99999)
                {
                  pid2[it] = id;
                }
              }

              if(fabs(id)==PID_ && trk.numberOfDaughters()==2)
              {
                const reco::Candidate * Dd1 = trk.daughter(0);
                const reco::Candidate * Dd2 = trk.daughter(1);
                TVector3 d1vect(Dd1->px(),Dd1->py(),Dd1->pz());
                TVector3 d2vect(Dd2->px(),Dd2->py(),Dd2->pz());
                int id1 = Dd1->pdgId();
                int id2 = Dd2->pdgId();
               
                double deltaR = d1vect.DeltaR(dauvec1);
                if(deltaR < deltaR_ && fabs((Dd1->pt()-pt1[it])/pt1[it]) < 0.5 && Dd1->charge()==charge1[it] && pid1[it]==-99999)
                {
                  pid1[it] = id1;
                }
                deltaR = d2vect.DeltaR(dauvec1);
                if(deltaR < deltaR_ && fabs((Dd2->pt()-pt1[it])/pt1[it]) < 0.5 && Dd2->charge()==charge1[it] && pid1[it]==-99999)
                {
                  pid1[it] = id1;
                }

                deltaR = d1vect.DeltaR(dauvec2);
                if(deltaR < deltaR_ && fabs((Dd1->pt()-pt2[it])/pt2[it]) < 0.5 && Dd1->charge()==charge2[it] && pid2[it]==-99999)
                {
                  pid2[it] = id2;
                }
                deltaR = d2vect.DeltaR(dauvec2);
                if(deltaR < deltaR_ && fabs((Dd2->pt()-pt2[it])/pt2[it]) < 0.5 && Dd2->charge()==charge2[it] && pid2[it]==-99999)
                {
                  pid2[it] = id2;
                }
              }

              if(pid1[it]!=-99999 && pid2[it]!=-99999) break;
          }
        }

        //vtxChi2
        //vtxChi2[it] = trk.vertexChi2();
        //ndf[it] = trk.vertexNdof();
        vtxChi2[it] = trk.userFloat("vertexChi2");
        ndf[it] = trk.userFloat("vertexNdof");
        VtxProb[it] = TMath::Prob(vtxChi2[it],ndf[it]);
        
        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
        TVector3 secvec2D(px,py,0);
        
        agl[it] = cos(secvec.Angle(ptosvec));
        agl_abs[it] = secvec.Angle(ptosvec);
        
        agl2D[it] = cos(secvec2D.Angle(ptosvec2D));
        agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
        
        //Decay length 3D
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        typedef ROOT::Math::SVector<double, 6> SVector6;
        
        //SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SMatrixSym3D totalCov = vtx.covariance() + *trk.userData<reco::Vertex::CovarianceMatrix>("vertexCovariance");
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        dl[it] = ROOT::Math::Mag(distanceVector);
        dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl[it];
        
        dlos[it] = dl[it]/dlerror[it];
        
        //Decay length 2D
        SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);

        const SMatrixSym3D& trkCovMat = *trk.userData<reco::Vertex::CovarianceMatrix>("vertexCovariance");
        SVector6 v2(trkCovMat(0,0), trkCovMat(0,1), trkCovMat(1,1), 0, 0, 0);
        //SVector6 v2(trk.userData<reco::Vertex::CovarianceMatrix>("vertexCovariance")(0,0), 
        //            trk.userData<reco::Vertex::CovarianceMatrix>("vertexCovariance")(0,1),
        //            trk.userData<reco::Vertex::CovarianceMatrix>("vertexCovariance")(1,1),0,0,0);
        
        SMatrixSym3D sv1(v1);
        SMatrixSym3D sv2(v2);
        
        SMatrixSym3D totalCov2D = sv1 + sv2;
        SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
        
        dl2D[it] = ROOT::Math::Mag(distanceVector2D);
        double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D[it];
        
        dlos2D[it] = dl2D[it]/dl2Derror;

        //trk info
        auto dau1 = d1->get<reco::TrackRef>();
        if(!twoLayerDecay_)
        {
            //trk quality
            trkquality1[it] = dau1->quality(reco::TrackBase::highPurity);
            
            //trk dEdx
            H2dedx1[it] = -999.9;
            
            if(dEdxHandle1.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                H2dedx1[it] = dEdxTrack[dau1].dEdx();
            }
            
            T4dedx1[it] = -999.9;
            
            if(dEdxHandle2.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                T4dedx1[it] = dEdxTrack[dau1].dEdx();
            }
            
            //track Chi2
            trkChi1[it] = dau1->normalizedChi2();
            
            //track pT error
            ptErr1[it] = dau1->ptError();
            
            //vertexCovariance 00-xError 11-y 22-z
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            
            //trkNHits
            nhit1[it] = dau1->numberOfValidHits();
            
            //DCA
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzbest1 = dau1->dz(bestvtx);
            double dxybest1 = dau1->dxy(bestvtx);
            double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
            double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
            
            dzos1[it] = dzbest1/dzerror1;
            dxyos1[it] = dxybest1/dxyerror1;
        }
        
        auto dau2 = d2->get<reco::TrackRef>();
        
        //trk quality
        trkquality2[it] = dau2->quality(reco::TrackBase::highPurity);
        
        //trk dEdx
        H2dedx2[it] = -999.9;
        
        if(dEdxHandle1.isValid()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
            H2dedx2[it] = dEdxTrack[dau2].dEdx();
        }
        
        T4dedx2[it] = -999.9;
        
        if(dEdxHandle2.isValid()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
            T4dedx2[it] = dEdxTrack[dau2].dEdx();
        }
        
        //track Chi2
        trkChi2[it] = dau2->normalizedChi2();
        
        //track pT error
        ptErr2[it] = dau2->ptError();
        
        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        
        //trkNHits
        nhit2[it] = dau2->numberOfValidHits();
        
        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        dzos2[it] = dzbest2/dzerror2;
        dxyos2[it] = dxybest2/dxyerror2;
        
        if(doMuon_)
        {
          edm::Handle<reco::MuonCollection> theMuonHandle;
          iEvent.getByToken(tok_muon_, theMuonHandle);
            
          nmatchedch1[it] = -1;
          nmatchedst1[it] = -1;
          matchedenergy1[it] = -1;
          nmatchedch2[it] = -1;
          nmatchedst2[it] = -1;
          matchedenergy2[it] = -1;
            
          double x_exp = -999.;
          double y_exp = -999.;
          double xerr_exp = -999.;
          double yerr_exp = -999.;
          double dxdz_exp = -999.;
          double dydz_exp = -999.;
          double dxdzerr_exp = -999.;
          double dydzerr_exp = -999.;
            
          double x_seg = -999.;
          double y_seg = -999.;
          double xerr_seg = -999.;
          double yerr_seg = -999.;
          double dxdz_seg = -999.;
          double dydz_seg = -999.;
          double dxdzerr_seg = -999.;
          double dydzerr_seg = -999.;
            
          double dx_seg = 999.;
          double dy_seg = 999.;
          double dxerr_seg = 999.;
          double dyerr_seg = 999.;
          double dxSig_seg = 999.;
          double dySig_seg = 999.;
          double ddxdz_seg = 999.;
          double ddydz_seg = 999.;
          double ddxdzerr_seg = 999.;
          double ddydzerr_seg = 999.;
          double ddxdzSig_seg = 999.;
          double ddydzSig_seg = 999.;
            
          onestmuon1[it] = false;
          pfmuon1[it] = false;
          glbmuon1[it] = false;
          trkmuon1[it] = false;
          calomuon1[it] = false; 
          onestmuon2[it] = false;
          pfmuon2[it] = false;
          glbmuon2[it] = false;
          trkmuon2[it] = false;
          calomuon2[it] = false;

          const int muId1 = muAssocToTrack( dau1, theMuonHandle );
          const int muId2 = muAssocToTrack( dau2, theMuonHandle );

          if( muId1 != -1 )
          {
            const reco::Muon& cand = (*theMuonHandle)[muId1];

            onestmuon1[it] = muon::isGoodMuon(cand, muon::selectionTypeFromString("TMOneStationTight"));
            pfmuon1[it] =  cand.isPFMuon();
            glbmuon1[it] =  cand.isGlobalMuon();
            trkmuon1[it] =  cand.isTrackerMuon();
            calomuon1[it] =  cand.isCaloMuon();
          }

          if( muId2 != -1 )
          {
            const reco::Muon& cand = (*theMuonHandle)[muId2];

            onestmuon2[it] = muon::isGoodMuon(cand, muon::selectionTypeFromString("TMOneStationTight"));
            pfmuon2[it] =  cand.isPFMuon();
            glbmuon2[it] =  cand.isGlobalMuon();
            trkmuon2[it] =  cand.isTrackerMuon();
            calomuon2[it] =  cand.isCaloMuon();
          }

          if(doMuonFull_)
          {

          if( muId1 != -1 )
          {
            const reco::Muon& cand = (*theMuonHandle)[muId1];

            nmatchedch1[it] = cand.numberOfMatches();
            nmatchedst1[it] = cand.numberOfMatchedStations();
                   
            reco::MuonEnergy muenergy = cand.calEnergy();
            matchedenergy1[it] = muenergy.hadMax;
                    
            const std::vector<reco::MuonChamberMatch>& muchmatches = cand.matches();
                   
            for(unsigned int ich=0;ich<muchmatches.size();ich++)
            {
              x_exp = muchmatches[ich].x;
              y_exp = muchmatches[ich].y;
              xerr_exp = muchmatches[ich].xErr;
              yerr_exp = muchmatches[ich].yErr;
              dxdz_exp = muchmatches[ich].dXdZ;
              dydz_exp = muchmatches[ich].dYdZ;
              dxdzerr_exp = muchmatches[ich].dXdZErr;
              dydzerr_exp = muchmatches[ich].dYdZErr;
                        
              std::vector<reco::MuonSegmentMatch> musegmatches = muchmatches[ich].segmentMatches;
                        
              if(!musegmatches.size()) continue;
              for(unsigned int jseg=0;jseg<musegmatches.size();jseg++)
              {
                x_seg = musegmatches[jseg].x;
                y_seg = musegmatches[jseg].y;
                xerr_seg = musegmatches[jseg].xErr;
                yerr_seg = musegmatches[jseg].yErr;
                dxdz_seg = musegmatches[jseg].dXdZ;
                dydz_seg = musegmatches[jseg].dYdZ;
                dxdzerr_seg = musegmatches[jseg].dXdZErr;
                dydzerr_seg = musegmatches[jseg].dYdZErr;
                            
                if(sqrt((x_seg-x_exp)*(x_seg-x_exp)+(y_seg-y_exp)*(y_seg-y_exp))<sqrt(dx_seg*dx_seg+dy_seg*dy_seg))
                {
                  dx_seg = x_seg - x_exp;
                  dy_seg = y_seg - y_exp;
                  dxerr_seg = sqrt(xerr_seg*xerr_seg+xerr_exp*xerr_exp);
                  dyerr_seg = sqrt(yerr_seg*yerr_seg+yerr_exp*yerr_exp);
                  dxSig_seg = dx_seg / dxerr_seg;
                  dySig_seg = dy_seg / dyerr_seg;
                  ddxdz_seg = dxdz_seg - dxdz_exp;
                  ddydz_seg = dydz_seg - dydz_exp;
                  ddxdzerr_seg = sqrt(dxdzerr_seg*dxdzerr_seg+dxdzerr_exp*dxdzerr_exp);
                  ddydzerr_seg = sqrt(dydzerr_seg*dydzerr_seg+dydzerr_exp*dydzerr_exp);
                  ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                  ddydzSig_seg = ddydz_seg / ddydzerr_seg;
                }
              }
                      
              dx1_seg_[it]=dx_seg;
              dy1_seg_[it]=dy_seg;
              dxSig1_seg_[it]=dxSig_seg;
              dySig1_seg_[it]=dySig_seg;
              ddxdz1_seg_[it]=ddxdz_seg;
              ddydz1_seg_[it]=ddydz_seg;
              ddxdzSig1_seg_[it]=ddxdzSig_seg;
              ddydzSig1_seg_[it]=ddydzSig_seg;
            }
          } 

          if( muId2 != -1 )
          {
            const reco::Muon& cand = (*theMuonHandle)[muId2];

            nmatchedch2[it] = cand.numberOfMatches();
            nmatchedst2[it] = cand.numberOfMatchedStations();
                    
            reco::MuonEnergy muenergy = cand.calEnergy();
            matchedenergy2[it] = muenergy.hadMax;
                    
            const std::vector<reco::MuonChamberMatch>& muchmatches = cand.matches();
            for(unsigned int ich=0;ich<muchmatches.size();ich++)
                        //                        for(unsigned int ich=0;ich<1;ich++)
            {
              x_exp = muchmatches[ich].x;
              y_exp = muchmatches[ich].y;
              xerr_exp = muchmatches[ich].xErr;
              yerr_exp = muchmatches[ich].yErr;
              dxdz_exp = muchmatches[ich].dXdZ;
              dydz_exp = muchmatches[ich].dYdZ;
              dxdzerr_exp = muchmatches[ich].dXdZErr;
              dydzerr_exp = muchmatches[ich].dYdZErr;
                        
              std::vector<reco::MuonSegmentMatch> musegmatches = muchmatches[ich].segmentMatches;
                        
              if(!musegmatches.size()) continue;
              for(unsigned int jseg=0;jseg<musegmatches.size();jseg++)
              {
                x_seg = musegmatches[jseg].x;
                y_seg = musegmatches[jseg].y;
                xerr_seg = musegmatches[jseg].xErr;
                yerr_seg = musegmatches[jseg].yErr;
                dxdz_seg = musegmatches[jseg].dXdZ;
                dydz_seg = musegmatches[jseg].dYdZ;
                dxdzerr_seg = musegmatches[jseg].dXdZErr;
                dydzerr_seg = musegmatches[jseg].dYdZErr;
                            
                if(sqrt((x_seg-x_exp)*(x_seg-x_exp)+(y_seg-y_exp)*(y_seg-y_exp))<sqrt(dx_seg*dx_seg+dy_seg*dy_seg))
                {
                  dx_seg = x_seg - x_exp;
                  dy_seg = y_seg - y_exp;
                  dxerr_seg = sqrt(xerr_seg*xerr_seg+xerr_exp*xerr_exp);
                  dyerr_seg = sqrt(yerr_seg*yerr_seg+yerr_exp*yerr_exp);
                  dxSig_seg = dx_seg / dxerr_seg;
                  dySig_seg = dy_seg / dyerr_seg;
                  ddxdz_seg = dxdz_seg - dxdz_exp;
                  ddydz_seg = dydz_seg - dydz_exp;
                  ddxdzerr_seg = sqrt(dxdzerr_seg*dxdzerr_seg+dxdzerr_exp*dxdzerr_exp);
                  ddydzerr_seg = sqrt(dydzerr_seg*dydzerr_seg+dydzerr_exp*dydzerr_exp);
                  ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                  ddydzSig_seg = ddydz_seg / ddydzerr_seg;
                }
              }
                        
              dx2_seg_[it]=dx_seg;
              dy2_seg_[it]=dy_seg;
              dxSig2_seg_[it]=dxSig_seg;
              dySig2_seg_[it]=dySig_seg;
              ddxdz2_seg_[it]=ddxdz_seg;
              ddydz2_seg_[it]=ddydz_seg;
              ddxdzSig2_seg_[it]=ddxdzSig_seg;
              ddydzSig2_seg_[it]=ddydzSig_seg;
            }
          }
          } // doMuonFull
        }
        
        if(twoLayerDecay_)
        {
            grand_mass[it] = d1->mass();
            
            const reco::Candidate * gd1 = d1->daughter(0);
            const reco::Candidate * gd2 = d1->daughter(1);
            
            double gpxd1 = gd1->px();
            double gpyd1 = gd1->py();
            double gpzd1 = gd1->pz();
            double gpxd2 = gd2->px();
            double gpyd2 = gd2->py();
            double gpzd2 = gd2->pz();
            
            TVector3 gdauvec1(gpxd1,gpyd1,gpzd1);
            TVector3 gdauvec2(gpxd2,gpyd2,gpzd2);
            
            auto gdau1 = gd1->get<reco::TrackRef>();
            auto gdau2 = gd2->get<reco::TrackRef>();
            
            //trk quality
            
            grand_trkquality1[it] = gdau1->quality(reco::TrackBase::highPurity);
            grand_trkquality2[it] = gdau2->quality(reco::TrackBase::highPurity);
            
            //trk dEdx
            grand_H2dedx1[it] = -999.9;
            grand_H2dedx2[it] = -999.9;
            
            if(dEdxHandle1.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                grand_H2dedx1[it] = dEdxTrack[gdau1].dEdx();
                grand_H2dedx2[it] = dEdxTrack[gdau2].dEdx();
            }
            
            grand_T4dedx1[it] = -999.9;
            grand_T4dedx2[it] = -999.9;
            
            if(dEdxHandle2.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                grand_T4dedx1[it] = dEdxTrack[gdau1].dEdx();
                grand_T4dedx2[it] = dEdxTrack[gdau2].dEdx();
            }
            
            //track pt
            grand_pt1[it] = gd1->pt();
            grand_pt2[it] = gd2->pt();
            
            //track momentum
            grand_p1[it] = gd1->p();
            grand_p2[it] = gd2->p();
            
            //track eta
            grand_eta1[it] = gd1->eta();
            grand_eta2[it] = gd2->eta();
            
            //track charge
            grand_charge1[it] = gd1->charge();
            grand_charge2[it] = gd2->charge();
            
            //track Chi2
            grand_trkChi1[it] = gdau1->normalizedChi2();
            grand_trkChi2[it] = gdau2->normalizedChi2();
            
            //track pT error
            grand_ptErr1[it] = gdau1->ptError();
            grand_ptErr2[it] = gdau2->ptError();
            
            //vertexCovariance 00-xError 11-y 22-z
            secvz = d1->vz(); secvx = d1->vx(); secvy = d1->vy();
            
            //trkNHits
            grand_nhit1[it] = gdau1->numberOfValidHits();
            grand_nhit2[it] = gdau2->numberOfValidHits();
            
            //DCA
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double gdzbest1 = gdau1->dz(bestvtx);
            double gdxybest1 = gdau1->dxy(bestvtx);
            double gdzerror1 = sqrt(gdau1->dzError()*gdau1->dzError()+bestvzError*bestvzError);
            double gdxyerror1 = sqrt(gdau1->d0Error()*gdau1->d0Error()+bestvxError*bestvyError);
            
            grand_dzos1[it] = gdzbest1/gdzerror1;
            grand_dxyos1[it] = gdxybest1/gdxyerror1;
            
            double gdzbest2 = gdau2->dz(bestvtx);
            double gdxybest2 = gdau2->dxy(bestvtx);
            double gdzerror2 = sqrt(gdau2->dzError()*gdau2->dzError()+bestvzError*bestvzError);
            double gdxyerror2 = sqrt(gdau2->d0Error()*gdau2->d0Error()+bestvxError*bestvyError);
            
            grand_dzos2[it] = gdzbest2/gdzerror2;
            grand_dxyos2[it] = gdxybest2/gdxyerror2;
            
            //vtxChi2
            grand_vtxChi2[it] = d1->vertexChi2();
            grand_ndf[it] = d1->vertexNdof();
            grand_VtxProb[it] = TMath::Prob(grand_vtxChi2[it],grand_ndf[it]);
            
            //PAngle
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(d1->px(),d1->py(),d1->pz());
            
            TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
            TVector3 secvec2D(d1->px(),d1->py(),0);
            
            grand_agl[it] = cos(secvec.Angle(ptosvec));
            grand_agl_abs[it] = secvec.Angle(ptosvec);
            
            grand_agl2D[it] = cos(secvec2D.Angle(ptosvec2D));
            grand_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
            
            //Decay length 3D
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            typedef ROOT::Math::SVector<double, 6> SVector6;
            
            SMatrixSym3D totalCov = vtx.covariance() + d1->vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            grand_dl[it] = ROOT::Math::Mag(distanceVector);
            grand_dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/grand_dl[it];
            
            grand_dlos[it] = grand_dl[it]/grand_dlerror[it];
            
            //Decay length 2D
            SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
            SVector6 v2(d1->vertexCovariance(0,0), d1->vertexCovariance(0,1),d1->vertexCovariance(1,1),0,0,0);
            
            SMatrixSym3D sv1(v1);
            SMatrixSym3D sv2(v2);
            
            SMatrixSym3D totalCov2D = sv1 + sv2;
            SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
            
            double gdl2D = ROOT::Math::Mag(distanceVector2D);
            double gdl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/gdl2D;
            
            grand_dlos2D[it] = gdl2D/gdl2Derror;
        }

        if(saveHistogram_)
        {
          for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
            for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
            {
              if(pt[it]<pTBins_[ipt+1] && pt[it]>pTBins_[ipt] && y[it]<yBins_[iy+1] && y[it]>yBins_[iy])
              {
                hMassVsMVA[iy][ipt]->Fill(mva[it],mass[it]);
//                h3DDCAVsMVA[iy][ipt]->Fill(mva[it],dl[it]*sin(agl_abs[it]));
//                h2DDCAVsMVA[iy][ipt]->Fill(mva[it],dl2D[it]*sin(agl2D_abs[it]));

                if(saveAllHistogram_)
                {
                hpTVsMVA[iy][ipt]->Fill(mva[it],pt[it]);
                hetaVsMVA[iy][ipt]->Fill(mva[it],eta[it]);
                hyVsMVA[iy][ipt]->Fill(mva[it],y[it]);
                hVtxProbVsMVA[iy][ipt]->Fill(mva[it],VtxProb[it]);
                h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it],agl[it]);
                h3DPointingAngleVsMVA[iy][ipt]->Fill(mva[it],agl_abs[it]);
                h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it],agl2D[it]);
                h2DPointingAngleVsMVA[iy][ipt]->Fill(mva[it],agl2D_abs[it]);
                h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it],dlos[it]);
                h3DDecayLengthVsMVA[iy][ipt]->Fill(mva[it],dl[it]);
                h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it],dlos2D[it]);
                h2DDecayLengthVsMVA[iy][ipt]->Fill(mva[it],dl2D[it]);
                hzDCASignificanceDaugther1VsMVA[iy][ipt]->Fill(mva[it],dzos1[it]);
                hxyDCASignificanceDaugther1VsMVA[iy][ipt]->Fill(mva[it],dxyos1[it]);
                hNHitD1VsMVA[iy][ipt]->Fill(mva[it],nhit1[it]);
                hpTD1VsMVA[iy][ipt]->Fill(mva[it],pt1[it]);
                hpTerrD1VsMVA[iy][ipt]->Fill(mva[it],ptErr1[it]/pt1[it]);
                hEtaD1VsMVA[iy][ipt]->Fill(mva[it],eta1[it]);
                hdedxHarmonic2D1VsMVA[iy][ipt]->Fill(mva[it],H2dedx1[it]);
                hdedxHarmonic2D1VsP[iy][ipt]->Fill(p1[it],H2dedx1[it]);
                hzDCASignificanceDaugther2VsMVA[iy][ipt]->Fill(mva[it],dzos2[it]);
                hxyDCASignificanceDaugther2VsMVA[iy][ipt]->Fill(mva[it],dxyos2[it]);
                hNHitD2VsMVA[iy][ipt]->Fill(mva[it],nhit2[it]);
                hpTD2VsMVA[iy][ipt]->Fill(mva[it],pt2[it]);
                hpTerrD2VsMVA[iy][ipt]->Fill(mva[it],ptErr2[it]/pt2[it]);
                hEtaD2VsMVA[iy][ipt]->Fill(mva[it],eta2[it]);
                hdedxHarmonic2D2VsMVA[iy][ipt]->Fill(mva[it],H2dedx2[it]);
                hdedxHarmonic2D2VsP[iy][ipt]->Fill(p2[it],H2dedx2[it]);
                if(threeProngDecay_)
                {
                  hzDCASignificanceDaugther3VsMVA[iy][ipt]->Fill(mva[it],dzos3[it]);
                  hxyDCASignificanceDaugther3VsMVA[iy][ipt]->Fill(mva[it],dxyos3[it]);
                  hNHitD3VsMVA[iy][ipt]->Fill(mva[it],nhit3[it]);
                  hpTD3VsMVA[iy][ipt]->Fill(mva[it],pt3[it]);
                  hpTerrD3VsMVA[iy][ipt]->Fill(mva[it],ptErr3[it]/pt3[it]);
                  hEtaD3VsMVA[iy][ipt]->Fill(mva[it],eta3[it]);
                  hdedxHarmonic2D3VsMVA[iy][ipt]->Fill(mva[it],H2dedx3[it]);
                  hdedxHarmonic2D3VsP[iy][ipt]->Fill(p1[it],H2dedx3[it]);
                }

                }
              }
            }
        }

    }
}

void
VertexCompositeTreeProducer::fillGEN(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(tok_genParticle_,genpars);

    candSize_gen = 0;
    for(unsigned it=0; it<genpars->size(); ++it){

        const reco::GenParticle & trk = (*genpars)[it];

        int id = trk.pdgId();
        if(fabs(id)!=PID_) continue; //check is target
        if(decayInGen_ && (trk.numberOfDaughters()!=2 || trk.numberOfDaughters()!=3)) continue; //check 2-pron decay if target decays in Gen

        candSize_gen+=1;

        pt_gen[candSize_gen-1] = trk.pt();
        eta_gen[candSize_gen-1] = trk.eta();
        status_gen[candSize_gen-1] = trk.status();
        idmom[candSize_gen-1] = -77;
        y_gen[candSize_gen-1] = trk.rapidity();

        if(trk.numberOfMothers()!=0)
        {
            const reco::Candidate * mom = trk.mother();
            idmom[candSize_gen-1] = mom->pdgId();
        }

        if(!decayInGen_) continue;

        const reco::Candidate * Dd1 = trk.daughter(0);
        const reco::Candidate * Dd2 = trk.daughter(1);
        const reco::Candidate * Dd3 = trk.daughter(2);

        iddau1[candSize_gen-1] = fabs(Dd1->pdgId());
        iddau2[candSize_gen-1] = fabs(Dd2->pdgId());
        if(Dd3) iddau3[candSize_gen-1] = fabs(Dd3->pdgId());
    }
}

// ------------ method called once each job just before starting event
//loop  ------------
void
VertexCompositeTreeProducer::beginJob()
{
    TH1D::SetDefaultSumw2();
    
    if(!doRecoNtuple_ && !doGenNtuple_)
    {
        cout<<"No output for either RECO or GEN!! Fix config!!"<<endl; return;
    }

    if(twoLayerDecay_ && doMuon_)
    {
        cout<<"Muons cannot be coming from two layer decay!! Fix config!!"<<endl; return;
    }
    
    if(saveHistogram_) initHistogram();
    if(saveTree_) initTree();
}

void
VertexCompositeTreeProducer::initHistogram()
{
  for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
  {
    for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
  {
   hMassVsMVA[iy][ipt] = fs->make<TH2F>(Form("hMassVsMVA_y%d_pt%d",iy,ipt),";mva;mass(GeV)",100,-1.,1.,massHistBins_,massHistPeak_-massHistWidth_,massHistPeak_+massHistWidth_);
//   h3DDCAVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;3D DCA;",100,-1.,1.,1000,0,10);
//   h2DDCAVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;2D DCA;",100,-1.,1.,1000,0,10);

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
   hzDCASignificanceDaugther1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificanceDaugther1VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificanceDaugther1;",100,-1.,1.,100,-10,10);
   hxyDCASignificanceDaugther1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificanceDaugther1VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificanceDaugther1;",100,-1.,1.,100,-10,10);
   hNHitD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hNHitD1VsMVA_y%d_pt%d",iy,ipt),";mva;NHitD1;",100,-1.,1.,100,0,100);
   hpTD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTD1VsMVA_y%d_pt%d",iy,ipt),";mva;pTD1;",100,-1.,1.,100,0,10);
   hpTerrD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTerrD1VsMVA_y%d_pt%d",iy,ipt),";mva;pTerrD1;",100,-1.,1.,50,0,0.5);
   hEtaD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hEtaD1VsMVA_y%d_pt%d",iy,ipt),";mva;EtaD1;",100,-1.,1.,40,-4,4);
   hdedxHarmonic2D1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D1VsMVA_y%d_pt%d",iy,ipt),";mva;dedxHarmonic2D1;",100,-1.,1.,100,0,10);
   hdedxHarmonic2D1VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D1VsP_y%d_pt%d",iy,ipt),";p (GeV);dedxHarmonic2D1",100,0,10,100,0,10);
   hzDCASignificanceDaugther2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificanceDaugther2VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificanceDaugther2;",100,-1.,1.,100,-10,10);
   hxyDCASignificanceDaugther2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificanceDaugther2VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificanceDaugther2;",100,-1.,1.,100,-10,10);
   hNHitD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hNHitD2VsMVA_y%d_pt%d",iy,ipt),";mva;NHitD2;",100,-1.,1.,100,0,100);
   hpTD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTD2VsMVA_y%d_pt%d",iy,ipt),";mva;pTD2;",100,-1.,1.,100,0,10);
   hpTerrD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTerrD2VsMVA_y%d_pt%d",iy,ipt),";mva;pTerrD2;",100,-1.,1.,50,0,0.5);
   hEtaD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hEtaD2VsMVA_y%d_pt%d",iy,ipt),";mva;EtaD2;",100,-1.,1.,40,-4,4);
   hdedxHarmonic2D2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D2VsMVA_y%d_pt%d",iy,ipt),";mva;dedxHarmonic2D2;",100,-1.,1.,100,0,10);
   hdedxHarmonic2D2VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D2VsP_y%d_pt%d",iy,ipt),";p (GeV);dedxHarmonic2D2",100,0,10,100,0,10);

   if(threeProngDecay_)
   {
     hzDCASignificanceDaugther3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificanceDaugther3VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificanceDaugther3;",100,-1.,1.,100,-10,10);
     hxyDCASignificanceDaugther3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificanceDaugther3VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificanceDaugther3;",100,-1.,1.,100,-10,10);
     hNHitD3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hNHitD3VsMVA_y%d_pt%d",iy,ipt),";mva;NHitD3;",100,-1.,1.,100,0,100);
     hpTD3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTD3VsMVA_y%d_pt%d",iy,ipt),";mva;pTD3;",100,-1.,1.,100,0,10);
     hpTerrD3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTerrD3VsMVA_y%d_pt%d",iy,ipt),";mva;pTerrD3;",100,-1.,1.,50,0,0.5);
     hEtaD3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hEtaD3VsMVA_y%d_pt%d",iy,ipt),";mva;EtaD3;",100,-1.,1.,40,-4,4);
     hdedxHarmonic2D3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D3VsMVA_y%d_pt%d",iy,ipt),";mva;dedxHarmonic2D3;",100,-1.,1.,100,0,10);
     hdedxHarmonic2D3VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D3VsP_y%d_pt%d",iy,ipt),";p (GeV);dedxHarmonic2D3",100,0,10,100,0,10);
   }

   }
  }
 }
}

void 
VertexCompositeTreeProducer::initTree()
{ 
    VertexCompositeNtuple = fs->make< TTree>("VertexCompositeNtuple","VertexCompositeNtuple");
    
    if(doRecoNtuple_) 
    { 
  
    // Event info
    VertexCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
    VertexCompositeNtuple->Branch("Npixel",&Npixel,"Npixel/I");
    VertexCompositeNtuple->Branch("HFsumET",&HFsumET,"HFsumET/F");
    VertexCompositeNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
    VertexCompositeNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
    VertexCompositeNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");
    VertexCompositeNtuple->Branch("candSize",&candSize,"candSize/I");
    if(isCentrality_) VertexCompositeNtuple->Branch("centrality",&centrality,"centrality/I");

    // particle info
    VertexCompositeNtuple->Branch("pT",&pt,"pT[candSize]/F");
    VertexCompositeNtuple->Branch("y",&y,"y[candSize]/F");
    VertexCompositeNtuple->Branch("mass",&mass,"mass[candSize]/F");
    if(useAnyMVA_) VertexCompositeNtuple->Branch("mva",&mva,"mva[candSize]/F");

    if(!isSkimMVA_)  
    {
        //Composite candidate info RECO
        VertexCompositeNtuple->Branch("flavor",&flavor,"flavor[candSize]/F");
        VertexCompositeNtuple->Branch("eta",&eta,"eta[candSize]/F");
        VertexCompositeNtuple->Branch("VtxProb",&VtxProb,"VtxProb[candSize]/F");
//        VertexCompositeNtuple->Branch("VtxChi2",&vtxChi2,"VtxChi2[candSize]/F");
//        VertexCompositeNtuple->Branch("VtxNDF",&ndf,"VtxNDF[candSize]/F");
        VertexCompositeNtuple->Branch("3DCosPointingAngle",&agl,"3DCosPointingAngle[candSize]/F");
        VertexCompositeNtuple->Branch("3DPointingAngle",&agl_abs,"3DPointingAngle[candSize]/F");
        VertexCompositeNtuple->Branch("2DCosPointingAngle",&agl2D,"2DCosPointingAngle[candSize]/F");
        VertexCompositeNtuple->Branch("2DPointingAngle",&agl2D_abs,"2DPointingAngle[candSize]/F");
        VertexCompositeNtuple->Branch("3DDecayLengthSignificance",&dlos,"3DDecayLengthSignificance[candSize]/F");
        VertexCompositeNtuple->Branch("3DDecayLength",&dl,"3DDecayLength[candSize]/F");
//        VertexCompositeNtuple->Branch("3DDecayLengthError",&dlerror,"3DDecayLengthError[candSize]/F");
        VertexCompositeNtuple->Branch("2DDecayLengthSignificance",&dlos2D,"2DDecayLengthSignificance[candSize]/F");
        VertexCompositeNtuple->Branch("2DDecayLength",&dl2D,"2DDecayLength[candSize]/F");
    
        if(doGenMatching_)
        {
            VertexCompositeNtuple->Branch("isSwap",&isSwap,"isSwap[candSize]/O");
            VertexCompositeNtuple->Branch("idmom_reco",&idmom_reco,"idmom_reco[candSize]/I");
            VertexCompositeNtuple->Branch("matchGEN",&matchGEN,"matchGEN[candSize]/O");
        }
        
        if(doGenMatchingTOF_)
        {
          VertexCompositeNtuple->Branch("PIDD1",&pid1,"PIDD1[candSize]/I");
          VertexCompositeNtuple->Branch("PIDD2",&pid1,"PIDD2[candSize]/I");
          VertexCompositeNtuple->Branch("TOFD1",&tof1,"TOFD1[candSize]/F");
          VertexCompositeNtuple->Branch("TOFD2",&tof1,"TOFD2[candSize]/F");
        }

        //daugther & grand daugther info
        if(twoLayerDecay_)
        {
            VertexCompositeNtuple->Branch("massDaugther1",&grand_mass,"massDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("pTD1",&pt1,"pTD1[candSize]/F");
            VertexCompositeNtuple->Branch("EtaD1",&eta1,"EtaD1[candSize]/F");
//            VertexCompositeNtuple->Branch("PhiD1",&phi1,"PhiD1[candSize]/F");
            VertexCompositeNtuple->Branch("VtxProbDaugther1",&grand_VtxProb,"VtxProbDaugther1[candSize]/F");
//            VertexCompositeNtuple->Branch("VtxChi2Daugther1",&grand_vtxChi2,"VtxChi2Daugther1[candSize]/F");
//            VertexCompositeNtuple->Branch("VtxNDFDaugther1",&grand_ndf,"VtxNDFDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("3DCosPointingAngleDaugther1",&grand_agl,"3DCosPointingAngleDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("3DPointingAngleDaugther1",&grand_agl_abs,"3DPointingAngleDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("2DCosPointingAngleDaugther1",&grand_agl2D,"2DCosPointingAngleDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("2DPointingAngleDaugther1",&grand_agl2D_abs,"2DPointingAngleDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthSignificanceDaugther1",&grand_dlos,"3DDecayLengthSignificanceDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthDaugther1",&grand_dl,"3DDecayLengthDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthErrorDaugther1",&grand_dlerror,"3DDecayLengthErrorDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("2DDecayLengthSignificanceDaugther1",&grand_dlos2D,"2DDecayLengthSignificanceDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("zDCASignificanceDaugther2",&dzos2,"zDCASignificanceDaugther2[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2",&dxyos2,"xyDCASignificanceDaugther2[candSize]/F");
            VertexCompositeNtuple->Branch("NHitD2",&nhit2,"NHitD2[candSize]/F");
            VertexCompositeNtuple->Branch("HighPurityDaugther2",&trkquality2,"HighPurityDaugther2[candSize]/O");
            VertexCompositeNtuple->Branch("pTD2",&pt2,"pTD2[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrD2",&ptErr2,"pTerrD2[candSize]/F");
//            VertexCompositeNtuple->Branch("pD2",&p2,"pD2[candSize]/F");
            VertexCompositeNtuple->Branch("EtaD2",&eta2,"EtaD2[candSize]/F");
//            VertexCompositeNtuple->Branch("PhiD2",&phi2,"PhiD2[candSize]/F");
//            VertexCompositeNtuple->Branch("chargeD2",&charge2,"chargeD2[candSize]/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2D2",&H2dedx2,"dedxHarmonic2D2[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40Daugther2",&T4dedx2,"dedxTruncated40Daugther2[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2Daugther2",&trkChi2,"normalizedChi2Daugther2[candSize]/F");
            VertexCompositeNtuple->Branch("zDCASignificanceGrandDaugther1",&grand_dzos1,"zDCASignificanceGrandDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("zDCASignificanceGrandDaugther2",&grand_dzos2,"zDCASignificanceGrandDaugther2[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceGrandDaugther1",&grand_dxyos1,"xyDCASignificanceGrandDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceGrandDaugther2",&grand_dxyos2,"xyDCASignificanceGrandDaugther2[candSize]/F");
            VertexCompositeNtuple->Branch("NHitGrandD1",&grand_nhit1,"NHitGrandD1[candSize]/F");
            VertexCompositeNtuple->Branch("NHitGrandD2",&grand_nhit2,"NHitGrandD2[candSize]/F");
            VertexCompositeNtuple->Branch("HighPurityGrandDaugther1",&grand_trkquality1,"HighPurityGrandDaugther1[candSize]/O");
            VertexCompositeNtuple->Branch("HighPurityGrandDaugther2",&grand_trkquality2,"HighPurityGrandDaugther2[candSize]/O");
            VertexCompositeNtuple->Branch("pTGrandD1",&grand_pt1,"pTGrandD1[candSize]/F");
            VertexCompositeNtuple->Branch("pTGrandD2",&grand_pt2,"pTGrandD2[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrGrandD1",&grand_ptErr1,"pTerrGrandD1[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrGrandD2",&grand_ptErr2,"pTerrGrandD2[candSize]/F");
//            VertexCompositeNtuple->Branch("pGrandD1",&grand_p1,"pGrandD1[candSize]/F");
//            VertexCompositeNtuple->Branch("pGrandD2",&grand_p2,"pGrandD2[candSize]/F");
            VertexCompositeNtuple->Branch("EtaGrandD1",&grand_eta1,"EtaGrandD1[candSize]/F");
            VertexCompositeNtuple->Branch("EtaGrandD2",&grand_eta2,"EtaGrandD2[candSize]/F");
//            VertexCompositeNtuple->Branch("chargeGrandD1",&grand_charge1,"chargeGrandD1[candSize]/I");
//            VertexCompositeNtuple->Branch("chargeGrandD2",&grand_charge2,"chargeGrandD2[candSize]/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2GrandD1",&grand_H2dedx1,"dedxHarmonic2GrandD1[candSize]/F");
            VertexCompositeNtuple->Branch("dedxHarmonic2GrandD2",&grand_H2dedx2,"dedxHarmonic2GrandD2[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40GrandDaugther1",&grand_T4dedx1,"dedxTruncated40GrandDaugther1[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40GrandDaugther2",&grand_T4dedx2,"dedxTruncated40GrandDaugther2[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2GrandDaugther1",&grand_trkChi1,"normalizedChi2GrandDaugther1[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2GrandDaugther2",&grand_trkChi2,"normalizedChi2GrandDaugther2[candSize]/F");
        }
        else
        {
            VertexCompositeNtuple->Branch("zDCASignificanceDaugther1",&dzos1,"zDCASignificanceDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceDaugther1",&dxyos1,"xyDCASignificanceDaugther1[candSize]/F");
            VertexCompositeNtuple->Branch("NHitD1",&nhit1,"NHitD1[candSize]/F");
            VertexCompositeNtuple->Branch("HighPurityDaugther1",&trkquality1,"HighPurityDaugther1[candSize]/O");
            VertexCompositeNtuple->Branch("pTD1",&pt1,"pTD1[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrD1",&ptErr1,"pTerrD1[candSize]/F");
//            VertexCompositeNtuple->Branch("pD1",&p1,"pD1[candSize]/F");
            VertexCompositeNtuple->Branch("EtaD1",&eta1,"EtaD1[candSize]/F");
//            VertexCompositeNtuple->Branch("PhiD1",&eta1,"PhiD1[candSize]/F");
//            VertexCompositeNtuple->Branch("chargeD1",&charge1,"chargeD1[candSize]/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2D1",&H2dedx1,"dedxHarmonic2D1[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40Daugther1",&T4dedx1,"dedxTruncated40Daugther1[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2Daugther1",&trkChi1,"normalizedChi2Daugther1[candSize]/F");
            VertexCompositeNtuple->Branch("zDCASignificanceDaugther2",&dzos2,"zDCASignificanceDaugther2[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2",&dxyos2,"xyDCASignificanceDaugther2[candSize]/F");
            VertexCompositeNtuple->Branch("NHitD2",&nhit2,"NHitD2[candSize]/F");
            VertexCompositeNtuple->Branch("HighPurityDaugther2",&trkquality2,"HighPurityDaugther2[candSize]/O");
            VertexCompositeNtuple->Branch("pTD2",&pt2,"pTD2[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrD2",&ptErr2,"pTerrD2[candSize]/F");
//            VertexCompositeNtuple->Branch("pD2",&p2,"pD2[candSize]/F");
            VertexCompositeNtuple->Branch("EtaD2",&eta2,"EtaD2[candSize]/F");
//            VertexCompositeNtuple->Branch("PhiD2",&eta2,"PhiD2[candSize]/F");
//            VertexCompositeNtuple->Branch("chargeD2",&charge2,"chargeD2[candSize]/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2D2",&H2dedx2,"dedxHarmonic2D2[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40Daugther2",&T4dedx2,"dedxTruncated40Daugther2[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2Daugther2",&trkChi2,"normalizedChi2Daugther2[candSize]/F");
            if(threeProngDecay_)
            {
              VertexCompositeNtuple->Branch("zDCASignificanceDaugther3",&dzos3,"zDCASignificanceDaugther3[candSize]/F");
              VertexCompositeNtuple->Branch("xyDCASignificanceDaugther3",&dxyos3,"xyDCASignificanceDaugther3[candSize]/F");
              VertexCompositeNtuple->Branch("NHitD3",&nhit3,"NHitD3[candSize]/F");
              VertexCompositeNtuple->Branch("HighPurityDaugther3",&trkquality3,"HighPurityDaugther3[candSize]/O");
              VertexCompositeNtuple->Branch("pTD3",&pt1,"pTD3[candSize]/F");
              VertexCompositeNtuple->Branch("pTerrD3",&ptErr3,"pTerrD3[candSize]/F");
              VertexCompositeNtuple->Branch("EtaD3",&eta1,"EtaD3[candSize]/F");
              VertexCompositeNtuple->Branch("dedxHarmonic2D3",&H2dedx1,"dedxHarmonic2D3[candSize]/F");
            }
        }
        
        if(doMuon_)
        {
            VertexCompositeNtuple->Branch("OneStMuon1",&onestmuon1,"OneStMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("OneStMuon2",&onestmuon2,"OneStMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("PFMuon1",&pfmuon1,"PFMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("PFMuon2",&pfmuon2,"PFMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("GlbMuon1",&glbmuon1,"GlbMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("GlbMuon2",&glbmuon2,"GlbMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("trkMuon1",&trkmuon1,"trkMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("trkMuon2",&trkmuon2,"trkMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("caloMuon1",&calomuon1,"caloMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("caloMuon2",&calomuon2,"caloMuon2[candSize]/O");
            if(doMuonFull_)
            {
              VertexCompositeNtuple->Branch("nMatchedChamberD1",&nmatchedch1,"nMatchedChamberD1[candSize]/F");
              VertexCompositeNtuple->Branch("nMatchedStationD1",&nmatchedst1,"nMatchedStationD1[candSize]/F");
              VertexCompositeNtuple->Branch("EnergyDepositionD1",&matchedenergy1,"EnergyDepositionD1[candSize]/F");
              VertexCompositeNtuple->Branch("nMatchedChamberD2",&nmatchedch2,"nMatchedChamberD2[candSize]/F");
              VertexCompositeNtuple->Branch("nMatchedStationD2",&nmatchedst2,"nMatchedStationD2[candSize]/F");
              VertexCompositeNtuple->Branch("EnergyDepositionD2",&matchedenergy2,"EnergyDepositionD2[candSize]/F");
              VertexCompositeNtuple->Branch("dx1_seg",        &dx1_seg_, "dx1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dy1_seg",        &dy1_seg_, "dy1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dxSig1_seg",     &dxSig1_seg_, "dxSig1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dySig1_seg",     &dySig1_seg_, "dySig1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddxdz1_seg",     &ddxdz1_seg_, "ddxdz1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddydz1_seg",     &ddydz1_seg_, "ddydz1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddxdzSig1_seg",  &ddxdzSig1_seg_, "ddxdzSig1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddydzSig1_seg",  &ddydzSig1_seg_, "ddydzSig1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dx2_seg",        &dx2_seg_, "dx2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dy2_seg",        &dy2_seg_, "dy2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dxSig2_seg",     &dxSig2_seg_, "dxSig2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dySig2_seg",     &dySig2_seg_, "dySig2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddxdz2_seg",     &ddxdz2_seg_, "ddxdz2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddydz2_seg",     &ddydz2_seg_, "ddydz2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddxdzSig2_seg",  &ddxdzSig2_seg_, "ddxdzSig2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddydzSig2_seg",  &ddydzSig2_seg_, "ddydzSig2_seg[candSize]/F");
           }
        }
    }

    } // doRecoNtuple_

    if(doGenNtuple_)
    {
        VertexCompositeNtuple->Branch("candSize_gen",&candSize_gen,"candSize_gen/I");
        VertexCompositeNtuple->Branch("pT_gen",&pt_gen,"pT_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("eta_gen",&eta_gen,"eta_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("y_gen",&y_gen,"y_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("status_gen",&status_gen,"status_gen[candSize_gen]/I");
        VertexCompositeNtuple->Branch("MotherID_gen",&idmom,"MotherID_gen[candSize_gen]/I");

        if(decayInGen_)
        {

            VertexCompositeNtuple->Branch("DauID1_gen",&iddau1,"DauID1_gen[candSize_gen]/I");
            VertexCompositeNtuple->Branch("DauID2_gen",&iddau2,"DauID2_gen[candSize_gen]/I");
            VertexCompositeNtuple->Branch("DauID3_gen",&iddau3,"DauID3_gen[candSize_gen]/I");
        }
    }
}

int VertexCompositeTreeProducer::
muAssocToTrack( const reco::TrackRef& trackref,
                const edm::Handle<reco::MuonCollection>& muonh) const {
  auto muon = std::find_if(muonh->cbegin(),muonh->cend(),
                           [&](const reco::Muon& m) {
                             return ( m.track().isNonnull() &&
                                      m.track() == trackref    );
                           });
  return ( muon != muonh->cend() ? std::distance(muonh->cbegin(),muon) : -1 );
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
VertexCompositeTreeProducer::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexCompositeTreeProducer);
