// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
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



#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

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
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define PI 3.1416

using namespace std;

class VertexCompositeSelector : public edm::EDProducer {
public:
  explicit VertexCompositeSelector(const edm::ParameterSet&);
  ~VertexCompositeSelector();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;

  double GetMVACut(double y, double pt);
  int muAssocToTrack( const reco::TrackRef& trackref, const edm::Handle<reco::MuonCollection>& muonh) const;

  // ----------member data ---------------------------
    
    //options
    bool doGenMatching_;
    bool hasSwap_;
    bool decayInGen_;
    bool twoLayerDecay_;
    bool threeProngDecay_;
    bool doMuon_;
    bool selectGenMatch_;
    bool selectGenUnMatch_;
    bool selectGenMatchSwap_;
    bool selectGenMatchUnSwap_;

    int PID_;
    int PID_dau1_;
    int PID_dau2_;
    int PID_dau3_;
    
    //cut variables
    double multMax_;
    double multMin_;
    double deltaR_; //deltaR for Gen matching
    bool   trkHighPurity_;
    double trkPMin_;
    double trkPtMin_;
    double trkEtaMax_;
    double trkPSumMin_;
    double trkPtSumMin_;
    double trkPtAsymMin_;
    double trkEtaDiffMax_;
    double trkPtErrMax_;
    int    trkNHitMin_;
    double candpTMin_;
    double candpTMax_;
    double candYMin_;
    double candYMax_;
    double cand3DDecayLengthSigMin_;
    double cand2DDecayLengthSigMin_;
    double cand3DPointingAngleMax_;
    double cand2DPointingAngleMax_;
    double cand3DDCAMin_;
    double cand3DDCAMax_;
    double cand2DDCAMin_;
    double cand2DDCAMax_;
    double candVtxProbMin_;

    //tree branches
    //event info
    int centrality;
    int Ntrkoffline;
    float bestvx;
    float bestvy;
    float bestvz;
    
    //Composite candidate info
    float mva;
    float pt;
    float eta;
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
    float dzos1;
    float dzos2;
    float dzos3;
    float dxyos1;
    float dxyos2;
    float dxyos3;
    float nhit1;
    float nhit2;
    float nhit3;
    bool trkquality1;
    bool trkquality2;
    bool trkquality3;
    float pt1;
    float pt2;
    float pt3;
    float ptErr1;
    float ptErr2;
    float ptErr3;
    float p1;
    float p2;
    float p3;
    float eta1;
    float eta2;
    float eta3;
    float phi1;
    float phi2;
    float phi3;
//    int charge1;
//    int charge2;
//    int charge3;
    float H2dedx1;
    float H2dedx2;
    float H2dedx3;
    float T4dedx1;
    float T4dedx2;
    float T4dedx3;
    float trkChi1;
    float trkChi2;
    float trkChi3;
    bool  isPionD1;
    bool  isPionD2;
    bool  isPionD3;
    bool  isKaonD1;
    bool  isKaonD2;
    bool  isKaonD3;
    
    //grand-dau info
    float grand_dzos1;
    float grand_dzos2;
    float grand_dxyos1;
    float grand_dxyos2;
    float grand_nhit1;
    float grand_nhit2;
    bool grand_trkquality1;
    bool grand_trkquality2;
    float grand_pt1;
    float grand_pt2;
    float grand_ptErr1;
    float grand_ptErr2;
    float grand_p1;
    float grand_p2;
    float grand_eta1;
    float grand_eta2;
//    int grand_charge1;
//    int grand_charge2;
    float grand_H2dedx1;
    float grand_H2dedx2;
    float grand_T4dedx1;
    float grand_T4dedx2;
    float grand_trkChi1;
    float grand_trkChi2;
    
    //dau muon info
    float nmatchedst1;
    float nmatchedch1;
    float ntrackerlayer1;
    float npixellayer1;
    float matchedenergy1;
    float nmatchedst2;
    float nmatchedch2;
    float ntrackerlayer2;
    float npixellayer2;
    float matchedenergy2;
    float dx1_seg_;
    float dy1_seg_;
    float dxSig1_seg_;
    float dySig1_seg_;
    float ddxdz1_seg_;
    float ddydz1_seg_;
    float ddxdzSig1_seg_;
    float ddydzSig1_seg_;
    float dx2_seg_;
    float dy2_seg_;
    float dxSig2_seg_;
    float dySig2_seg_;
    float ddxdz2_seg_;
    float ddydz2_seg_;
    float ddxdzSig2_seg_;
    float ddydzSig2_seg_;
    
    //vector for gen match
    vector< vector<double> > *pVect;
    vector<double> *Dvector1;
    vector<double> *Dvector2;
    vector<double> *Dvector3;
    vector<int> *pVectIDmom;
    
    int  selectFlavor_;
    bool usePID_;
    bool useAnyMVA_;
    bool useExistingMVA_;

    std::string mvaType_;
    std::string forestLabel_;
    GBRForest * forest_;
    bool useForestFromDB_;
    std::vector<float> mvaVals_;
    std::string dbFileName_;

    TF2* func_mva;
    std::vector<double> mvaCuts_;

    TH2D* hist_bdtcut;

    float mvaMin_;
    float mvaMax_;

    int   centMin_;
    int   centMax_;
    bool isCentrality_;

    edm::Handle<int> cbin_;

    //tokens
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
    edm::EDGetTokenT<MVACollection> MVAValues_Token_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
    edm::EDGetTokenT<reco::MuonCollection> tok_muon_;
    edm::EDGetTokenT<int> tok_centBinLabel_;
    edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

    std::string v0IDName_;

    reco::VertexCompositeCandidateCollection theVertexComps;
    MVACollection theMVANew;
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

VertexCompositeSelector::VertexCompositeSelector(const edm::ParameterSet& iConfig)
{
    //options
    twoLayerDecay_ = iConfig.getUntrackedParameter<bool>("twoLayerDecay");
    threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
    doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
    hasSwap_ = iConfig.getUntrackedParameter<bool>("hasSwap");
    decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
    doMuon_ = iConfig.getUntrackedParameter<bool>("doMuon");
    selectGenMatch_ = iConfig.getUntrackedParameter<bool>("selectGenMatch");
    selectGenUnMatch_ = iConfig.getUntrackedParameter<bool>("selectGenUnMatch");
    selectGenMatchSwap_ = iConfig.getUntrackedParameter<bool>("selectGenMatchSwap");
    selectGenMatchUnSwap_ = iConfig.getUntrackedParameter<bool>("selectGenMatchUnSwap");

    PID_ = iConfig.getUntrackedParameter<int>("PID");
    PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
    PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
    if(threeProngDecay_) PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");
    
    //cut variables
    centMin_ = iConfig.getUntrackedParameter<int>("centMin", 0);
    centMax_ = iConfig.getUntrackedParameter<int>("centMax", 10000);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
    deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);
    mvaMax_ = iConfig.getUntrackedParameter<double>("mvaMax", 999.9);
    mvaMin_ = iConfig.getUntrackedParameter<double>("mvaMin", -999.9);

    trkHighPurity_ = iConfig.getUntrackedParameter<bool>("trkHighPurity");
    trkPMin_ = iConfig.getUntrackedParameter<double>("trkPMin");
    trkPtMin_ = iConfig.getUntrackedParameter<double>("trkPtMin");
    trkEtaMax_ = iConfig.getUntrackedParameter<double>("trkEtaMax");
    trkPSumMin_ = iConfig.getUntrackedParameter<double>("trkPSumMin");
    trkPtSumMin_ = iConfig.getUntrackedParameter<double>("trkPtSumMin");
    trkPtAsymMin_ = iConfig.getUntrackedParameter<double>("trkPtAsymMin");
    trkEtaDiffMax_ = iConfig.getUntrackedParameter<double>("trkEtaDiffMax");
    trkPtErrMax_ = iConfig.getUntrackedParameter<double>("trkPtErrMax");
    trkNHitMin_ = iConfig.getUntrackedParameter<int>("trkNHitMin");
    candpTMin_ = iConfig.getUntrackedParameter<double>("candpTMin");
    candpTMax_ = iConfig.getUntrackedParameter<double>("candpTMax");
    candYMin_ = iConfig.getUntrackedParameter<double>("candYMin");
    candYMax_ = iConfig.getUntrackedParameter<double>("candYMax");
    cand3DDecayLengthSigMin_ = iConfig.getUntrackedParameter<double>("cand3DDecayLengthSigMin");
    cand2DDecayLengthSigMin_ = iConfig.getUntrackedParameter<double>("cand2DDecayLengthSigMin");
    cand3DPointingAngleMax_ = iConfig.getUntrackedParameter<double>("cand3DPointingAngleMax");
    cand2DPointingAngleMax_ = iConfig.getUntrackedParameter<double>("cand2DPointingAngleMax");
    cand3DDCAMin_ = iConfig.getUntrackedParameter<double>("cand3DDCAMin");
    cand3DDCAMax_ = iConfig.getUntrackedParameter<double>("cand3DDCAMax");
    cand2DDCAMin_ = iConfig.getUntrackedParameter<double>("cand2DDCAMin");
    cand2DDCAMax_ = iConfig.getUntrackedParameter<double>("cand2DDCAMax");
    candVtxProbMin_ = iConfig.getUntrackedParameter<double>("candVtxProbMin");

    //input tokens
    tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
    recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"));
    tok_muon_ = consumes<reco::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection"));
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
    tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

    usePID_ = false;
    selectFlavor_ = 0;
    if(iConfig.exists("usePID")) usePID_ = iConfig.getParameter<bool>("usePID");
    if(iConfig.exists("useFlavor")) selectFlavor_ = iConfig.getUntrackedParameter<int>("selectFlavor");
 
    // Loading TMVA
    useAnyMVA_ = false;
    useExistingMVA_ = false;

    forestLabel_ = "D0InpPb";
    std::string type = "BDT";
    useForestFromDB_ = true;
    dbFileName_ = "";

    forest_ = nullptr;

    isCentrality_ = false;
    if(iConfig.exists("isCentrality")) isCentrality_ = iConfig.getParameter<bool>("isCentrality");
    if(isCentrality_)
    {
      tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
      tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
    }

    if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
    if(iConfig.exists("useExistingMVA")) useExistingMVA_ = iConfig.getParameter<bool>("useExistingMVA");

    if(useAnyMVA_){
      if(useExistingMVA_ && iConfig.exists("MVACollection")){ 
        MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
      }
      else{
        if(iConfig.exists("mvaType"))type = iConfig.getParameter<std::string>("mvaType");
        if(iConfig.exists("GBRForestLabel"))forestLabel_ = iConfig.getParameter<std::string>("GBRForestLabel");
        if(iConfig.exists("GBRForestFileName")){
          dbFileName_ = iConfig.getParameter<std::string>("GBRForestFileName");
          useForestFromDB_ = false;
        }
      }
   
      mvaCuts_ = iConfig.getParameter< std::vector<double> >("mvaCuts");

      TString bdtcut_filename;
      if(iConfig.exists("BDTCutFileName")) bdtcut_filename = iConfig.getParameter<string>("BDTCutFileName"); 
      hist_bdtcut = 0;
      if(!bdtcut_filename.IsNull()) 
      {
        edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeAnalyzer/data/%s",bdtcut_filename.Data()));
        TFile ff(fip.fullPath().c_str(),"READ");
        hist_bdtcut = (TH2D*)ff.Get("hist_bdtcut");
        ff.Close();
      }

      func_mva = new TF2("func_mva","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)*(1+[4]*y+[5]*y*y+[6]*y*y*y+[7]*y*y*y*y)",0,5.0,0,100);
      func_mva->SetParameters(mvaCuts_[0],mvaCuts_[1],mvaCuts_[2],mvaCuts_[3],mvaCuts_[4],mvaCuts_[5],mvaCuts_[6],mvaCuts_[7]);
    }

    if(!useForestFromDB_){
      edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeAnalyzer/data/%s",dbFileName_.c_str()));
      TFile gbrfile(fip.fullPath().c_str(),"READ");
      forest_ = (GBRForest*)gbrfile.Get(forestLabel_.c_str());
      gbrfile.Close();
    }

    mvaType_ = type;

    v0IDName_ = (iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection")).instance();

    produces< reco::VertexCompositeCandidateCollection >(v0IDName_);
    produces<MVACollection>(Form("MVAValuesNew%s",v0IDName_.c_str()));

    isPionD1 = true;
    isPionD2 = true;
    isKaonD1 = false;
    isKaonD2 = false;
}


VertexCompositeSelector::~VertexCompositeSelector()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VertexCompositeSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using std::vector;
    using namespace edm;

    fillRECO(iEvent,iSetup);

//    std::make_unique< reco::VertexCompositeCandidateCollection > theNewV0Cands( new reco::VertexCompositeCandidateCollection );
    auto theNewV0Cands = std::make_unique<reco::VertexCompositeCandidateCollection>();

    theNewV0Cands->reserve( theVertexComps.size() );

    std::copy( theVertexComps.begin(),
               theVertexComps.end(),
               std::back_inserter(*theNewV0Cands) );

    iEvent.put(std::move(theNewV0Cands), v0IDName_);
    theVertexComps.clear();

    if(useAnyMVA_)
    {
      auto mvas = std::make_unique<MVACollection>(theMVANew.begin(),theMVANew.end());
      iEvent.put(std::move(mvas), Form("MVAValuesNew%s",v0IDName_.c_str()));
      theMVANew.clear();
    }
}

void
VertexCompositeSelector::fillRECO(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //get collections
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
    iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,v0candidates);
    const reco::VertexCompositeCandidateCollection * v0candidates_ = v0candidates.product();
    
    edm::Handle<MVACollection> mvavalues;
    if(useAnyMVA_ && useExistingMVA_)
    {
      iEvent.getByToken(MVAValues_Token_,mvavalues);
      assert( (*mvavalues).size() == v0candidates->size() );
    }

    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
    if(usePID_) iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
    
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
    if(usePID_) iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
    
    centrality=-1;
    if(isCentrality_)
    {
      edm::Handle<reco::Centrality> cent;
      iEvent.getByToken(tok_centSrc_, cent);

      iEvent.getByToken(tok_centBinLabel_,cbin_);
      centrality = *cbin_;

//      HFsumET = cent->EtHFtowerSum();
//      Npixel = cent->multiplicityPixel();
//      int ntrk = cent->Ntracks();
    }
    if(centrality!=-1 && (centrality >= centMax_ || centrality < centMin_)) return;

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

      if(Ntrkoffline >= multMax_ || Ntrkoffline < multMin_) return;
    }

    //Gen info for matching
    if(doGenMatching_)
    {
        pVect = new vector< vector<double>>;
        pVectIDmom = new vector<int>;
        
        edm::Handle<reco::GenParticleCollection> genpars;
        iEvent.getByToken(tok_genParticle_,genpars);
        
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
    for(unsigned it=0; it<v0candidates_->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();

        eta = trk.eta();
        y = trk.rapidity();
        pt = trk.pt();
        flavor = trk.pdgId()/421;

        // select particle vs antiparticle  
        if(usePID_ && selectFlavor_ && (int)flavor!=selectFlavor_) continue; 

        // select on pT and y
        if(pt<candpTMin_ || pt>candpTMax_) continue;
        if(y<candYMin_ || y>candYMax_) continue;

        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        mass = trk.mass();
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        const reco::Candidate * d3 = 0;
        if(threeProngDecay_) d3 = trk.daughter(2);        

        //Gen match
        if(doGenMatching_)
        {
            matchGEN = false;
            int nGenDau = (int)pVect->size();
            isSwap = false;
            idmom_reco = -77;

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

                    matchGEN = true; //matched gen
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

                    matchGEN = true; //matched gen
                  }

                  if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01) isSwap = true;

                  idmom_reco = pVectIDmom->at(i/2);
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

                    matchGEN = true; //matched gen
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

                    matchGEN = true; //matched gen
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

                    matchGEN = true; //matched gen
                  }

                  if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01 || abs(d3massGEN - d3mass)>0.01) isSwap = true;

                  idmom_reco = pVectIDmom->at(i/3);
                }

            }

            if(selectGenMatch_ && !matchGEN) continue;
            if(selectGenUnMatch_ && matchGEN) continue;
            if(selectGenMatchSwap_)
            {
              if(!matchGEN) continue;
              else if(!isSwap) continue;
            }
            if(selectGenMatchUnSwap_)
            {
              if(!matchGEN) continue;
              else if(isSwap) continue;
            }
        }

/*
        double pxd1 = d1->px();
        double pyd1 = d1->py();
        double pzd1 = d1->pz();
        double pxd2 = d2->px();
        double pyd2 = d2->py();
        double pzd2 = d2->pz();

        TVector3 dauvec1(pxd1,pyd1,pzd1);
        TVector3 dauvec2(pxd2,pyd2,pzd2);
*/
        //pt
        pt1 = d1->pt();
        pt2 = d2->pt();

        if(pt1 < trkPtMin_ || pt2 < trkPtMin_) continue;
        if((pt1+pt2) < trkPtSumMin_) continue;
                
        if(pt2/pt1 < trkPtAsymMin_ || pt1/pt2 < trkPtAsymMin_) continue;

        //momentum
        p1 = d1->p();
        p2 = d2->p();
        
        if(p1 < trkPMin_ || p2 < trkPMin_) continue;
        if((p1+p2) < trkPSumMin_) continue;

        //eta
        eta1 = d1->eta();
        eta2 = d2->eta();
        
        if(fabs(eta1) > trkEtaMax_ || fabs(eta2) > trkEtaMax_) continue;
        if(fabs(eta1-eta2) > trkEtaDiffMax_) continue;

        //phi
        phi1 = d1->phi();
        phi2 = d2->phi();
        
        //charge
//        charge1 = d1->charge();
//        charge2 = d2->charge();
        
        if(threeProngDecay_)
        {
          pt3 = d3->pt();
          if(pt3 < trkPtMin_) continue;
          p3 = d3->p();
          if(p3 < trkPMin_) continue;
          eta3 = d3->eta();
          if(fabs(eta3) > trkEtaMax_) continue;
        }

        //vtxChi2
        vtxChi2 = trk.vertexChi2();
        ndf = trk.vertexNdof();
        VtxProb = TMath::Prob(vtxChi2,ndf);

        if(VtxProb < candVtxProbMin_) continue;       

        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
        TVector3 secvec2D(px,py,0);
        
        agl = cos(secvec.Angle(ptosvec));
        agl_abs = secvec.Angle(ptosvec);
        if(agl_abs > cand3DPointingAngleMax_) continue;

        agl2D = cos(secvec2D.Angle(ptosvec2D));
        agl2D_abs = secvec2D.Angle(ptosvec2D);
        if(agl2D_abs > cand2DPointingAngleMax_) continue;
        
        //Decay length 3D
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        typedef ROOT::Math::SVector<double, 6> SVector6;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        dl = ROOT::Math::Mag(distanceVector);
        dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
        
        dlos = dl/dlerror;
        if(dlos < cand3DDecayLengthSigMin_ || dlos > 1000.) continue;
 
        //Decay length 2D
        SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
        SVector6 v2(trk.vertexCovariance(0,0), trk.vertexCovariance(0,1),trk.vertexCovariance(1,1),0,0,0);
        
        SMatrixSym3D sv1(v1);
        SMatrixSym3D sv2(v2);
        
        SMatrixSym3D totalCov2D = sv1 + sv2;
        SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
        
        double dl2D = ROOT::Math::Mag(distanceVector2D);
        double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D;
        
        dlos2D = dl2D/dl2Derror;
        if(dlos2D < cand2DDecayLengthSigMin_ || dlos2D > 1000.) continue;

        double dca3D = dl*sin(agl_abs);
        if(dca3D < cand3DDCAMin_ || dca3D > cand3DDCAMax_) continue;

        double dca2D = dl2D*sin(agl2D_abs);
        if(dca2D < cand2DDCAMin_ || dca2D > cand2DDCAMax_) continue;

        //trk info
        auto dau1 = d1->get<reco::TrackRef>();
        if(!twoLayerDecay_)
        {
            //trk quality
            trkquality1 = dau1->quality(reco::TrackBase::highPurity);
            if(trkHighPurity_ && !trkquality1) continue;
            
            //trk dEdx
            H2dedx1 = -999.9;
            T4dedx1 = -999.9;
            if(usePID_)
            {
               if(dEdxHandle1.isValid())
               {
                  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                  H2dedx1 = dEdxTrack[dau1].dEdx();
                  if(H2dedx1 > (2.8/pow(pt1*cosh(eta1),0.4)+0.2) && H2dedx1 < (2.8/pow(pt1*cosh(eta1),0.9)+1.8) && H2dedx1> (2.8/pow(0.75,0.4)+0.2)) { isKaonD1 = true; isPionD1 = false; }
                  if((H2dedx1 < (2.8/pow(pt1*cosh(eta1),0.4)+0.2) || H2dedx1< (2.8/pow(0.75,0.4)+0.2)) && H2dedx1>0) { isPionD1 = true; isKaonD1 = false;}
               }
            
               if(dEdxHandle2.isValid())
               {
                  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                  T4dedx1 = dEdxTrack[dau1].dEdx();
               }
            }
            
            //track Chi2
            trkChi1 = dau1->normalizedChi2();
            
            //track pT error
            ptErr1 = dau1->ptError();
            if(ptErr1/dau1->pt() > trkPtErrMax_) continue;
    
            //vertexCovariance 00-xError 11-y 22-z
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
           
            //trkNHits
            nhit1 = dau1->numberOfValidHits();
            if(nhit1 < trkNHitMin_) continue;
            
            //DCA
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzbest1 = dau1->dz(bestvtx);
            double dxybest1 = dau1->dxy(bestvtx);
            double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
            double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
            
            dzos1 = dzbest1/dzerror1;
            dxyos1 = dxybest1/dxyerror1;
        }
        
        auto dau2 = d2->get<reco::TrackRef>();
        
        trkquality2 = dau2->quality(reco::TrackBase::highPurity);
        if(trkHighPurity_ && !trkquality2) continue;

        //trk dEdx
        H2dedx2 = -999.9;
        T4dedx2 = -999.9;
        
        if(usePID_)
        {
          if(dEdxHandle1.isValid())
          {
             const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
             H2dedx2 = dEdxTrack[dau2].dEdx();
        
             if(H2dedx2 > (2.8/pow(pt2*cosh(eta2),0.4)+0.2) && H2dedx2 < (2.8/pow(pt2*cosh(eta2),0.9)+1.8) && H2dedx2> (2.8/pow(0.75,0.4)+0.2)) { isKaonD2 = true; isPionD2 = false; }
             if((H2dedx2 < (2.8/pow(pt2*cosh(eta2),0.4)+0.2) || H2dedx2< (2.8/pow(0.75,0.4)+0.2)) && H2dedx2>0) { isPionD2 = true; isKaonD2 = false; }
          }

          if(dEdxHandle2.isValid()){
             const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
             T4dedx2 = dEdxTrack[dau2].dEdx();
          }

          if(flavor>0 && (!isPionD1 || !isKaonD2)) continue;
          if(flavor<0 && (!isPionD2 || !isKaonD1)) continue;
        }

        //track Chi2
        trkChi2 = dau2->normalizedChi2();
        
        //track pT error
        ptErr2 = dau2->ptError();
        if(ptErr2/dau2->pt() > trkPtErrMax_) continue;

        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        
        //trkNHits
        nhit2 = dau2->numberOfValidHits();
        if(nhit2 < trkNHitMin_) continue;
        
        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        dzos2 = dzbest2/dzerror2;
        dxyos2 = dxybest2/dxyerror2;
        
        if(threeProngDecay_)
        {
          auto dau3 = d3->get<reco::TrackRef>();

          trkquality3 = dau3->quality(reco::TrackBase::highPurity);
          if(trkHighPurity_ && !trkquality3) continue;
          trkChi3 = dau3->normalizedChi2();
          ptErr3 = dau3->ptError();
          if(ptErr3/dau3->pt() > trkPtErrMax_) continue;
          nhit3 = dau3->numberOfValidHits();
          if(nhit3 < trkNHitMin_) continue;

          double dzbest3 = dau3->dz(bestvtx);
          double dxybest3 = dau3->dxy(bestvtx);
          double dzerror3 = sqrt(dau3->dzError()*dau3->dzError()+bestvzError*bestvzError);
          double dxyerror3 = sqrt(dau3->d0Error()*dau3->d0Error()+bestvxError*bestvyError);
          dzos3 = dzbest3/dzerror3;
          dxyos3 = dxybest3/dxyerror3;
        }

        if(doMuon_)
        {
            edm::Handle<reco::MuonCollection> theMuonHandle;
            iEvent.getByToken(tok_muon_, theMuonHandle);
            
            nmatchedch1 = -1;
            nmatchedst1 = -1;
            matchedenergy1 = -1;
            nmatchedch2 = -1;
            nmatchedst2 = -1;
            matchedenergy2 = -1;
            
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
            

            const int muId1 = muAssocToTrack( dau1, theMuonHandle );
            if( muId1 != -1 )
            {
              const reco::Muon& cand = (*theMuonHandle)[muId1];
              nmatchedch1 = cand.numberOfMatches();
              nmatchedst1 = cand.numberOfMatchedStations();
                    
              reco::MuonEnergy muenergy = cand.calEnergy();
              matchedenergy1 = muenergy.hadMax;
                    
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
                        
                dx1_seg_=dx_seg;
                dy1_seg_=dy_seg;
                dxSig1_seg_=dxSig_seg;
                dySig1_seg_=dySig_seg;
                ddxdz1_seg_=ddxdz_seg;
                ddydz1_seg_=ddydz_seg;
                ddxdzSig1_seg_=ddxdzSig_seg;
                ddydzSig1_seg_=ddydzSig_seg;
              }
            }
                

            const int muId2 = muAssocToTrack( dau2, theMuonHandle );
            if( muId2 != -1 )
            {
              const reco::Muon& cand = (*theMuonHandle)[muId2];                

              nmatchedch2 = cand.numberOfMatches();
              nmatchedst2 = cand.numberOfMatchedStations();
                    
              reco::MuonEnergy muenergy = cand.calEnergy();
              matchedenergy2 = muenergy.hadMax;
                    
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
                        
                dx2_seg_=dx_seg;
                dy2_seg_=dy_seg;
                dxSig2_seg_=dxSig_seg;
                dySig2_seg_=dySig_seg;
                ddxdz2_seg_=ddxdz_seg;
                ddydz2_seg_=ddydz_seg;
                ddxdzSig2_seg_=ddxdzSig_seg;
                ddydzSig2_seg_=ddydzSig_seg;
              }
            }
        }
        
        if(twoLayerDecay_)
        {
            grand_mass = d1->mass();
            
            const reco::Candidate * gd1 = d1->daughter(0);
            const reco::Candidate * gd2 = d1->daughter(1);
/*            
            double gpxd1 = gd1->px();
            double gpyd1 = gd1->py();
            double gpzd1 = gd1->pz();
            double gpxd2 = gd2->px();
            double gpyd2 = gd2->py();
            double gpzd2 = gd2->pz();
            
            TVector3 gdauvec1(gpxd1,gpyd1,gpzd1);
            TVector3 gdauvec2(gpxd2,gpyd2,gpzd2);
*/
            auto gdau1 = gd1->get<reco::TrackRef>();
            auto gdau2 = gd2->get<reco::TrackRef>();
            
            //trk quality
            
            grand_trkquality1 = gdau1->quality(reco::TrackBase::highPurity);
            grand_trkquality2 = gdau2->quality(reco::TrackBase::highPurity);
            
            //trk dEdx
            grand_H2dedx1 = -999.9;
            grand_H2dedx2 = -999.9;
            
            if(dEdxHandle1.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                grand_H2dedx1 = dEdxTrack[gdau1].dEdx();
                grand_H2dedx2 = dEdxTrack[gdau2].dEdx();
            }
            
            grand_T4dedx1 = -999.9;
            grand_T4dedx2 = -999.9;
            
            if(dEdxHandle2.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                grand_T4dedx1 = dEdxTrack[gdau1].dEdx();
                grand_T4dedx2 = dEdxTrack[gdau2].dEdx();
            }
            
            //track pt
            grand_pt1 = gd1->pt();
            grand_pt2 = gd2->pt();
            
            //track momentum
            grand_p1 = gd1->p();
            grand_p2 = gd2->p();
            
            //track eta
            grand_eta1 = gd1->eta();
            grand_eta2 = gd2->eta();
            
            //track charge
//            grand_charge1 = gd1->charge();
//            grand_charge2 = gd2->charge();
            
            //track Chi2
//            grand_trkChi1 = gdau1->normalizedChi2();
//            grand_trkChi2 = gdau2->normalizedChi2();
            
            //track pT error
            grand_ptErr1 = gdau1->ptError();
            grand_ptErr2 = gdau2->ptError();
            
            //vertexCovariance 00-xError 11-y 22-z
            secvz = d1->vz(); secvx = d1->vx(); secvy = d1->vy();
            
            //trkNHits
            grand_nhit1 = gdau1->numberOfValidHits();
            grand_nhit2 = gdau2->numberOfValidHits();
            
            //DCA
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double gdzbest1 = gdau1->dz(bestvtx);
            double gdxybest1 = gdau1->dxy(bestvtx);
            double gdzerror1 = sqrt(gdau1->dzError()*gdau1->dzError()+bestvzError*bestvzError);
            double gdxyerror1 = sqrt(gdau1->d0Error()*gdau1->d0Error()+bestvxError*bestvyError);
            
            grand_dzos1 = gdzbest1/gdzerror1;
            grand_dxyos1 = gdxybest1/gdxyerror1;
            
            double gdzbest2 = gdau2->dz(bestvtx);
            double gdxybest2 = gdau2->dxy(bestvtx);
            double gdzerror2 = sqrt(gdau2->dzError()*gdau2->dzError()+bestvzError*bestvzError);
            double gdxyerror2 = sqrt(gdau2->d0Error()*gdau2->d0Error()+bestvxError*bestvyError);
            
            grand_dzos2 = gdzbest2/gdzerror2;
            grand_dxyos2 = gdxybest2/gdxyerror2;
            
            //vtxChi2
            grand_vtxChi2 = d1->vertexChi2();
            grand_ndf = d1->vertexNdof();
            grand_VtxProb = TMath::Prob(grand_vtxChi2,grand_ndf);
            
            //PAngle
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(d1->px(),d1->py(),d1->pz());
            
            TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
            TVector3 secvec2D(d1->px(),d1->py(),0);
            
            grand_agl = cos(secvec.Angle(ptosvec));
            grand_agl_abs = secvec.Angle(ptosvec);
            
            grand_agl2D = cos(secvec2D.Angle(ptosvec2D));
            grand_agl2D_abs = secvec2D.Angle(ptosvec2D);
            
            //Decay length 3D
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            typedef ROOT::Math::SVector<double, 6> SVector6;
            
            SMatrixSym3D totalCov = vtx.covariance() + d1->vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            grand_dl = ROOT::Math::Mag(distanceVector);
            grand_dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/grand_dl;
            
            grand_dlos = grand_dl/grand_dlerror;
            
            //Decay length 2D
            SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
            SVector6 v2(d1->vertexCovariance(0,0), d1->vertexCovariance(0,1),d1->vertexCovariance(1,1),0,0,0);
            
            SMatrixSym3D sv1(v1);
            SMatrixSym3D sv2(v2);
            
            SMatrixSym3D totalCov2D = sv1 + sv2;
            SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
            
            double gdl2D = ROOT::Math::Mag(distanceVector2D);
            double gdl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/gdl2D;
            
            grand_dlos2D = gdl2D/gdl2Derror;
        }

        // select MVA value
        mva=0;
        if(useAnyMVA_ && useExistingMVA_)
        {
          mva = (*mvavalues)[it];
          if(mva < mvaMin_ || mva > mvaMax_) continue;

          if(mva<GetMVACut(y,pt)) continue;

          theVertexComps.push_back( trk );
          theMVANew.push_back( mva );
          continue;
        }
        else if(useAnyMVA_ && !useExistingMVA_)
        {
          float gbrVals_[50]={0};
          if(forestLabel_ == "D0InpPb" || forestLabel_ == "D0Inpp" || forestLabel_ == "D0InPbPb")
          { 
/*
            gbrVals_[0] = pt;
            gbrVals_[1] = y;
            gbrVals_[2] = VtxProb;
            gbrVals_[3] = dlos;
            gbrVals_[4] = dlos2D;
            gbrVals_[5] = dl;
            gbrVals_[6] = agl_abs;
            gbrVals_[7] = agl2D_abs;
            gbrVals_[8] = dzos1;
            gbrVals_[9] = dzos2;
            gbrVals_[10] = dxyos1;
            gbrVals_[11] = dxyos2;
            gbrVals_[12] = pt1;
            gbrVals_[13] = pt2;
            gbrVals_[14] = eta1;
            gbrVals_[15] = eta2;
            gbrVals_[16] = nhit1;
            gbrVals_[17] = nhit2;
            gbrVals_[18] = ptErr1;
            gbrVals_[19] = ptErr2;
*/
            gbrVals_[0] = pt;
            gbrVals_[1] = y;
            gbrVals_[2] = VtxProb;
            gbrVals_[3] = dlos;
            gbrVals_[4] = dl;
            gbrVals_[5] = agl_abs;
            gbrVals_[6] = dzos1;
            gbrVals_[7] = dzos2;
            gbrVals_[8] = dxyos1;
            gbrVals_[9] = dxyos2;
            gbrVals_[10] = pt1;
            gbrVals_[11] = pt2;
            gbrVals_[12] = eta1;
            gbrVals_[13] = eta2;
            gbrVals_[14] = nhit1;
            gbrVals_[15] = nhit2;
            gbrVals_[16] = ptErr1;
            gbrVals_[17] = ptErr2;
            gbrVals_[18] = dlos2D;
            gbrVals_[19] = agl2D_abs;
          }

          if(forestLabel_ == "DsInpPb" || forestLabel_ == "DsInpp" || forestLabel_ == "DsInPbPb")
          {
            gbrVals_[0] = pt;
            gbrVals_[1] = y;
            gbrVals_[2] = VtxProb;
            gbrVals_[3] = dlos;
            gbrVals_[4] = dlos2D;
            gbrVals_[5] = dl;
            gbrVals_[6] = agl_abs;
            gbrVals_[7] = agl2D_abs;
            gbrVals_[8] = dzos2;
            gbrVals_[9] = dxyos2;
            gbrVals_[10] = pt1;
            gbrVals_[11] = pt2;
            gbrVals_[12] = eta1;
            gbrVals_[13] = eta2;
            gbrVals_[14] = nhit2;
            gbrVals_[15] = ptErr2;
            gbrVals_[16] = H2dedx2;
          }

          if(forestLabel_ == "DiMuInpPb" || forestLabel_ == "DiMuInpp" || forestLabel_ == "DiMuInPbPb")
          {
            gbrVals_[0] = pt;
            gbrVals_[1] = y;
            gbrVals_[2] = VtxProb;
            gbrVals_[3] = dlos;
            gbrVals_[4] = dlos2D;
            gbrVals_[5] = dl;
            gbrVals_[6] = agl_abs;
            gbrVals_[7] = agl2D_abs;
            gbrVals_[8] = dzos1;
            gbrVals_[9] = dzos2;
            gbrVals_[10] = dxyos1;
            gbrVals_[11] = dxyos2;
            gbrVals_[12] = nhit1;
            gbrVals_[13] = nhit2;
            gbrVals_[14] = nmatchedch1;
            gbrVals_[15] = nmatchedst1;
            gbrVals_[16] = matchedenergy1;
            gbrVals_[17] = nmatchedch2;
            gbrVals_[18] = nmatchedst2;
            gbrVals_[19] = matchedenergy2;
            gbrVals_[20] = dxSig1_seg_;
            gbrVals_[21] = dySig1_seg_;
            gbrVals_[22] = ddxdzSig1_seg_;
            gbrVals_[23] = ddydzSig1_seg_;
            gbrVals_[24] = dxSig2_seg_;
            gbrVals_[25] = dySig2_seg_;
            gbrVals_[26] = ddxdzSig2_seg_;
            gbrVals_[27] = ddydzSig2_seg_;
            gbrVals_[28] = pt1;
            gbrVals_[29] = pt2;
            gbrVals_[30] = eta1;
            gbrVals_[31] = eta2;
          }

          GBRForest const * forest = forest_;
          if(useForestFromDB_){
            edm::ESHandle<GBRForest> forestHandle;
            iSetup.get<GBRWrapperRcd>().get(forestLabel_,forestHandle);
            forest = forestHandle.product();
          }

          auto gbrVal = forest->GetClassifier(gbrVals_);

          if(gbrVal < mvaMin_ || gbrVal > mvaMax_) continue;
          if(gbrVal < GetMVACut(y,pt)) continue;

          theMVANew.push_back( gbrVal );
        } 
        theVertexComps.push_back( trk );
    }
}

double
VertexCompositeSelector::GetMVACut(double y, double pt)
{
  double mvacut = -1.0;
  if(fabs(y)>2.4) return 1.0;

  //temporary
  if(pt<4) return 0.45;
  else if(pt>4 && pt<6) return 0.25;
  else if(pt>6 && pt<8) return -0.2;
  else return -1.0;

  if(!hist_bdtcut) return mvacut;

  mvacut = hist_bdtcut->GetBinContent(hist_bdtcut->GetXaxis()->FindBin(y),hist_bdtcut->GetYaxis()->FindBin(pt));
  if(pt>7.4) mvacut = hist_bdtcut->GetBinContent(hist_bdtcut->GetXaxis()->FindBin(y),hist_bdtcut->GetYaxis()->FindBin(7.4));
  if(pt<1.37) mvacut = hist_bdtcut->GetBinContent(hist_bdtcut->GetXaxis()->FindBin(y),hist_bdtcut->GetYaxis()->FindBin(1.37));

  return mvacut;
}

int VertexCompositeSelector::
muAssocToTrack( const reco::TrackRef& trackref,
                const edm::Handle<reco::MuonCollection>& muonh) const {
  auto muon = std::find_if(muonh->cbegin(),muonh->cend(),
                           [&](const reco::Muon& m) {
                             return ( m.track().isNonnull() &&
                                      m.track() == trackref    );
                           });
  return ( muon != muonh->cend() ? std::distance(muonh->cbegin(),muon) : -1 );
}

// ------------ method called once each job just before starting event
//loop  ------------
void
VertexCompositeSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
VertexCompositeSelector::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_FWK_MODULE(VertexCompositeSelector);
