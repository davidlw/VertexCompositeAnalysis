// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

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

  // ----------member data ---------------------------
    
    //options
    bool doGenMatching_;
    bool hasSwap_;
    bool decayInGen_;
    bool twoLayerDecay_;
    bool doMuon_;
    int PID_;
    int PID_dau1_;
    int PID_dau2_;
    
    //cut variables
    double multMax_;
    double multMin_;
    double deltaR_; //deltaR for Gen matching

    //tree branches
    //event info
    int Ntrkoffline;
    float bestvx;
    float bestvy;
    float bestvz;
    
    //Composite candidate info
    float mva;
    float pt;
    float eta;
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
    float dxyos1;
    float dxyos2;
    float nhit1;
    float nhit2;
    bool trkquality1;
    bool trkquality2;
    float pt1;
    float pt2;
    float ptErr1;
    float ptErr2;
    float p1;
    float p2;
    float eta1;
    float eta2;
    float phi1;
    float phi2;
    int charge1;
    int charge2;
    float H2dedx1;
    float H2dedx2;
    float T4dedx1;
    float T4dedx2;
    float trkChi1;
    float trkChi2;
    
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
    int grand_charge1;
    int grand_charge2;
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
    vector<int> *pVectIDmom;
    
    bool useAnyMVA_;
    float mvaMin_;
    float mvaMax_;

    //tokens
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
    edm::EDGetTokenT<MVACollection> MVAValues_Token_;

    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
    edm::EDGetTokenT<reco::MuonCollection> tok_muon_;

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
    doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
    hasSwap_ = iConfig.getUntrackedParameter<bool>("hasSwap");
    decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
    doMuon_ = iConfig.getUntrackedParameter<bool>("doMuon");
    PID_ = iConfig.getUntrackedParameter<int>("PID");
    PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
    PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
    
    useAnyMVA_ = iConfig.getUntrackedParameter<bool>("useAnyMVA");

    //cut variables
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 0.0);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 999.9);
    deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);
    mvaMax_ = iConfig.getUntrackedParameter<double>("mvaMax", 999.9);
    mvaMin_ = iConfig.getUntrackedParameter<double>("mvaMin", -999.9);

    //input tokens
    tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
    recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"));
    MVAValues_Token_ = consumes<MVACollection>(iConfig.getUntrackedParameter<edm::InputTag>("MVACollection"));
    tok_muon_ = consumes<reco::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection"));
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
    
    tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

    v0IDName_ = (iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection")).instance();

    produces< reco::VertexCompositeCandidateCollection >(v0IDName_);
    produces<MVACollection>("MVAValuesNew");
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

    std::auto_ptr< reco::VertexCompositeCandidateCollection > theNewV0Cands( new reco::VertexCompositeCandidateCollection );
    theNewV0Cands->reserve( theVertexComps.size() );

    std::copy( theVertexComps.begin(),
               theVertexComps.end(),
               std::back_inserter(*theNewV0Cands) );

    iEvent.put(theNewV0Cands, v0IDName_);
    theVertexComps.clear();

    if(useAnyMVA_)
    {
      auto mvas = std::make_unique<MVACollection>(theMVANew.begin(),theMVANew.end());
      iEvent.put(std::move(mvas), std::string("MVAValuesNew"));
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
    if(useAnyMVA_)
    {
      iEvent.getByToken(MVAValues_Token_,mvavalues);
      assert( (*mvavalues).size() == v0candidates->size() );
    }

    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
    iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
    
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
    iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
    
    //best vertex
    bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    //Ntrkoffline
    Ntrkoffline = 0;
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
            if(decayInGen_ && trk.numberOfDaughters()!=2) continue; //check 2-pron decay if target decays in Gen
            
            int idmom_tmp = -77;
            
            if(trk.numberOfMothers()!=0)
            {
                const reco::Candidate * mom = trk.mother();
                idmom_tmp = mom->pdgId();
            }
            
            const reco::Candidate * Dd1 = trk.daughter(0);
            const reco::Candidate * Dd2 = trk.daughter(1);
            
            if(!(fabs(Dd1->pdgId())==PID_dau1_ && fabs(Dd2->pdgId())==PID_dau2_) && !(fabs(Dd2->pdgId())==PID_dau1_ && fabs(Dd1->pdgId())==PID_dau2_)) continue; //check daughter id
                
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

        mva=0;
        if(useAnyMVA_) 
        {
          mva = (*mvavalues)[it];
          if(mva < mvaMin_ || mva > mvaMax_) continue;
          theVertexComps.push_back( trk );
          theMVANew.push_back( mva );
          continue;
        }

        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        mass = trk.mass();
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        
        if(doGenMatching_)
        {
            //Gen match
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
                //check dau2
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
                
                //check swap
                if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01) isSwap = true;
                
                //check prompt & record mom id
                idmom_reco = pVectIDmom->at(i/2);
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
        pt1 = d1->pt();
        pt2 = d2->pt();
        
        //momentum
        p1 = d1->p();
        p2 = d2->p();
        
        //eta
        eta1 = d1->eta();
        eta2 = d2->eta();
        
        //phi
        phi1 = d1->phi();
        phi2 = d2->phi();
        
        //charge
        charge1 = d1->charge();
        charge2 = d2->charge();
        
        //vtxChi2
        vtxChi2 = trk.vertexChi2();
        ndf = trk.vertexNdof();
        VtxProb = TMath::Prob(vtxChi2,ndf);
        
        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
        TVector3 secvec2D(px,py,0);
        
        agl = cos(secvec.Angle(ptosvec));
        agl_abs = secvec.Angle(ptosvec);
        
        agl2D = cos(secvec2D.Angle(ptosvec2D));
        agl2D_abs = secvec2D.Angle(ptosvec2D);
        
        //Decay length 3D
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        typedef ROOT::Math::SVector<double, 6> SVector6;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        dl = ROOT::Math::Mag(distanceVector);
        dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
        
        dlos = dl/dlerror;
        
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

        //trk info
        if(!twoLayerDecay_)
        {
            auto dau1 = d1->get<reco::TrackRef>();
            
            //trk quality
            trkquality1 = dau1->quality(reco::TrackBase::highPurity);
            
            //trk dEdx
            H2dedx1 = -999.9;
            
            if(dEdxHandle1.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                H2dedx1 = dEdxTrack[dau1].dEdx();
            }
            
            T4dedx1 = -999.9;
            
            if(dEdxHandle2.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                T4dedx1 = dEdxTrack[dau1].dEdx();
            }
            
            //track Chi2
            trkChi1 = dau1->normalizedChi2();
            
            //track pT error
            ptErr1 = dau1->ptError();
            
            //vertexCovariance 00-xError 11-y 22-z
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            
            //trkNHits
            nhit1 = dau1->numberOfValidHits();
            
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
        
        //trk quality
        trkquality2 = dau2->quality(reco::TrackBase::highPurity);
        
        //trk dEdx
        H2dedx2 = -999.9;
        
        if(dEdxHandle1.isValid()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
            H2dedx2 = dEdxTrack[dau2].dEdx();
        }
        
        T4dedx2 = -999.9;
        
        if(dEdxHandle2.isValid()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
            T4dedx2 = dEdxTrack[dau2].dEdx();
        }
        
        //track Chi2
        trkChi2 = dau2->normalizedChi2();
        
        //track pT error
        ptErr2 = dau2->ptError();
        
        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        
        //trkNHits
        nhit2 = dau2->numberOfValidHits();
        
        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        dzos2 = dzbest2/dzerror2;
        dxyos2 = dxybest2/dxyerror2;
        
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
            
            for( unsigned id = 0; id < theMuonHandle->size(); id++ )
            {
                const reco::Muon& cand = (*theMuonHandle)[id];
                reco::TrackRef trackReftmp = cand.track();
                
                if(trackReftmp.isNull()) continue;
                
                if(fabs(trackReftmp->pt()-pt1)<0.0001 && fabs(trackReftmp->eta()-eta1)<0.001 && fabs(trackReftmp->phi()-phi1)<0.001
                   && trackReftmp->numberOfValidHits() == nhit1 )
                {
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
                        //std::cout<<musegmatches.size()<<std::endl;
                        
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
                
                
                
                if(fabs(trackReftmp->pt()-pt2)<0.0001 && fabs(trackReftmp->eta()-eta2)<0.001 && fabs(trackReftmp->phi()-phi2)<0.001
                   && trackReftmp->numberOfValidHits() == nhit2 )
                {
                    nmatchedch2 = cand.numberOfMatches();
                    nmatchedst2 = cand.numberOfMatchedStations();
                    
                    reco::MuonEnergy muenergy = cand.calEnergy();
                    matchedenergy2 = muenergy.hadMax;
                    
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
                        //std::cout<<musegmatches.size()<<std::endl;
                        
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
        }
        
        if(twoLayerDecay_)
        {
            grand_mass = d1->mass();
            
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
            grand_charge1 = gd1->charge();
            grand_charge2 = gd2->charge();
            
            //track Chi2
            grand_trkChi1 = gdau1->normalizedChi2();
            grand_trkChi2 = gdau2->normalizedChi2();
            
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
        theVertexComps.push_back( trk );
    }
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
