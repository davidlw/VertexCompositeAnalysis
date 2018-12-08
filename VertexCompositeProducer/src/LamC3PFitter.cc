// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      LamC3PFitter
// 
/**\class LamC3PFitter LamC3PFitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/LamC3PFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/LamC3PFitter.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>
#include <TVector3.h>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

const float piMassLamC3P = 0.13957018;
const float piMassLamC3PSquared = piMassLamC3P*piMassLamC3P;
const float kaonMassLamC3P = 0.493677;
const float kaonMassLamC3PSquared = kaonMassLamC3P*kaonMassLamC3P;
const float protonMassLamC3P = 0.938272013; 
const float protonMassLamC3PSquared = protonMassLamC3P*protonMassLamC3P;
const float lamCMassLamC3P = 2.28646;
float piMassLamC3P_sigma = 3.5E-7f;
float kaonMassLamC3P_sigma = 1.6E-5f;
float protonMassLamC3P_sigma = 1.6E-5f;
float lamCMassLamC3P_sigma = lamCMassLamC3P*1.e-6;

float cand1Mass[2] = {kaonMassLamC3P, protonMassLamC3P};
float cand2Mass[2] = {protonMassLamC3P, kaonMassLamC3P};
float cand1Mass_sigma[2] = {kaonMassLamC3P_sigma, protonMassLamC3P_sigma};
float cand2Mass_sigma[2] = {protonMassLamC3P_sigma, kaonMassLamC3P_sigma};

// Constructor and (empty) destructor
LamC3PFitter::LamC3PFitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) {
//		   const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::ConsumesCollector && iC) {
  using std::string;

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_dedx = iC.consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));

  // Second, initialize post-fit cuts
  mKPCutMin = theParameters.getParameter<double>(string("mKPCutMin"));
  mKPCutMax = theParameters.getParameter<double>(string("mKPCutMax"));
  mPiKPCutMin = theParameters.getParameter<double>(string("mPiKPCutMin"));
  mPiKPCutMax = theParameters.getParameter<double>(string("mPiKPCutMax"));
  tkDCACut = theParameters.getParameter<double>(string("tkDCACut"));
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  tkPtCut = theParameters.getParameter<double>(string("tkPtCut"));
  tkPtErrCut = theParameters.getParameter<double>(string("tkPtErrCut"));
  tkEtaCut = theParameters.getParameter<double>(string("tkEtaCut"));
//  tkPtSumCut = theParameters.getParameter<double>(string("tkPtSumCut"));
//  tkEtaDiffCut = theParameters.getParameter<double>(string("tkEtaDiffCut"));
  chi2Cut = theParameters.getParameter<double>(string("vtxChi2Cut"));
  rVtxCut = theParameters.getParameter<double>(string("rVtxCut"));
  rVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance2DCut"));
  lVtxCut = theParameters.getParameter<double>(string("lVtxCut"));
  lVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance3DCut"));
  collinCut2D = theParameters.getParameter<double>(string("collinearityCut2D"));
  collinCut3D = theParameters.getParameter<double>(string("collinearityCut3D"));
  lamCMassCut = theParameters.getParameter<double>(string("lamCMassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  VtxChiProbCut = theParameters.getParameter<double>(string("VtxChiProbCut"));
  dPt3Cut = theParameters.getParameter<double>(string("dPt3Cut"));
  alphaCut = theParameters.getParameter<double>(string("alphaCut"));
  alpha2DCut = theParameters.getParameter<double>(string("alpha2DCut"));
  isWrongSign = theParameters.getParameter<bool>(string("isWrongSign"));


  useAnyMVA_ = false;
  forestLabel_ = "LamC3PInpPb";
  std::string type = "BDT";
  useForestFromDB_ = true;
  dbFileName_ = "";

  forest_ = nullptr;

  if(theParameters.exists("useAnyMVA")) useAnyMVA_ = theParameters.getParameter<bool>("useAnyMVA");

  if(useAnyMVA_){
    if(theParameters.exists("mvaType"))type = theParameters.getParameter<std::string>("mvaType");
    if(theParameters.exists("GBRForestLabel"))forestLabel_ = theParameters.getParameter<std::string>("GBRForestLabel");
    if(theParameters.exists("GBRForestFileName")){
      dbFileName_ = theParameters.getParameter<std::string>("GBRForestFileName");
      useForestFromDB_ = false;
    }

    if(!useForestFromDB_){
      edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeProducer/data/%s",dbFileName_.c_str()));
      TFile gbrfile(fip.fullPath().c_str(),"READ");
      forest_ = (GBRForest*)gbrfile.Get(forestLabel_.c_str());
      gbrfile.Close();
    }

    mvaType_ = type;
  }

  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }
}

LamC3PFitter::~LamC3PFitter() {
  delete forest_;
}

// Method containing the algorithm for vertex reconstruction
void LamC3PFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;
  using namespace std; 

  // Create std::vectors for Tracks and TrackRefs (required for
  //  passing to the KalmanVertexFitter)
  std::vector<TrackRef> theTrackRefs_pos;
  std::vector<TrackRef> theTrackRefs_neg;
  std::vector<TransientTrack> theTransTracks_pos;
  std::vector<TransientTrack> theTransTracks_neg;

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;
  Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks, theTrackHandle); 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);  
  iEvent.getByToken(token_dedx, dEdxHandle);

  if( !theTrackHandle->size() ) return;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  magField = bFieldHandle.product();

  // Setup TMVA
//  mvaValValueMap = auto_ptr<edm::ValueMap<float> >(new edm::ValueMap<float>);
//  edm::ValueMap<float>::Filler mvaFiller(*mvaValValueMap);

  bool isVtxPV = 0;
  double xVtx=-99999.0;
  double yVtx=-99999.0;
  double zVtx=-99999.0;
  double xVtxError=-999.0;
  double yVtxError=-999.0;
  double zVtxError=-999.0;
  const reco::VertexCollection vtxCollection = *(theVertexHandle.product());
  reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();
  if(vtxCollection.size()>0 && !vtxPrimary->isFake() && vtxPrimary->tracksSize()>=2)
  {
    isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  }
  else {
    isVtxPV = 0;
    xVtx = theBeamSpotHandle->position().x();
    yVtx = theBeamSpotHandle->position().y();
    zVtx = 0.0;
    xVtxError = theBeamSpotHandle->BeamWidthX();
    yVtxError = theBeamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }
  math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
  math::XYZPoint bestvtxError(xVtxError,yVtxError,zVtxError);

  // Fill vectors of TransientTracks and TrackRefs after applying preselection cuts.
  for(unsigned int indx = 0; indx < theTrackHandle->size(); indx++) {
    TrackRef tmpRef( theTrackHandle, indx );
    bool quality_ok = true;
    if (qualities.size()!=0) {
      quality_ok = false;
      for (unsigned int ndx_ = 0; ndx_ < qualities.size(); ndx_++) {
	if (tmpRef->quality(qualities[ndx_])){
	  quality_ok = true;
	  break;          
	}
      }
    }
    if( !quality_ok ) continue;

    if( tmpRef->normalizedChi2() < tkChi2Cut &&
        tmpRef->numberOfValidHits() >= tkNhitsCut &&
        tmpRef->ptError() / tmpRef->pt() < tkPtErrCut &&
        tmpRef->pt() > tkPtCut && fabs(tmpRef->eta()) < tkEtaCut ) {
//      TransientTrack tmpTk( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
      TransientTrack tmpTk( *tmpRef, magField );

      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauLongImpactSig = dzvtx/dzerror;
      double dauTransImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut ) {
        if(tmpRef->charge()>0.0)
        {
          theTrackRefs_pos.push_back( tmpRef );
          theTransTracks_pos.push_back( tmpTk );
        }
        if(tmpRef->charge()<0.0)
        {
          theTrackRefs_neg.push_back( tmpRef );
          theTransTracks_neg.push_back( tmpTk );
        }
      }
    }
  }

  if(!isWrongSign)
  {
    fitLamCCandidates(theTrackRefs_pos,theTrackRefs_neg,theTransTracks_pos,theTransTracks_neg,isVtxPV,vtxPrimary,theBeamSpotHandle,bestvtx,bestvtxError,4122);
    fitLamCCandidates(theTrackRefs_neg,theTrackRefs_pos,theTransTracks_neg,theTransTracks_pos,isVtxPV,vtxPrimary,theBeamSpotHandle,bestvtx,bestvtxError,-4122);
  }
  else 
  {
    fitLamCCandidates(theTrackRefs_pos,theTrackRefs_pos,theTransTracks_pos,theTransTracks_pos,isVtxPV,vtxPrimary,theBeamSpotHandle,bestvtx,bestvtxError,4122);
    fitLamCCandidates(theTrackRefs_neg,theTrackRefs_neg,theTransTracks_neg,theTransTracks_neg,isVtxPV,vtxPrimary,theBeamSpotHandle,bestvtx,bestvtxError,-4122);    
  }
}

void LamC3PFitter::fitLamCCandidates(
                                  std::vector<reco::TrackRef> theTrackRefs_sgn1,
                                  std::vector<reco::TrackRef> theTrackRefs_sgn2,
                                  std::vector<reco::TransientTrack> theTransTracks_sgn1,
                                  std::vector<reco::TransientTrack> theTransTracks_sgn2,
                                  bool isVtxPV,
                                  reco::VertexCollection::const_iterator vtxPrimary, edm::Handle<reco::BeamSpot> theBeamSpotHandle,
                                  math::XYZPoint bestvtx, math::XYZPoint bestvtxError,
                                  int pdg_id
                                )
{
  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;
  using namespace std;

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  int lamCCharge = pdg_id/abs(pdg_id);

  // Loop over tracks and vertex good charged track pairs
  for(unsigned int trdx1 = 0; trdx1 < theTrackRefs_sgn1.size(); trdx1++) {

    for(unsigned int trdx2 = trdx1 + 1; trdx2 < theTrackRefs_sgn1.size(); trdx2++) {
//      if( (theTrackRefs[trdx1]->pt() + theTrackRefs[trdx2]->pt()) < tkPtSumCut) continue;
//      if( abs(theTrackRefs[trdx1]->eta() - theTrackRefs[trdx2]->eta()) > tkEtaDiffCut) continue;

      //This vector holds the pair of oppositely-charged tracks to be vertexed
      std::vector<TransientTrack> transTracks;

      TrackRef trackRef1 = theTrackRefs_sgn1[trdx1];
      TrackRef trackRef2 = theTrackRefs_sgn1[trdx2];
      TransientTrack* transTkPtr1 = &theTransTracks_sgn1[trdx1];
      TransientTrack* transTkPtr2 = &theTransTracks_sgn1[trdx2];

      double dzvtx1 = trackRef1->dz(bestvtx);
      double dxyvtx1 = trackRef1->dxy(bestvtx);
      double dzerror1 = sqrt(trackRef1->dzError()*trackRef1->dzError()+bestvtxError.z()*bestvtxError.z());
      double dxyerror1 = sqrt(trackRef1->d0Error()*trackRef1->d0Error()+bestvtxError.x()*bestvtxError.y());
      double dauLongImpactSig1 = dzvtx1/dzerror1;
      double dauTransImpactSig1 = dxyvtx1/dxyerror1;

      double dzvtx2 = trackRef2->dz(bestvtx);
      double dxyvtx2 = trackRef2->dxy(bestvtx);
      double dzerror2 = sqrt(trackRef2->dzError()*trackRef2->dzError()+bestvtxError.z()*bestvtxError.z());
      double dxyerror2 = sqrt(trackRef2->d0Error()*trackRef2->d0Error()+bestvtxError.x()*bestvtxError.y());
      double dauLongImpactSig2 = dzvtx2/dzerror2;
      double dauTransImpactSig2 = dxyvtx2/dxyerror2;

      double nhits1 = trackRef1->numberOfValidHits();
      double nhits2 = trackRef2->numberOfValidHits(); 
    
      double ptErr1 = trackRef1->ptError();
      double ptErr2 = trackRef2->ptError();

      // Fill the vector of TransientTracks to send to KVF
      transTracks.push_back(*transTkPtr1);
      transTracks.push_back(*transTkPtr2);

      // Trajectory states to calculate DCA for the 2 tracks
      FreeTrajectoryState trkState1 = transTkPtr1->impactPointTSCP().theState();
      FreeTrajectoryState trkState2 = transTkPtr2->impactPointTSCP().theState();

      if( !transTkPtr1->impactPointTSCP().isValid() || !transTkPtr2->impactPointTSCP().isValid() ) continue;

      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(trkState1, trkState2);
      if( !cApp.status() ) continue;
      float dca = fabs( cApp.distance() );
      GlobalPoint cxPt = cApp.crossingPoint();

      if (dca < 0. || dca > tkDCACut) continue;

      // Get trajectory states for the tracks at POCA for later cuts
      TrajectoryStateClosestToPoint trkTSCP1 =
        transTkPtr1->trajectoryStateClosestToPoint( cxPt );
      TrajectoryStateClosestToPoint trkTSCP2 =
        transTkPtr2->trajectoryStateClosestToPoint( cxPt );

      if( !trkTSCP1.isValid() || !trkTSCP2.isValid() ) continue;

      double totalE1 = sqrt( trkTSCP1.momentum().mag2() + protonMassLamC3PSquared ) +
                      sqrt( trkTSCP2.momentum().mag2() + kaonMassLamC3PSquared );
      double totalE1Sq = totalE1*totalE1;

      double totalE2 = sqrt( trkTSCP1.momentum().mag2() + kaonMassLamC3PSquared ) +
                      sqrt( trkTSCP2.momentum().mag2() + protonMassLamC3PSquared );
      double totalE2Sq = totalE2*totalE2;

      double totalPSq =
        ( trkTSCP1.momentum() + trkTSCP2.momentum() ).mag2();

      double totalPt =
        ( trkTSCP1.momentum() + trkTSCP2.momentum() ).perp();

      double mass1 = sqrt( totalE1Sq - totalPSq);
      double mass2 = sqrt( totalE2Sq - totalPSq);

      if( (mass1 > mKPCutMax || mass1 < mKPCutMin) && (mass2 > mKPCutMax || mass2 < mKPCutMin)) continue;
//      if( totalPt < dPtCut ) continue;

      for(unsigned int trdx3 = 0; trdx3 < theTrackRefs_sgn2.size(); trdx3++) {

        TrackRef trackRef3 = theTrackRefs_sgn2[trdx3];
        TransientTrack* transTkPtr3 = &theTransTracks_sgn2[trdx3];
  
        double dzvtx3 = trackRef3->dz(bestvtx);
        double dxyvtx3 = trackRef3->dxy(bestvtx);
        double dzerror3 = sqrt(trackRef3->dzError()*trackRef3->dzError()+bestvtxError.z()*bestvtxError.z());
        double dxyerror3 = sqrt(trackRef3->d0Error()*trackRef3->d0Error()+bestvtxError.x()*bestvtxError.y());
        double dauLongImpactSig3 = dzvtx3/dzerror3;
        double dauTransImpactSig3 = dxyvtx3/dxyerror3;

        double nhits3 = trackRef3->numberOfValidHits();
        double ptErr3 = trackRef3->ptError();

        transTracks.push_back(*transTkPtr3);
        FreeTrajectoryState trkState3 = transTkPtr3->impactPointTSCP().theState();
        if( !transTkPtr3->impactPointTSCP().isValid() ) continue;

        // Measure distance between tracks at their closest approach
        ClosestApproachInRPhi cApp13;
        cApp13.calculate(trkState1, trkState3);
        if( !cApp13.status() ) continue;
        float dca13 = fabs( cApp13.distance() );
        GlobalPoint cxPt13 = cApp13.crossingPoint();
        if (dca13 < 0. || dca13 > tkDCACut) continue;

        // Get trajectory states for the tracks at POCA for later cuts
        TrajectoryStateClosestToPoint trkTSCP31 =
          transTkPtr3->trajectoryStateClosestToPoint( cxPt13 );

        if( !trkTSCP31.isValid() ) continue;

        double totalE31 = sqrt( trkTSCP1.momentum().mag2() + protonMassLamC3PSquared ) +
                          sqrt( trkTSCP2.momentum().mag2() + kaonMassLamC3PSquared ) + 
                          sqrt( trkTSCP31.momentum().mag2() + piMassLamC3PSquared );
        double totalE31Sq = totalE31*totalE31;

        double totalE32 = sqrt( trkTSCP1.momentum().mag2() + kaonMassLamC3PSquared ) +
                          sqrt( trkTSCP2.momentum().mag2() + protonMassLamC3PSquared ) + 
                          sqrt( trkTSCP31.momentum().mag2() + piMassLamC3PSquared );
        double totalE32Sq = totalE32*totalE32;

        double totalP3Sq =
          ( trkTSCP1.momentum() + trkTSCP2.momentum() + trkTSCP31.momentum()).mag2();

        double totalPt3 =
          ( trkTSCP1.momentum() + trkTSCP2.momentum() + trkTSCP31.momentum()).perp();

        double mass31 = sqrt( totalE31Sq - totalP3Sq);
        double mass32 = sqrt( totalE32Sq - totalP3Sq);

        if( (mass31 > mPiKPCutMax || mass31 < mPiKPCutMin) && (mass32 > mPiKPCutMax || mass32 < mPiKPCutMin)) continue;
        if( totalPt3 < dPt3Cut ) continue;

        // Create the vertex fitter object and vertex the tracks
        float cand1TotalE[2]={0.0};
        float cand2TotalE[2]={0.0};
        float lamCTotalE[2]={0.0};

        for(int i=0;i<2;i++)
        {
          //Creating a KinematicParticleFactory
          KinematicParticleFactoryFromTransientTrack pFactory;
        
          float chi = 0.0;
          float ndf = 0.0;

          vector<RefCountedKinematicParticle> lamCParticles;
          lamCParticles.push_back(pFactory.particle(*transTkPtr1,cand1Mass[i],chi,ndf,cand1Mass_sigma[i]));
          lamCParticles.push_back(pFactory.particle(*transTkPtr2,cand2Mass[i],chi,ndf,cand2Mass_sigma[i]));
          lamCParticles.push_back(pFactory.particle(*transTkPtr3,piMassLamC3P,chi,ndf,piMassLamC3P_sigma));

          KinematicParticleVertexFitter lamCFitter;
          RefCountedKinematicTree lamCVertex;
          lamCVertex = lamCFitter.fit(lamCParticles);

          if( !lamCVertex->isValid() ) continue;

          lamCVertex->movePointerToTheTop();
          RefCountedKinematicParticle lamCCand = lamCVertex->currentParticle();
          if(!lamCCand->currentState().isValid()) continue;

          RefCountedKinematicVertex lamCDecayVertex = lamCVertex->currentDecayVertex();
          if(!lamCDecayVertex->vertexIsValid()) continue;

  	  float lamCC2Prob = TMath::Prob(lamCDecayVertex->chiSquared(),lamCDecayVertex->degreesOfFreedom());
  	  if (lamCC2Prob < VtxChiProbCut) continue;

          lamCVertex->movePointerToTheFirstChild();
          RefCountedKinematicParticle cand1 = lamCVertex->currentParticle();
          lamCVertex->movePointerToTheNextChild();
          RefCountedKinematicParticle cand2 = lamCVertex->currentParticle();
          lamCVertex->movePointerToTheNextChild();
          RefCountedKinematicParticle cand3 = lamCVertex->currentParticle();

          if(!cand1->currentState().isValid() || !cand2->currentState().isValid() || !cand3->currentState().isValid()) continue;

          KinematicParameters cand1KP = cand1->currentState().kinematicParameters();
          KinematicParameters cand2KP = cand2->currentState().kinematicParameters();
          KinematicParameters cand3KP = cand3->currentState().kinematicParameters();

          GlobalVector lamCTotalP = GlobalVector (lamCCand->currentState().globalMomentum().x(),
                                                  lamCCand->currentState().globalMomentum().y(),
                                                  lamCCand->currentState().globalMomentum().z());

          GlobalVector cand1TotalP = GlobalVector(cand1KP.momentum().x(),cand1KP.momentum().y(),cand1KP.momentum().z());
          GlobalVector cand2TotalP = GlobalVector(cand2KP.momentum().x(),cand2KP.momentum().y(),cand2KP.momentum().z());
          GlobalVector cand3TotalP = GlobalVector(cand3KP.momentum().x(),cand3KP.momentum().y(),cand3KP.momentum().z());

          cand1TotalE[i] = sqrt( cand1TotalP.mag2() + cand1Mass[i]*cand1Mass[i] );
          cand2TotalE[i] = sqrt( cand2TotalP.mag2() + cand2Mass[i]*cand2Mass[i] );
          float cand3TotalE = sqrt( cand3TotalP.mag2() + piMassLamC3PSquared );

          lamCTotalE[i] = cand1TotalE[i] + cand2TotalE[i] + cand3TotalE;

          const Particle::LorentzVector lamCP4(lamCTotalP.x(), lamCTotalP.y(), lamCTotalP.z(), lamCTotalE[i]);

          Particle::Point lamCVtx((*lamCDecayVertex).position().x(), (*lamCDecayVertex).position().y(), (*lamCDecayVertex).position().z());
          std::vector<double> lamCVtxEVec;
          lamCVtxEVec.push_back( lamCDecayVertex->error().cxx() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().cyx() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().cyy() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().czx() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().czy() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().czz() );
          SMatrixSym3D lamCVtxCovMatrix(lamCVtxEVec.begin(), lamCVtxEVec.end());
          const Vertex::CovarianceMatrix lamCVtxCov(lamCVtxCovMatrix);
          double lamCVtxChi2(lamCDecayVertex->chiSquared());
          double lamCVtxNdof(lamCDecayVertex->degreesOfFreedom());
          double lamCNormalizedChi2 = lamCVtxChi2/lamCVtxNdof;

          double rVtxMag = 99999.0; 
          double lVtxMag = 99999.0;
          double sigmaRvtxMag = 999.0;
          double sigmaLvtxMag = 999.0;
          double lamCAngle3D = -100.0;
          double lamCAngle2D = -100.0;

          GlobalVector lamCLineOfFlight = GlobalVector (lamCVtx.x() - bestvtx.x(),
                                                        lamCVtx.y() - bestvtx.y(),
                                                        lamCVtx.z() - bestvtx.z());

          SMatrixSym3D lamCTotalCov;
          if(isVtxPV) lamCTotalCov = lamCVtxCovMatrix + vtxPrimary->covariance();
          else lamCTotalCov = lamCVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

          SVector3 distanceVector3D(lamCLineOfFlight.x(), lamCLineOfFlight.y(), lamCLineOfFlight.z());
          SVector3 distanceVector2D(lamCLineOfFlight.x(), lamCLineOfFlight.y(), 0.0);

          lamCAngle3D = angle(lamCLineOfFlight.x(), lamCLineOfFlight.y(), lamCLineOfFlight.z(),
                          lamCTotalP.x(), lamCTotalP.y(), lamCTotalP.z());
          lamCAngle2D = angle(lamCLineOfFlight.x(), lamCLineOfFlight.y(), (float)0.0,
                          lamCTotalP.x(), lamCTotalP.y(), (float)0.0);

          lVtxMag = lamCLineOfFlight.mag();
          rVtxMag = lamCLineOfFlight.perp();
          sigmaLvtxMag = sqrt(ROOT::Math::Similarity(lamCTotalCov, distanceVector3D)) / lVtxMag;
          sigmaRvtxMag = sqrt(ROOT::Math::Similarity(lamCTotalCov, distanceVector2D)) / rVtxMag;

          if( lamCNormalizedChi2 > chi2Cut ||
              rVtxMag < rVtxCut ||
              rVtxMag / sigmaRvtxMag < rVtxSigCut ||
              lVtxMag < lVtxCut ||
              lVtxMag / sigmaLvtxMag < lVtxSigCut ||
              cos(lamCAngle3D) < collinCut3D || cos(lamCAngle2D) < collinCut2D || lamCAngle3D > alphaCut || lamCAngle2D > alpha2DCut
          ) continue;

          VertexCompositeCandidate* theLamC3P = 0;
          theLamC3P = new VertexCompositeCandidate(lamCCharge, lamCP4, lamCVtx, lamCVtxCov, lamCVtxChi2, lamCVtxNdof);

          RecoChargedCandidate
            theCand1(lamCCharge, Particle::LorentzVector(cand1TotalP.x(),
                                                  cand1TotalP.y(), cand1TotalP.z(),
                                                  cand1TotalE[i]), lamCVtx);
            theCand1.setTrack(trackRef1);

          RecoChargedCandidate
            theCand2(lamCCharge, Particle::LorentzVector(cand2TotalP.x(),
                                                   cand2TotalP.y(), cand2TotalP.z(),
                                                   cand2TotalE[i]), lamCVtx);
            theCand2.setTrack(trackRef2);

          RecoChargedCandidate
            theCand3(-lamCCharge, Particle::LorentzVector(cand3TotalP.x(),
                                                   cand3TotalP.y(), cand3TotalP.z(),
                                                   cand3TotalE), lamCVtx);
            theCand3.setTrack(trackRef3);

          AddFourMomenta addp4;
          theLamC3P->addDaughter(theCand1);
          theLamC3P->addDaughter(theCand2);
          theLamC3P->addDaughter(theCand3);
          theLamC3P->setPdgId(pdg_id);
          addp4.set( *theLamC3P );
          if( theLamC3P->mass() < lamCMassLamC3P + lamCMassCut &&
              theLamC3P->mass() > lamCMassLamC3P - lamCMassCut ) //&&
	     // theLamC3P->pt() > dPtCut ) {
          {
            theLamC3Ps.push_back( *theLamC3P );
          }
// perform MVA evaluation
/*
          if(useAnyMVA_)
          {
            float gbrVals_[20];
            gbrVals_[0] = lamCP4.Pt();
            gbrVals_[1] = lamCP4.Eta();
            gbrVals_[2] = lamCC2Prob;
            gbrVals_[3] = lVtxMag / sigmaLvtxMag;
            gbrVals_[4] = rVtxMag / sigmaRvtxMag;
            gbrVals_[5] = lVtxMag;
            gbrVals_[6] = lamCAngle3D;
            gbrVals_[7] = lamCAngle2D;
            gbrVals_[8] = dauLongImpactSig_pos;
            gbrVals_[9] = dauLongImpactSig_neg;
            gbrVals_[10] = dauTransImpactSig_pos;
            gbrVals_[11] = dauTransImpactSig_neg;
            gbrVals_[12] = nhits_pos;
            gbrVals_[13] = nhits_neg;
            gbrVals_[14] = ptErr_pos;
            gbrVals_[15] = ptErr_neg;
            gbrVals_[16] = posCandTotalP.perp();
            gbrVals_[17] = negCandTotalP.perp();
            gbrVals_[18] = posCandTotalP.eta();
            gbrVals_[19] = negCandTotalP.eta();

            GBRForest const * forest = forest_;
            if(useForestFromDB_){
              edm::ESHandle<GBRForest> forestHandle;
              iSetup.get<GBRWrapperRcd>().get(forestLabel_,forestHandle);
              forest = forestHandle.product();
            }

            auto gbrVal = forest->GetClassifier(gbrVals_);
            mvaVals_.push_back(gbrVal);
          }
*/
          if(theLamC3P) delete theLamC3P;
        } // swap mass
      } // trk3 
    }  // trk2
  } // trk1

//  mvaFiller.insert(theLamC3Ps,mvaVals_.begin(),mvaVals_.end());
//  mvaFiller.fill();
//  mvas = std::make_unique<MVACollection>(mvaVals_.begin(),mvaVals_.end());

}
// Get methods

const reco::VertexCompositeCandidateCollection& LamC3PFitter::getLamC3P() const {
  return theLamC3Ps;
}

const std::vector<float>& LamC3PFitter::getMVAVals() const {
  return mvaVals_;
}

/*
auto_ptr<edm::ValueMap<float> > LamC3PFitter::getMVAMap() const {
  return mvaValValueMap;
}
*/

void LamC3PFitter::resetAll() {
    theLamC3Ps.clear();
    mvaVals_.clear();
}
