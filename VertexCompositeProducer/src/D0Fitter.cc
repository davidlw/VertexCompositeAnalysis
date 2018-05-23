// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      D0Fitter
// 
/**\class D0Fitter D0Fitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/D0Fitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/D0Fitter.h"
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

#include <typeinfo>
#include <memory>
#include <vector>

const float piMassD0 = 0.13957018;
const float piMassD0Squared = piMassD0*piMassD0;
const float kaonMassD0 = 0.493677;
const float kaonMassD0Squared = kaonMassD0*kaonMassD0;
const float d0MassD0 = 1.86484;
float piMassD0_sigma = 3.5E-7f;
float kaonMassD0_sigma = 1.6E-5f;
float d0MassD0_sigma = d0MassD0*1.e-6;

// Constructor and (empty) destructor
D0Fitter::D0Fitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) {
//		   const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::ConsumesCollector && iC) {
  using std::string;

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));

  // Second, initialize post-fit cuts
  mPiKCutMin = theParameters.getParameter<double>(string("mPiKCutMin"));
  mPiKCutMax = theParameters.getParameter<double>(string("mPiKCutMax"));
  tkDCACut = theParameters.getParameter<double>(string("tkDCACut"));
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  tkPtCut = theParameters.getParameter<double>(string("tkPtCut"));
  tkEtaCut = theParameters.getParameter<double>(string("tkEtaCut"));
  chi2Cut = theParameters.getParameter<double>(string("vtxChi2Cut"));
  rVtxCut = theParameters.getParameter<double>(string("rVtxCut"));
  rVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance2DCut"));
  lVtxCut = theParameters.getParameter<double>(string("lVtxCut"));
  lVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance3DCut"));
  collinCut = theParameters.getParameter<double>(string("collinearityCut"));
  d0MassCut = theParameters.getParameter<double>(string("d0MassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  VtxChiProbCut = theParameters.getParameter<double>(string("VtxChiProbCut"));
  dPtCut = theParameters.getParameter<double>(string("dPtCut"));
  alphaCut = theParameters.getParameter<double>(string("alphaCut"));
  isWrongSign = theParameters.getParameter<bool>(string("isWrongSign"));

  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }
}

D0Fitter::~D0Fitter() {
}

// Method containing the algorithm for vertex reconstruction
void D0Fitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  // Create std::vectors for Tracks and TrackRefs (required for
  //  passing to the KalmanVertexFitter)
  std::vector<TrackRef> theTrackRefs;
  std::vector<TransientTrack> theTransTracks;

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks, theTrackHandle); 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);  

  if( !theTrackHandle->size() ) return;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  magField = bFieldHandle.product();

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
        tmpRef->pt() > tkPtCut && fabs(tmpRef->eta()) < tkEtaCut ) {
//      TransientTrack tmpTk( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
      TransientTrack tmpTk( *tmpRef, magField );

      math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauLongImpactSig = dzvtx/dzerror;
      double dauTransImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut ) {
        theTrackRefs.push_back( tmpRef );
        theTransTracks.push_back( tmpTk );
      }
    }
  }

  float posCandMass[2] = {piMassD0, kaonMassD0};
  float negCandMass[2] = {kaonMassD0, piMassD0};
  float posCandMass_sigma[2] = {piMassD0_sigma, kaonMassD0_sigma};
  float negCandMass_sigma[2] = {kaonMassD0_sigma, piMassD0_sigma};
  int   pdg_id[2] = {421, -421};

  // Loop over tracks and vertex good charged track pairs
  for(unsigned int trdx1 = 0; trdx1 < theTrackRefs.size(); trdx1++) {

    for(unsigned int trdx2 = trdx1 + 1; trdx2 < theTrackRefs.size(); trdx2++) {

      //This vector holds the pair of oppositely-charged tracks to be vertexed
      std::vector<TransientTrack> transTracks;

      TrackRef positiveTrackRef;
      TrackRef negativeTrackRef;
      TransientTrack* posTransTkPtr = 0;
      TransientTrack* negTransTkPtr = 0;

      // Look at the two tracks we're looping over.  If they're oppositely
      //  charged, load them into the hypothesized positive and negative tracks
      //  and references to be sent to the KalmanVertexFitter
      if(!isWrongSign && theTrackRefs[trdx1]->charge() < 0. && 
	 theTrackRefs[trdx2]->charge() > 0.) {
	negativeTrackRef = theTrackRefs[trdx1];
	positiveTrackRef = theTrackRefs[trdx2];
	negTransTkPtr = &theTransTracks[trdx1];
	posTransTkPtr = &theTransTracks[trdx2];
      }
      else if(!isWrongSign && theTrackRefs[trdx1]->charge() > 0. &&
	      theTrackRefs[trdx2]->charge() < 0.) {
	negativeTrackRef = theTrackRefs[trdx2];
	positiveTrackRef = theTrackRefs[trdx1];
	negTransTkPtr = &theTransTracks[trdx2];
	posTransTkPtr = &theTransTracks[trdx1];
      }
      else if(isWrongSign && theTrackRefs[trdx1]->charge()*theTrackRefs[trdx2]->charge() > 0.0) { 
        negativeTrackRef = theTrackRefs[trdx2];
        positiveTrackRef = theTrackRefs[trdx1];
        negTransTkPtr = &theTransTracks[trdx2];
        posTransTkPtr = &theTransTracks[trdx1];
      }
      // If they're not 2 oppositely charged tracks, loop back to the
      //  beginning and try the next pair.
      else continue;

      // Fill the vector of TransientTracks to send to KVF
      transTracks.push_back(*posTransTkPtr);
      transTracks.push_back(*negTransTkPtr);

      // Trajectory states to calculate DCA for the 2 tracks
      FreeTrajectoryState posState = posTransTkPtr->impactPointTSCP().theState();
      FreeTrajectoryState negState = negTransTkPtr->impactPointTSCP().theState();

      if( !posTransTkPtr->impactPointTSCP().isValid() || !negTransTkPtr->impactPointTSCP().isValid() ) continue;

      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(posState, negState);
      if( !cApp.status() ) continue;
      float dca = fabs( cApp.distance() );
      GlobalPoint cxPt = cApp.crossingPoint();

      if (dca < 0. || dca > tkDCACut) continue;
//      if (sqrt( cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y() ) > 120. 
//          || std::abs(cxPt.z()) > 300.) continue;

      // Get trajectory states for the tracks at POCA for later cuts
      TrajectoryStateClosestToPoint posTSCP =
        posTransTkPtr->trajectoryStateClosestToPoint( cxPt );
      TrajectoryStateClosestToPoint negTSCP =
        negTransTkPtr->trajectoryStateClosestToPoint( cxPt );

      if( !posTSCP.isValid() || !negTSCP.isValid() ) continue;

      double totalE1 = sqrt( posTSCP.momentum().mag2() + kaonMassD0Squared ) +
                      sqrt( negTSCP.momentum().mag2() + piMassD0Squared );
      double totalE1Sq = totalE1*totalE1;

      double totalE2 = sqrt( posTSCP.momentum().mag2() + piMassD0Squared ) +
                      sqrt( negTSCP.momentum().mag2() + kaonMassD0Squared );
      double totalE2Sq = totalE2*totalE2;

      double totalPSq =
        ( posTSCP.momentum() + negTSCP.momentum() ).mag2();

      double mass1 = sqrt( totalE1Sq - totalPSq);
      double mass2 = sqrt( totalE2Sq - totalPSq);

      if( (mass1 > mPiKCutMax || mass1 < mPiKCutMin) && (mass2 > mPiKCutMax || mass2 < mPiKCutMin)) continue;

      // Create the vertex fitter object and vertex the tracks
    
      float posCandTotalE[2]={0.0};
      float negCandTotalE[2]={0.0};
      float d0TotalE[2]={0.0};

      for(int i=0;i<2;i++)
      {
        //Creating a KinematicParticleFactory
        KinematicParticleFactoryFromTransientTrack pFactory;
        
        float chi = 0.0;
        float ndf = 0.0;

        vector<RefCountedKinematicParticle> d0Particles;
        d0Particles.push_back(pFactory.particle(*posTransTkPtr,posCandMass[i],chi,ndf,posCandMass_sigma[i]));
        d0Particles.push_back(pFactory.particle(*negTransTkPtr,negCandMass[i],chi,ndf,negCandMass_sigma[i]));

        KinematicParticleVertexFitter d0Fitter;
        RefCountedKinematicTree d0Vertex;
        d0Vertex = d0Fitter.fit(d0Particles);

        if( !d0Vertex->isValid() ) continue;

        d0Vertex->movePointerToTheTop();
        RefCountedKinematicParticle d0Cand = d0Vertex->currentParticle();
        if (!d0Cand->currentState().isValid()) continue;

        RefCountedKinematicVertex d0DecayVertex = d0Vertex->currentDecayVertex();
        if (!d0DecayVertex->vertexIsValid()) continue;

        //if ( d0DecayVertex->chiSquared()<0 || d0DecayVertex->chiSquared()>1000 ) continue;

        //float d0C2Prob =
        //   ChiSquaredProbability((double)(d0DecayVertex->chiSquared()),(double)(d0DecayVertex->degreesOfFreedom()));
        //if (d0C2Prob < 0.0001) continue;

	float d0C2Prob = TMath::Prob(d0DecayVertex->chiSquared(),d0DecayVertex->degreesOfFreedom());
	if (d0C2Prob < VtxChiProbCut) continue;

        //if ( d0Cand->currentState().mass() > 2.5 || d0Cand->currentState().mass() < 1.0) continue;

        d0Vertex->movePointerToTheFirstChild();
        RefCountedKinematicParticle posCand = d0Vertex->currentParticle();
        d0Vertex->movePointerToTheNextChild();
        RefCountedKinematicParticle negCand = d0Vertex->currentParticle();

        if(!posCand->currentState().isValid() || !negCand->currentState().isValid()) continue;

        KinematicParameters posCandKP = posCand->currentState().kinematicParameters();
        KinematicParameters negCandKP = negCand->currentState().kinematicParameters();

        GlobalVector d0TotalP = GlobalVector (d0Cand->currentState().globalMomentum().x(),
                                              d0Cand->currentState().globalMomentum().y(),
                                              d0Cand->currentState().globalMomentum().z());

        GlobalVector posCandTotalP = GlobalVector(posCandKP.momentum().x(),posCandKP.momentum().y(),posCandKP.momentum().z());
        GlobalVector negCandTotalP = GlobalVector(negCandKP.momentum().x(),negCandKP.momentum().y(),negCandKP.momentum().z());

        posCandTotalE[i] = sqrt( posCandTotalP.mag2() + posCandMass[i]*posCandMass[i] );
        negCandTotalE[i] = sqrt( negCandTotalP.mag2() + negCandMass[i]*negCandMass[i] );
        d0TotalE[i] = posCandTotalE[i] + negCandTotalE[i];

        const Particle::LorentzVector d0P4(d0TotalP.x(), d0TotalP.y(), d0TotalP.z(), d0TotalE[i]);

        Particle::Point d0Vtx((*d0DecayVertex).position().x(), (*d0DecayVertex).position().y(), (*d0DecayVertex).position().z());
        std::vector<double> d0VtxEVec;
        d0VtxEVec.push_back( d0DecayVertex->error().cxx() );
        d0VtxEVec.push_back( d0DecayVertex->error().cyx() );
        d0VtxEVec.push_back( d0DecayVertex->error().cyy() );
        d0VtxEVec.push_back( d0DecayVertex->error().czx() );
        d0VtxEVec.push_back( d0DecayVertex->error().czy() );
        d0VtxEVec.push_back( d0DecayVertex->error().czz() );
        SMatrixSym3D d0VtxCovMatrix(d0VtxEVec.begin(), d0VtxEVec.end());
        const Vertex::CovarianceMatrix d0VtxCov(d0VtxCovMatrix);
        double d0VtxChi2(d0DecayVertex->chiSquared());
        double d0VtxNdof(d0DecayVertex->degreesOfFreedom());
        double d0NormalizedChi2 = d0VtxChi2/d0VtxNdof;

        double rVtxMag = 99999.0; 
        double lVtxMag = 99999.0;
        double sigmaRvtxMag = 999.0;
        double sigmaLvtxMag = 999.0;
        double d0Angle = -100.0;

        GlobalVector d0LineOfFlight = GlobalVector (d0Vtx.x() - xVtx,
                                                    d0Vtx.y() - yVtx,
                                                    d0Vtx.z() - zVtx);

        SMatrixSym3D d0TotalCov;
        if(isVtxPV) d0TotalCov = d0VtxCovMatrix + vtxPrimary->covariance();
        else d0TotalCov = d0VtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

        SVector3 distanceVector3D(d0LineOfFlight.x(), d0LineOfFlight.y(), d0LineOfFlight.z());
        SVector3 distanceVector2D(d0LineOfFlight.x(), d0LineOfFlight.y(), 0.0);

        d0Angle = angle(d0LineOfFlight.x(), d0LineOfFlight.y(), d0LineOfFlight.z(),
                        d0TotalP.x(), d0TotalP.y(), d0TotalP.z());
        lVtxMag = d0LineOfFlight.mag();
        rVtxMag = d0LineOfFlight.perp();
        sigmaLvtxMag = sqrt(ROOT::Math::Similarity(d0TotalCov, distanceVector3D)) / lVtxMag;
        sigmaRvtxMag = sqrt(ROOT::Math::Similarity(d0TotalCov, distanceVector2D)) / rVtxMag;

          TVector3 svpvVec;
          svpvVec.SetXYZ(d0Vtx.x() - xVtx,
                         d0Vtx.y() - yVtx,
                         d0Vtx.z() - zVtx);
          TVector3 dVec;
          dVec.SetXYZ(d0TotalP.x(), d0TotalP.y(), d0TotalP.z());
          double alpha = svpvVec.Angle(dVec);
          
        if( d0NormalizedChi2 > chi2Cut ||
            rVtxMag < rVtxCut ||
            rVtxMag / sigmaRvtxMag < rVtxSigCut ||
            lVtxMag < lVtxCut ||
            lVtxMag / sigmaLvtxMag < lVtxSigCut ||
            cos(d0Angle) < collinCut || alpha > alphaCut
        ) continue;

        VertexCompositeCandidate* theD0 = 0;
        theD0 = new VertexCompositeCandidate(0, d0P4, d0Vtx, d0VtxCov, d0VtxChi2, d0VtxNdof);

        RecoChargedCandidate
          thePosCand(1, Particle::LorentzVector(posCandTotalP.x(),
                                                   posCandTotalP.y(), posCandTotalP.z(),
                                                   posCandTotalE[i]), d0Vtx);
        thePosCand.setTrack(positiveTrackRef);

        RecoChargedCandidate
          theNegCand(-1, Particle::LorentzVector(negCandTotalP.x(),
                                                   negCandTotalP.y(), negCandTotalP.z(),
                                                   negCandTotalE[i]), d0Vtx);
        theNegCand.setTrack(negativeTrackRef);

        if(isWrongSign)
        {
          thePosCand.setCharge(theTrackRefs[trdx1]->charge());
          theNegCand.setCharge(theTrackRefs[trdx1]->charge());
        }

        AddFourMomenta addp4;
        theD0->addDaughter(thePosCand);
        theD0->addDaughter(theNegCand);
        theD0->setPdgId(pdg_id[i]);
        addp4.set( *theD0 );
        if( theD0->mass() < d0MassD0 + d0MassCut &&
            theD0->mass() > d0MassD0 - d0MassCut &&
	    theD0->pt() > dPtCut ) {
          theD0s.push_back( *theD0 );
        }

        if(theD0) delete theD0;
      }
    }
  }
}
// Get methods

const reco::VertexCompositeCandidateCollection& D0Fitter::getD0() const {
  return theD0s;
}

void D0Fitter::resetAll() {
    theD0s.clear();
}
