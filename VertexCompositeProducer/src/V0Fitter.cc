// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      V0Fitter
// 
/**\class V0Fitter V0Fitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/V0Fitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/V0Fitter.h"
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
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <typeinfo>
#include <memory>
#include <vector>

const ParticleMass piMass = 0.13957018;
const double piMassSquared = piMass*piMass;
const ParticleMass protonMass = 0.938272013;
const double protonMassSquared = protonMass*protonMass;
const ParticleMass kaonMass = 0.493677;
const double kaonMassSquared = kaonMass*kaonMass;
const ParticleMass kShortMass = 0.497614;
const ParticleMass phiMass = 1.019445;
const ParticleMass lambdaMass = 1.115683;
const ParticleMass d0Mass = 1.86484;
const ParticleMass dsMass = 1.96847;
const ParticleMass dpmMass = 1.86962;
const ParticleMass lambdaCMass = 2.28646;
const ParticleMass xiMass = 1.32171;
const ParticleMass omegaMass = 1.67245;
float piMass_sigma = piMass*1.e-6;
float protonMass_sigma = protonMass*1.e-6;
float kaonMass_sigma = kaonMass*1.e-6;
float kShortMass_sigma = kShortMass*1.e-6;
float lambdaMass_sigma = lambdaMass*1.e-6;
float d0Mass_sigma = d0Mass*1.e-6;
float dsMass_sigma = dsMass*1.e-6;
float dpmMass_sigma = dpmMass*1.e-6;
float lambdaCMass_sigma = lambdaCMass*1.e-6;
float xiMass_sigma = xiMass*1.e-6;
float omegaMass_sigma = omegaMass*1.e-6;
float phiMass_sigma = phiMass*1.e-6;

// Constructor and (empty) destructor
V0Fitter::V0Fitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) {
//		   const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::ConsumesCollector && iC) {
  using std::string;

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
//  recoAlg = theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm");
//  vtxAlg  = theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm");

  // ------> Initialize parameters from PSet. ALL TRACKED, so no defaults.
  // First set bits to do various things:
  //  -decide whether to use the KVF track smoother, and whether to store those
  //     tracks in the reco::Vertex
  useRefTrax = theParameters.getParameter<bool>(string("useSmoothing"));

  //  -whether to reconstruct K0s
  doKshorts = theParameters.getParameter<bool>(string("selectKshorts"));
  //  -whether to reconstruct Phis
  doPhis = theParameters.getParameter<bool>(string("selectPhis"));
  //  -whether to reconstruct Lambdas
  doLambdas = theParameters.getParameter<bool>(string("selectLambdas"));
  //  -whether to reconstruct Xi 
  doXis = theParameters.getParameter<bool>(string("selectXis"));
  //  -whether to reconstruct Omega
  doOmegas = theParameters.getParameter<bool>(string("selectOmegas"));
  //  -whether to reconstruct D0
  doD0s = theParameters.getParameter<bool>(string("selectD0s"));
  doDSToKsKs = theParameters.getParameter<bool>(string("selectDSToKsKs"));
  doDSToPhiPis = theParameters.getParameter<bool>(string("selectDSToPhiPis"));
  doDPMs = theParameters.getParameter<bool>(string("selectDPMs"));
  //  -whether to reconstruct D0
  doLambdaCToLamPis = theParameters.getParameter<bool>(string("selectLambdaCToLamPis"));
  doLambdaCToKsPs = theParameters.getParameter<bool>(string("selectLambdaCToKsPs"));

  doVertexFit = theParameters.getParameter<bool>(string("doVertexFit"));

  // Second, initialize post-fit cuts
  tkDCACut = theParameters.getParameter<double>(string("tkDCACut"));
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  tkPtCut = theParameters.getParameter<double>(string("tkPtCut"));
  chi2Cut = theParameters.getParameter<double>(string("vtxChi2Cut"));
  rVtxCut = theParameters.getParameter<double>(string("rVtxCut"));
  rVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance2DCut"));
  lVtxCut = theParameters.getParameter<double>(string("lVtxCut"));
  lVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance3DCut"));
  collinCut = theParameters.getParameter<double>(string("collinearityCut"));
  // Sames cuts for Xi
  xiChi2Cut = theParameters.getParameter<double>(string("xiVtxChi2Cut"));
  xiRVtxCut = theParameters.getParameter<double>(string("xiRVtxCut"));
  xiRVtxSigCut = theParameters.getParameter<double>(string("xiVtxSignificance2DCut"));
  xiLVtxCut = theParameters.getParameter<double>(string("xiLVtxCut"));
  xiLVtxSigCut = theParameters.getParameter<double>(string("xiVtxSignificance3DCut"));
  xiCollinCut = theParameters.getParameter<double>(string("xiCollinearityCut"));
  // other common cuts
  kShortMassCut = theParameters.getParameter<double>(string("kShortMassCut"));
  phiMassCut = theParameters.getParameter<double>(string("phiMassCut"));
  lambdaMassCut = theParameters.getParameter<double>(string("lambdaMassCut"));
  d0MassCut = theParameters.getParameter<double>(string("d0MassCut"));
  dsMassCut = theParameters.getParameter<double>(string("dsMassCut"));
  dpmMassCut = theParameters.getParameter<double>(string("dpmMassCut"));
  lambdaCMassCut = theParameters.getParameter<double>(string("lambdaCMassCut"));
  xiMassCut = theParameters.getParameter<double>(string("xiMassCut"));
  omegaMassCut = theParameters.getParameter<double>(string("omegaMassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  batDauTransImpactSigCut = theParameters.getParameter<double>(string("batDauTransImpactSigCut"));
  batDauLongImpactSigCut = theParameters.getParameter<double>(string("batDauLongImpactSigCut"));
  mPiPiCutMin = theParameters.getParameter<double>(string("mPiPiCutMin"));
  mPiPiCutMax = theParameters.getParameter<double>(string("mPiPiCutMax"));
  mKKCutMin = theParameters.getParameter<double>(string("mKKCutMin"));
  mKKCutMax = theParameters.getParameter<double>(string("mKKCutMax"));
  vtxFitter = theParameters.getParameter<edm::InputTag>("vertexFitter");
  innerHitPosCut = theParameters.getParameter<double>(string("innerHitPosCut"));
  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }

  //edm::LogInfo("V0Producer") << "Using " << vtxFitter << " to fit V0 vertices.\n";
  //std::cout << "Using " << vtxFitter << " to fit V0 vertices." << std::endl;
  // FOR DEBUG:
  //initFileOutput();
  //--------------------

  //std::cout << "Entering V0Producer" << std::endl;

//  fitAll(iEvent, iSetup);

  // FOR DEBUG:
  //cleanupFileOutput();
  //--------------------

}

V0Fitter::~V0Fitter() {
}

// Method containing the algorithm for vertex reconstruction
void V0Fitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
//  ESHandle<TrackerGeometry> trackerGeomHandle;
//  ESHandle<GlobalTrackingGeometry> globTkGeomHandle;

  //cout << "Check 0" << endl;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks, theTrackHandle); 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);  
//  iEvent.getByLabel(recoAlg, theTrackHandle);
//  iEvent.getByLabel(vtxAlg,  theVertexHandle);
//  iEvent.getByLabel(std::string("offlineBeamSpot"), theBeamSpotHandle);
  if( !theTrackHandle->size() ) return;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
//  iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeomHandle);
//  iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);

//  trackerGeom = trackerGeomHandle.product();
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
        tmpRef->pt() > tkPtCut ) {
//      TransientTrack tmpTk( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
      TransientTrack tmpTk( *tmpRef, magField );

      math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauLongImpactSig = dzvtx/dzerror;
      double dauTransImpactSig = dxyvtx/dxyerror;
//std::cout<<dxyvtx<<" "<<dauTransImpactSig<<"; "<<dzvtx<<" "<<dauLongImpactSig<<std::endl;
      if( fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut ) {
        theTrackRefs.push_back( tmpRef );
        theTransTracks.push_back( tmpTk );
      }
    }
  }

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
      if(theTrackRefs[trdx1]->charge() < 0. && 
	 theTrackRefs[trdx2]->charge() > 0.) {
	negativeTrackRef = theTrackRefs[trdx1];
	positiveTrackRef = theTrackRefs[trdx2];
	negTransTkPtr = &theTransTracks[trdx1];
	posTransTkPtr = &theTransTracks[trdx2];
      }
      else if(theTrackRefs[trdx1]->charge() > 0. &&
	      theTrackRefs[trdx2]->charge() < 0.) {
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
      if (sqrt( cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y() ) > 120. 
          || std::abs(cxPt.z()) > 300.) continue;

      // Get trajectory states for the tracks at POCA for later cuts
      TrajectoryStateClosestToPoint posTSCP = 
	posTransTkPtr->trajectoryStateClosestToPoint( cxPt );
      TrajectoryStateClosestToPoint negTSCP =
	negTransTkPtr->trajectoryStateClosestToPoint( cxPt );

      if( !posTSCP.isValid() || !negTSCP.isValid() ) continue;

      double totalE = sqrt( posTSCP.momentum().mag2() + piMassSquared ) +
	              sqrt( negTSCP.momentum().mag2() + piMassSquared );
      double totalESq = totalE*totalE;
      double totalPSq =
	( posTSCP.momentum() + negTSCP.momentum() ).mag2();
      double mass = sqrt( totalESq - totalPSq);

      if( mass > mPiPiCutMax || mass < mPiPiCutMin ) continue;

      totalE = sqrt( posTSCP.momentum().mag2() + kaonMassSquared ) +
               sqrt( negTSCP.momentum().mag2() + kaonMassSquared );
      totalESq = totalE*totalE;
      totalPSq =
        ( posTSCP.momentum() + negTSCP.momentum() ).mag2();
      mass = sqrt( totalESq - totalPSq);

      if( mass > mKKCutMax || mass < mKKCutMin ) continue;

      // Create the vertex fitter object and vertex the tracks
      TransientVertex theRecoVertex;
      if(vtxFitter == std::string("KalmanVertexFitter")) {
	KalmanVertexFitter theKalmanFitter(useRefTrax == 0 ? false : true);
        theRecoVertex = theKalmanFitter.vertex(transTracks);
/*
        vector<RefCountedKinematicParticle> dauParticles;
        dauParticles.push_back(pFactory.particle(,piMass,chi,ndf,piMass_sigma));
        dauParticles.push_back(pFactory.particle(protonTT,protonMass,chi,ndf,protonMass_sigma));

        KinematicParticleVertexFitter theKalmanFitter;
        RefCountedKinematicTree theRecoVertex;
        theRecoVertex = theKalmanFitter.fit(dauParticles);
*/
      }
      else if (vtxFitter == std::string("AdaptiveVertexFitter")) {
	useRefTrax = false;
	AdaptiveVertexFitter theAdaptiveFitter;
	theRecoVertex = theAdaptiveFitter.vertex(transTracks);
      }
    
      // Create reco::Vertex object for use in creating the Candidate
      reco::Vertex theVtx;
      // If do not fit vertex, assign a dummy vertex using the beam spot
      if(!doVertexFit) {
        theVtx = reco::Vertex(theBeamSpotHandle->position(), 
	                      theBeamSpotHandle->rotatedCovariance3D(),0.,0.,0);
      }
      // If the vertex is valid, make a VertexCompositeCandidate with it
      else if( theRecoVertex.isValid() && theRecoVertex.totalChiSquared() >= 0. ) {
        theVtx = theRecoVertex;
      }
      else continue;

      // Create and fill vector of refitted TransientTracks
      //  (iff they've been created by the KVF)
      std::vector<TransientTrack> refittedTrax;
      if( theRecoVertex.hasRefittedTracks() ) {
	refittedTrax = theRecoVertex.refittedTracks();
      }

      // Find the vertex d0 and its error

      GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z()); // for V0 candidate vertex
      GlobalPoint vtxPosPrimary(xVtx, yVtx, zVtx); // for primary vertex

      double rVtxMag = 99999.0; 
      double lVtxMag = 99999.0;
      double sigmaRvtxMag = 999.0;
      double sigmaLvtxMag = 999.0;
      double V0Angle = -100.0;

      GlobalVector V0LineOfFlight = GlobalVector (theVtx.x()  - xVtx,
                                                  theVtx.y()  - yVtx,
                                                  theVtx.z()  - zVtx);

      GlobalVector V0GlobalMomentum = posTSCP.momentum() + negTSCP.momentum();

      SMatrixSym3D totalCov;
      if(isVtxPV) totalCov = theVtx.covariance() + vtxPrimary->covariance(); 
      else totalCov = theVtx.covariance() + theBeamSpotHandle->rotatedCovariance3D();

      SVector3 distanceVector3D(V0LineOfFlight.x(), V0LineOfFlight.y(), V0LineOfFlight.z());
      SVector3 distanceVector2D(V0LineOfFlight.x(), V0LineOfFlight.y(), 0.0);
      
      V0Angle = angle(V0LineOfFlight.x(), V0LineOfFlight.y(), V0LineOfFlight.z(),
                      V0GlobalMomentum.x(), V0GlobalMomentum.y(), V0GlobalMomentum.z());
      lVtxMag = V0LineOfFlight.mag();
      rVtxMag = V0LineOfFlight.perp();
      sigmaLvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVector3D)) / lVtxMag;
      sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVector2D)) / rVtxMag;

      // The methods innerOk() and innerPosition() require TrackExtra, which
      // is only available in the RECO data tier, not AOD. Setting innerHitPosCut
      // to -1 avoids this problem and allows to run on AOD.
      if( innerHitPosCut > 0. && positiveTrackRef->innerOk() ) {
	reco::Vertex::Point posTkHitPos = positiveTrackRef->innerPosition();
	double posTkHitPosD2 = 
	  (posTkHitPos.x()-xVtx)*(posTkHitPos.x()-xVtx) +
	  (posTkHitPos.y()-yVtx)*(posTkHitPos.y()-yVtx);
	if( sqrt( posTkHitPosD2 ) < ( rVtxMag - sigmaRvtxMag*innerHitPosCut )
	    ) {
	  continue;
	}
      }
      if( innerHitPosCut > 0. && negativeTrackRef->innerOk() ) {
	reco::Vertex::Point negTkHitPos = negativeTrackRef->innerPosition();
	double negTkHitPosD2 = 
	  (negTkHitPos.x()-xVtx)*(negTkHitPos.x()-xVtx) +
	  (negTkHitPos.y()-yVtx)*(negTkHitPos.y()-yVtx);
	if( sqrt( negTkHitPosD2 ) < ( rVtxMag - sigmaRvtxMag*innerHitPosCut )
	    ) {
	  continue;
	}
      }
      
      if( theVtx.normalizedChi2() > chi2Cut ||
	  rVtxMag < rVtxCut ||
	  rVtxMag / sigmaRvtxMag < rVtxSigCut ||
          lVtxMag < lVtxCut ||
          lVtxMag / sigmaLvtxMag < lVtxSigCut ||
          cos(V0Angle) < collinCut
      )	continue;

      // Cuts finished, now we create the candidates and push them back into the collections.
      
      std::auto_ptr<TrajectoryStateClosestToPoint> trajPlus;
      std::auto_ptr<TrajectoryStateClosestToPoint> trajMins;

      if( useRefTrax && refittedTrax.size() > 1 ) {
	// Need an iterator over the refitted tracks for below
	std::vector<TransientTrack>::iterator traxIter = refittedTrax.begin(),
	  traxEnd = refittedTrax.end();

	// TransientTrack objects to hold the positive and negative
	//  refitted tracks
	TransientTrack* thePositiveRefTrack = 0;
	TransientTrack* theNegativeRefTrack = 0;
        
	for( ; traxIter != traxEnd; ++traxIter) {
	  if( traxIter->track().charge() > 0. ) {
	    thePositiveRefTrack = &*traxIter;
	  }
	  else if (traxIter->track().charge() < 0.) {
	    theNegativeRefTrack = &*traxIter;
	  }
	}
        if (thePositiveRefTrack == 0 || theNegativeRefTrack == 0) continue;
	trajPlus.reset(new TrajectoryStateClosestToPoint(thePositiveRefTrack->trajectoryStateClosestToPoint(vtxPos)));
	trajMins.reset(new TrajectoryStateClosestToPoint(theNegativeRefTrack->trajectoryStateClosestToPoint(vtxPos)));
      }
      else {
	trajPlus.reset(new TrajectoryStateClosestToPoint(posTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
	trajMins.reset(new TrajectoryStateClosestToPoint(negTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));

      }

      if( trajPlus.get() == 0 || trajMins.get() == 0 || !trajPlus->isValid() || !trajMins->isValid() ) continue;

      posTransTkPtr = negTransTkPtr = 0;

      GlobalVector positiveP(trajPlus->momentum());
      GlobalVector negativeP(trajMins->momentum());
      GlobalVector totalP(positiveP + negativeP);

      //cleanup stuff we don't need anymore
      trajPlus.reset();
      trajMins.reset();

      // calculate total energy of V0 3 ways:
      //  Assume it's a kShort, a Lambda, or a LambdaBar.
      double piPlusE = sqrt( positiveP.mag2() + piMassSquared );
      double piMinusE = sqrt( negativeP.mag2() + piMassSquared );
      double kaonPlusE = sqrt( positiveP.mag2() + kaonMassSquared );
      double kaonMinusE = sqrt( negativeP.mag2() + kaonMassSquared );
      double protonE = sqrt( positiveP.mag2() + protonMassSquared );
      double antiProtonE = sqrt( negativeP.mag2() + protonMassSquared );
      double d0ETot = kaonMinusE + piPlusE;
      double d0BarETot = kaonPlusE + piMinusE;
      double kShortETot = piPlusE + piMinusE;
      double phiETot = kaonPlusE + kaonMinusE;
      double lambdaEtot = protonE + piMinusE;
      double lambdaBarEtot = antiProtonE + piPlusE;

      using namespace reco;

      // Create momentum 4-vectors for the 3 candidate types
      const Particle::LorentzVector kShortP4(totalP.x(), 
					     totalP.y(), totalP.z(), 
					     kShortETot);
      const Particle::LorentzVector phiP4(totalP.x(),
                                             totalP.y(), totalP.z(),
                                             phiETot);
      const Particle::LorentzVector lambdaP4(totalP.x(), 
					     totalP.y(), totalP.z(), 
					     lambdaEtot);
      const Particle::LorentzVector lambdaBarP4(totalP.x(), 
						totalP.y(), totalP.z(), 
						lambdaBarEtot);
      const Particle::LorentzVector d0P4(totalP.x(),
                                             totalP.y(), totalP.z(),
                                             d0ETot);
      const Particle::LorentzVector d0BarP4(totalP.x(),
                                             totalP.y(), totalP.z(),
                                             d0BarETot);

      Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
      const Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
      double vtxChi2(theVtx.chi2());
      double vtxNdof(theVtx.ndof());

      // Create the VertexCompositeCandidate object that will be stored in the Event
      VertexCompositeCandidate* theKshort = 0;
      VertexCompositeCandidate* thePhi = 0;
      VertexCompositeCandidate* theLambda = 0;
      VertexCompositeCandidate* theLambdaBar = 0;
      VertexCompositeCandidate* theD0 = 0;
      VertexCompositeCandidate* theD0Bar = 0;

      if( doKshorts ) {
	theKshort = new VertexCompositeCandidate(0, kShortP4, vtx, vtxCov, vtxChi2, vtxNdof);
      }

      if( doPhis ) {
        thePhi = new VertexCompositeCandidate(0, phiP4, vtx, vtxCov, vtxChi2, vtxNdof);
      }

      if( doLambdas ) {
	if( positiveP.mag() > negativeP.mag() ) {
	  theLambda = 
	    new VertexCompositeCandidate(0, lambdaP4, vtx, vtxCov, vtxChi2, vtxNdof);
	}
	else {
	  theLambdaBar = 
	    new VertexCompositeCandidate(0, lambdaBarP4, vtx, vtxCov, vtxChi2, vtxNdof);
	}
      }
      if( doD0s ) {
          theD0 =
            new VertexCompositeCandidate(0, d0P4, vtx, vtxCov, vtxChi2, vtxNdof);
          theD0Bar =
            new VertexCompositeCandidate(0, d0BarP4, vtx, vtxCov, vtxChi2, vtxNdof);
      }

      // Create daughter candidates for the VertexCompositeCandidates
      RecoChargedCandidate 
	thePiPlusCand(1, Particle::LorentzVector(positiveP.x(), 
						 positiveP.y(), positiveP.z(),
						 piPlusE), vtx);
      thePiPlusCand.setTrack(positiveTrackRef);
      
      RecoChargedCandidate
	thePiMinusCand(-1, Particle::LorentzVector(negativeP.x(), 
						   negativeP.y(), negativeP.z(),
						   piMinusE), vtx);
      thePiMinusCand.setTrack(negativeTrackRef);
      
       RecoChargedCandidate
        theKaonPlusCand(1, Particle::LorentzVector(positiveP.x(),
                                                 positiveP.y(), positiveP.z(),
                                                 kaonPlusE), vtx);
      theKaonPlusCand.setTrack(positiveTrackRef);

      RecoChargedCandidate
        theKaonMinusCand(-1, Particle::LorentzVector(negativeP.x(),
                                                   negativeP.y(), negativeP.z(),
                                                  kaonMinusE), vtx);
      theKaonMinusCand.setTrack(negativeTrackRef);

      RecoChargedCandidate
	theProtonCand(1, Particle::LorentzVector(positiveP.x(),
						 positiveP.y(), positiveP.z(),
						 protonE), vtx);
      theProtonCand.setTrack(positiveTrackRef);

      RecoChargedCandidate
	theAntiProtonCand(-1, Particle::LorentzVector(negativeP.x(),
						      negativeP.y(), negativeP.z(),
						      antiProtonE), vtx);
      theAntiProtonCand.setTrack(negativeTrackRef);

      AddFourMomenta addp4;
      // Store the daughter Candidates in the VertexCompositeCandidates 
      //    if they pass mass cuts
      if( doKshorts ) {
	theKshort->addDaughter(thePiPlusCand);
	theKshort->addDaughter(thePiMinusCand);
	theKshort->setPdgId(310);
	addp4.set( *theKshort );
	if( theKshort->mass() < kShortMass + kShortMassCut &&
	    theKshort->mass() > kShortMass - kShortMassCut ) {
	  theKshorts.push_back( *theKshort );
	}
      }
      
      if( doPhis ) {
        thePhi->addDaughter(theKaonPlusCand);
        thePhi->addDaughter(theKaonMinusCand);
        thePhi->setPdgId(333);
        addp4.set( *thePhi );
        if( thePhi->mass() < phiMass + phiMassCut &&
            thePhi->mass() > phiMass - phiMassCut ) {
          thePhis.push_back( *thePhi );
        }
      }

      if( doLambdas && theLambda ) {
	theLambda->addDaughter(theProtonCand);
	theLambda->addDaughter(thePiMinusCand);
	theLambda->setPdgId(3122);
	addp4.set( *theLambda );
	if( theLambda->mass() < lambdaMass + lambdaMassCut &&
	    theLambda->mass() > lambdaMass - lambdaMassCut ) {
	  theLambdas.push_back( *theLambda );
	}
      }
      else if ( doLambdas && theLambdaBar ) {
	theLambdaBar->addDaughter(theAntiProtonCand);
	theLambdaBar->addDaughter(thePiPlusCand);
	theLambdaBar->setPdgId(-3122);
	addp4.set( *theLambdaBar );
	if( theLambdaBar->mass() < lambdaMass + lambdaMassCut &&
	    theLambdaBar->mass() > lambdaMass - lambdaMassCut ) {
	  theLambdas.push_back( *theLambdaBar );
	}
      }

      if( doD0s && theD0 ) {
//std::cout<<"do D0"<<std::endl;
        theD0->addDaughter(theKaonMinusCand);
        theD0->addDaughter(thePiPlusCand);
        theD0->setPdgId(421);
        addp4.set( *theD0 );
//std::cout<<"D0 mass="<<theD0->mass()<<std::endl;
        if( theD0->mass() < d0Mass + d0MassCut &&
            theD0->mass() > d0Mass - d0MassCut ) {
          theD0s.push_back( *theD0 );
//std::cout<<"add D0"<<std::endl;
        }
      }
      else if ( doD0s && theD0Bar ) {
        theD0Bar->addDaughter(theKaonPlusCand);
        theD0Bar->addDaughter(thePiMinusCand);
        theD0Bar->setPdgId(-421);
        addp4.set( *theD0Bar );
        if( theD0Bar->mass() < d0Mass + d0MassCut &&
            theD0Bar->mass() > d0Mass - d0MassCut ) {
          theD0s.push_back( *theD0Bar );
        }
      }

      if(theKshort) delete theKshort;
      if(thePhi) delete thePhi;
      if(theLambda) delete theLambda;
      if(theLambdaBar) delete theLambdaBar;
      if(theD0) delete theD0;
      if(theD0Bar) delete theD0Bar;
      theKshort = thePhi = theLambda = theLambdaBar = theD0 = theD0Bar = 0;
    }
  }

  if((doLambdaCToKsPs || doDSToKsKs || doDPMs) && theKshorts.size() > 0) 
  {
    for(unsigned it=0; it<theKshorts.size(); ++it){

      const reco::VertexCompositeCandidate & theKshort = theKshorts[it];

      float massWindow = 0.040;
      if(theKshort.mass() > kShortMass + massWindow || theKshort.mass() < kShortMass - massWindow) continue;

      vector<RecoChargedCandidate> v0daughters;
      vector<TrackRef> theDaughterTracks;

      v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                              (theKshort.daughter(0))) );
      v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                              (theKshort.daughter(1))) );

      for(unsigned int j = 0; j < v0daughters.size(); ++j) theDaughterTracks.push_back(v0daughters[j].track());

      vector<TransientTrack> tracksForKalmanFit;
      for (unsigned int ndx = 0; ndx < theDaughterTracks.size(); ndx++) {
          tracksForKalmanFit.push_back(TransientTrack(theDaughterTracks[ndx], &(*bFieldHandle)));
      }

      for(unsigned int trdx = 0; trdx < theTrackRefs.size(); trdx++) {

         if ( theTrackRefs[trdx].isNull() ) continue;

         float batPionMass = protonMass;
         float batPionMass_sigma = protonMass_sigma;

         if(doLambdaCToKsPs)
         {
           batPionMass = protonMass;
           batPionMass_sigma = protonMass_sigma;
         }

         if(doDSToKsKs)
         {
           batPionMass = kaonMass;
           batPionMass_sigma = kaonMass_sigma;
         }

        if(doDPMs)
         {
           batPionMass = piMass;
           batPionMass_sigma = piMass_sigma;
         }

         double totalE1 = theKshort.energy()+sqrt(theTrackRefs[trdx]->p()*theTrackRefs[trdx]->p()+batPionMass*batPionMass);
         double totalE1Sq = totalE1*totalE1;
         double totalPSq = ( theTrackRefs[trdx]->momentum() + theKshort.momentum() ).mag2();
         double mass = sqrt( totalE1Sq - totalPSq);

         if ( doDSToKsKs && (mass > dsMass + dsMassCut*1.3 ||  mass < dsMass - dsMassCut*1.3)) continue;
         if ( doDPMs && (mass > dpmMass + dpmMassCut*1.3 || mass < dpmMass - dpmMassCut*1.3)) continue;
         if ( doLambdaCToKsPs && (mass > lambdaCMass + lambdaCMassCut*1.3 ||  mass < lambdaCMass - lambdaCMassCut*1.3)) continue;

         math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
         double dzvtx = theTrackRefs[trdx]->dz(bestvtx);
         double dxyvtx = theTrackRefs[trdx]->dxy(bestvtx);
         if(fabs(dzvtx)>50 || fabs(dxyvtx)>50) continue;

         double dzerror = sqrt(theTrackRefs[trdx]->dzError()*theTrackRefs[trdx]->dzError()+zVtxError*zVtxError);
         double dxyerror = sqrt(theTrackRefs[trdx]->d0Error()*theTrackRefs[trdx]->d0Error()+xVtxError*yVtxError);
         double dauLongImpactSig = dzvtx/dzerror;
         double dauTransImpactSig = dxyvtx/dxyerror;
         if( fabs(dauTransImpactSig) < batDauTransImpactSigCut || fabs(dauLongImpactSig) < batDauLongImpactSigCut ) continue;

         bool match = false;

         for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
           if (theTrackRefs[trdx]->charge() == theDaughterTracks[j]->charge() &&
               theTrackRefs[trdx]->momentum() == theDaughterTracks[j]->momentum() ) match = true;
           if (match) break;
         }
         if (match) continue; // Track is already used in making the V0

         TransientTrack pion1TT(theDaughterTracks[0], &(*bFieldHandle) );
         TransientTrack pion2TT(theDaughterTracks[1], &(*bFieldHandle) );
         TransientTrack batPionTT(theTrackRefs[trdx], &(*bFieldHandle) );

         if (!batPionTT.isValid()) continue;
         if (!pion1TT.isValid()) continue;
         if (!pion2TT.isValid()) continue;

         //Creating a KinematicParticleFactory
         KinematicParticleFactoryFromTransientTrack pFactory;

         float chi = 0.;
         float ndf = 0.;
         vector<RefCountedKinematicParticle> ksParticles;
         ksParticles.push_back(pFactory.particle(pion1TT,piMass,chi,ndf,piMass_sigma));
         ksParticles.push_back(pFactory.particle(pion2TT,piMass,chi,ndf,piMass_sigma));

         KinematicParticleVertexFitter fitter;
         RefCountedKinematicTree ksVertexFitTree;
         ksVertexFitTree = fitter.fit(ksParticles);
         if (!ksVertexFitTree->isValid()) continue;

         ksVertexFitTree->movePointerToTheTop();

         RefCountedKinematicParticle ks_vFit = ksVertexFitTree->currentParticle();
         RefCountedKinematicVertex ks_vFit_vertex = ksVertexFitTree->currentDecayVertex();

         ksVertexFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle ksPi1 = ksVertexFitTree->currentParticle();
         ksVertexFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle ksPi2 = ksVertexFitTree->currentParticle();

         KinematicParticleFitter csFitterKs;
         KinematicConstraint * ks_c = new MassKinematicConstraint(kShortMass,kShortMass_sigma);

         ksVertexFitTree->movePointerToTheTop();
         ksVertexFitTree = csFitterKs.fit(ks_c,ksVertexFitTree);
         if (!ksVertexFitTree->isValid()) continue;
         ksVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle ks_vFit_withMC = ksVertexFitTree->currentParticle();

         if (!ks_vFit_withMC->currentState().isValid()) continue;
/*
         float batPionMass = protonMass;
         float batPionMass_sigma = protonMass_sigma;

         if(doLambdaCToKsPs)
         {
           batPionMass = protonMass;
           batPionMass_sigma = protonMass_sigma;
         }

         if(doDSToKsKs)
         {
           batPionMass = kaonMass;
           batPionMass_sigma = kaonMass_sigma;
         }

        if(doDPMs)
         {
           batPionMass = piMass;
           batPionMass_sigma = piMass_sigma;
         }
*/
         vector<RefCountedKinematicParticle> xiFitParticles;

         xiFitParticles.push_back(pFactory.particle(batPionTT,batPionMass,chi,ndf,batPionMass_sigma));
         xiFitParticles.push_back(ks_vFit_withMC);

         RefCountedKinematicTree xiFitTree = fitter.fit(xiFitParticles);
         if (!xiFitTree->isValid()) continue;

         xiFitTree->movePointerToTheTop();
         RefCountedKinematicParticle xiCand = xiFitTree->currentParticle();
         if (!xiCand->currentState().isValid()) continue;

         RefCountedKinematicVertex xiDecayVertex = xiFitTree->currentDecayVertex();
         if (!xiDecayVertex->vertexIsValid()) continue;

//         if ( doDSToKsKs && (xiCand->currentState().mass() > dsMass + dsMassCut*1.5 ||  xiCand->currentState().mass() < dsMass - dsMassCut*1.5)) continue;
//         if ( doDPMs && (xiCand->currentState().mass() > dpmMass + dpmMassCut*1.5 ||  xiCand->currentState().mass() < dpmMass - dpmMassCut*1.5)) continue;
//         if ( doLambdaCToKsPs && (xiCand->currentState().mass() > lambdaCMass + lambdaCMassCut*1.5 ||  xiCand->currentState().mass() < lambdaCMass - lambdaCMassCut*1.5)) continue;

         xiFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle batPionCand = xiFitTree->currentParticle();
         xiFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle ksCandMC = xiFitTree->currentParticle();

         if(!batPionCand->currentState().isValid() || !ksCandMC->currentState().isValid()) continue;

         KinematicParameters batPionKP = batPionCand->currentState().kinematicParameters();
         KinematicParameters ksCandKP = ksCandMC->currentState().kinematicParameters();

         GlobalVector xiTotalP = GlobalVector (xiCand->currentState().globalMomentum().x(),
                                               xiCand->currentState().globalMomentum().y(),
                                               xiCand->currentState().globalMomentum().z());

         GlobalVector batPionTotalP = GlobalVector(batPionKP.momentum().x(),batPionKP.momentum().y(),batPionKP.momentum().z());
         GlobalVector ksTotalP = GlobalVector(ksCandKP.momentum().x(),ksCandKP.momentum().y(),ksCandKP.momentum().z());

         double batPionTotalE = sqrt( batPionTotalP.mag2() + batPionMass*batPionMass );
         double ksTotalE = sqrt( ksTotalP.mag2() + kShortMass*kShortMass );
         double xiTotalE = batPionTotalE + ksTotalE;

         const Particle::LorentzVector xiP4(xiTotalP.x(), xiTotalP.y(), xiTotalP.z(), xiTotalE);

         Particle::Point xiVtx((*xiDecayVertex).position().x(), (*xiDecayVertex).position().y(), (*xiDecayVertex).position().z());
         std::vector<double> xiVtxEVec;
         xiVtxEVec.push_back( ks_vFit_vertex->error().cxx() );
         xiVtxEVec.push_back( ks_vFit_vertex->error().cyx() );
         xiVtxEVec.push_back( ks_vFit_vertex->error().cyy() );
         xiVtxEVec.push_back( ks_vFit_vertex->error().czx() );
         xiVtxEVec.push_back( ks_vFit_vertex->error().czy() );
         xiVtxEVec.push_back( ks_vFit_vertex->error().czz() );
         SMatrixSym3D xiVtxCovMatrix(xiVtxEVec.begin(), xiVtxEVec.end());
         const Vertex::CovarianceMatrix xiVtxCov(xiVtxCovMatrix);
         double xiVtxChi2(xiDecayVertex->chiSquared());
         double xiVtxNdof(xiDecayVertex->degreesOfFreedom());
         double xiNormalizedChi2 = xiVtxChi2/xiVtxNdof;

         double xiRVtxMag = 99999.0;
         double xiLVtxMag = 99999.0;
         double xiSigmaRvtxMag = 999.0;
         double xiSigmaLvtxMag = 999.0;
         double xiV0Angle = -100.0;

         GlobalVector xiV0LineOfFlight = GlobalVector (xiVtx.x() - xVtx,
                                                       xiVtx.y() - yVtx,
                                                       xiVtx.z() - zVtx);

         SMatrixSym3D xiTotalCov;
         if(isVtxPV) xiTotalCov = xiVtxCovMatrix + vtxPrimary->covariance();
         else xiTotalCov = xiVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

         SVector3 distanceVector3D(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), xiV0LineOfFlight.z());
         SVector3 distanceVector2D(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), 0.0);

         xiV0Angle = angle(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), xiV0LineOfFlight.z(),
                           xiTotalP.x(), xiTotalP.y(), xiTotalP.z());
         xiLVtxMag = xiV0LineOfFlight.mag();
         xiRVtxMag = xiV0LineOfFlight.perp();
         xiSigmaLvtxMag = sqrt(ROOT::Math::Similarity(xiTotalCov, distanceVector3D)) / xiLVtxMag;
         xiSigmaRvtxMag = sqrt(ROOT::Math::Similarity(xiTotalCov, distanceVector2D)) / xiRVtxMag;

         if( xiNormalizedChi2 > xiChi2Cut ||
             xiRVtxMag < xiRVtxCut ||
             xiRVtxMag / xiSigmaRvtxMag < xiRVtxSigCut ||
             xiLVtxMag < xiLVtxCut ||
             xiLVtxMag / xiSigmaLvtxMag < xiLVtxSigCut ||
             cos(xiV0Angle) < xiCollinCut
         ) continue;

         RecoChargedCandidate PionCand(theTrackRefs[trdx]->charge(), Particle::LorentzVector(batPionTotalP.x(),
                                                 batPionTotalP.y(), batPionTotalP.z(),
                                                 batPionTotalE), xiVtx);
         PionCand.setTrack(theTrackRefs[trdx]);
         AddFourMomenta addp4;

         if(doLambdaCToKsPs)
         {
           VertexCompositeCandidate* theLambdaCToKsP = 0;
           theLambdaCToKsP = new VertexCompositeCandidate(0, xiP4, xiVtx, xiVtxCov, xiVtxChi2, xiVtxNdof);
           theLambdaCToKsP->addDaughter(theKshort);
           theLambdaCToKsP->addDaughter(PionCand);
           if(theTrackRefs[trdx]->charge()>0) { theLambdaCToKsP->setPdgId(4122); theLambdaCToKsP->setCharge(1); }
           else { theLambdaCToKsP->setPdgId(-4122); theLambdaCToKsP->setCharge(-1); }
           addp4.set( *theLambdaCToKsP );
           if( theLambdaCToKsP->mass() < lambdaCMass + lambdaCMassCut &&
               theLambdaCToKsP->mass() > lambdaCMass - lambdaCMassCut ) theLambdaCToKsPs.push_back( *theLambdaCToKsP );
           if(theLambdaCToKsP) delete theLambdaCToKsP;
           theLambdaCToKsP = 0;
         }

         if(doDSToKsKs)
         {
           VertexCompositeCandidate* theDSToKsK = 0;
           theDSToKsK = new VertexCompositeCandidate(0, xiP4, xiVtx, xiVtxCov, xiVtxChi2, xiVtxNdof);
           theDSToKsK->addDaughter(theKshort);
           theDSToKsK->addDaughter(PionCand);
           if(theTrackRefs[trdx]->charge()>0) { theDSToKsK->setPdgId(431); theDSToKsK->setCharge(1); }
           else { theDSToKsK->setPdgId(-431); theDSToKsK->setCharge(-1); }
           addp4.set( *theDSToKsK );
           if( theDSToKsK->mass() < dsMass + dsMassCut &&
               theDSToKsK->mass() > dsMass - dsMassCut ) theDSToKsKs.push_back( *theDSToKsK );
           if(theDSToKsK) delete theDSToKsK;
           theDSToKsK = 0;
         }

         if(doDPMs)
         {
           VertexCompositeCandidate* theDPM = 0;
           theDPM = new VertexCompositeCandidate(0, xiP4, xiVtx, xiVtxCov, xiVtxChi2, xiVtxNdof);
           theDPM->addDaughter(theKshort);
           theDPM->addDaughter(PionCand);
           if(theTrackRefs[trdx]->charge()>0) { theDPM->setPdgId(411); theDPM->setCharge(1); }
           else { theDPM->setPdgId(-431); theDPM->setCharge(-1); }
           addp4.set( *theDPM );
           if( theDPM->mass() < dpmMass + dpmMassCut &&
               theDPM->mass() > dpmMass - dpmMassCut ) theDPMs.push_back( *theDPM );
           if(theDPM) delete theDPM;
           theDPM = 0;
         }

      }
    }

  };

  // DsToPhiPi reconstruction                                                                                                                                                                                                        
  if( doDSToPhiPis && thePhis.size() > 0) 
    {
      for(unsigned it=0; it<thePhis.size(); ++it){

	const reco::VertexCompositeCandidate & thePhi = thePhis[it];

	float massWindow = 0.010;
	if(thePhi.mass() > phiMass + massWindow || thePhi.mass() < phiMass - massWindow) continue;

	vector<RecoChargedCandidate> v0daughters;
	vector<TrackRef> theDaughterTracks;

	v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
				 (thePhi.daughter(0))) );
	v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
				 (thePhi.daughter(1))) );

	for(unsigned int j = 0; j < v0daughters.size(); ++j) theDaughterTracks.push_back(v0daughters[j].track());

	vector<TransientTrack> tracksForKalmanFit;
	for (unsigned int ndx = 0; ndx < theDaughterTracks.size(); ndx++) {
          tracksForKalmanFit.push_back(TransientTrack(theDaughterTracks[ndx], &(*bFieldHandle)));
	}

	for(unsigned int trdx = 0; trdx < theTrackRefs.size(); trdx++) {

	  if ( theTrackRefs[trdx].isNull() ) continue;

	  math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
	  double dzvtx = theTrackRefs[trdx]->dz(bestvtx);
	  double dxyvtx = theTrackRefs[trdx]->dxy(bestvtx);
	  if(fabs(dzvtx)>50 || fabs(dxyvtx)>50) continue;

	  double dzerror = sqrt(theTrackRefs[trdx]->dzError()*theTrackRefs[trdx]->dzError()+zVtxError*zVtxError);
	  double dxyerror = sqrt(theTrackRefs[trdx]->d0Error()*theTrackRefs[trdx]->d0Error()+xVtxError*yVtxError);
	  double dauLongImpactSig = dzvtx/dzerror;
	  double dauTransImpactSig = dxyvtx/dxyerror;
	  if( fabs(dauTransImpactSig) < batDauTransImpactSigCut || fabs(dauLongImpactSig) < batDauLongImpactSigCut ) continue;

	  bool match = false;

	  for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
	    if (theTrackRefs[trdx]->charge() == theDaughterTracks[j]->charge() &&
		theTrackRefs[trdx]->momentum() == theDaughterTracks[j]->momentum() ) match = true;
	    if (match) break;
	  }
	  if (match) continue; // Track is already used in making the V0

	  TransientTrack pion1TT(theDaughterTracks[0], &(*bFieldHandle) );
	  TransientTrack pion2TT(theDaughterTracks[1], &(*bFieldHandle) );
	  TransientTrack batPionTT(theTrackRefs[trdx], &(*bFieldHandle) );

	  if (!batPionTT.isValid()) continue;
	  if (!pion1TT.isValid()) continue;
	  if (!pion2TT.isValid()) continue;

	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;

	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> phiParticles;
	  phiParticles.push_back(pFactory.particle(pion1TT,kaonMass,chi,ndf,kaonMass_sigma));
	  phiParticles.push_back(pFactory.particle(pion2TT,kaonMass,chi,ndf,kaonMass_sigma));

	  KinematicParticleVertexFitter fitter;
	  RefCountedKinematicTree phiVertexFitTree;
	  phiVertexFitTree = fitter.fit(phiParticles);
	  if (!phiVertexFitTree->isValid()) continue;

	  phiVertexFitTree->movePointerToTheTop();

	  RefCountedKinematicParticle phi_vFit = phiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex phi_vFit_vertex = phiVertexFitTree->currentDecayVertex();

	  phiVertexFitTree->movePointerToTheFirstChild();
	  RefCountedKinematicParticle phiPi1 = phiVertexFitTree->currentParticle();
	  phiVertexFitTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle phiPi2 = phiVertexFitTree->currentParticle();

	  KinematicParticleFitter csFitterPhi;
	  KinematicConstraint * phi_c = new MassKinematicConstraint(phiMass,phiMass_sigma);

	  phiVertexFitTree->movePointerToTheTop();
	  phiVertexFitTree = csFitterPhi.fit(phi_c,phiVertexFitTree);
	  if (!phiVertexFitTree->isValid()) continue;
	  phiVertexFitTree->movePointerToTheTop();
	  RefCountedKinematicParticle phi_vFit_withMC = phiVertexFitTree->currentParticle();

	  if (!phi_vFit_withMC->currentState().isValid()) continue;

	  float batPionMass = piMass;
	  float batPionMass_sigma = piMass_sigma;

         vector<RefCountedKinematicParticle> xiFitParticles;

         xiFitParticles.push_back(pFactory.particle(batPionTT,batPionMass,chi,ndf,batPionMass_sigma));
         xiFitParticles.push_back(phi_vFit_withMC);

         RefCountedKinematicTree xiFitTree = fitter.fit(xiFitParticles);
         if (!xiFitTree->isValid()) continue;

         xiFitTree->movePointerToTheTop();
         RefCountedKinematicParticle xiCand = xiFitTree->currentParticle();
         if (!xiCand->currentState().isValid()) continue;

         RefCountedKinematicVertex xiDecayVertex = xiFitTree->currentDecayVertex();
         if (!xiDecayVertex->vertexIsValid()) continue;

         if ( doDSToPhiPis && (xiCand->currentState().mass() > dsMass + dsMassCut*1.5 ||  xiCand->currentState().mass() < dsMass - dsMassCut*1.5)) continue;

         xiFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle batPionCand = xiFitTree->currentParticle();
         xiFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle phiCandMC = xiFitTree->currentParticle();

	if(!batPionCand->currentState().isValid() || !phiCandMC->currentState().isValid()) continue;

         KinematicParameters batPionKP = batPionCand->currentState().kinematicParameters();
         KinematicParameters phiCandKP = phiCandMC->currentState().kinematicParameters();

         GlobalVector xiTotalP = GlobalVector (xiCand->currentState().globalMomentum().x(),
                                               xiCand->currentState().globalMomentum().y(),
                                               xiCand->currentState().globalMomentum().z());

         GlobalVector batPionTotalP = GlobalVector(batPionKP.momentum().x(),batPionKP.momentum().y(),batPionKP.momentum().z());
         GlobalVector phiTotalP = GlobalVector(phiCandKP.momentum().x(),phiCandKP.momentum().y(),phiCandKP.momentum().z());

         double batPionTotalE = sqrt( batPionTotalP.mag2() + batPionMass*batPionMass );
         double phiTotalE = sqrt( phiTotalP.mag2() + phiMass*phiMass );
         double xiTotalE = batPionTotalE + phiTotalE;

         const Particle::LorentzVector xiP4(xiTotalP.x(), xiTotalP.y(), xiTotalP.z(), xiTotalE);

         Particle::Point xiVtx((*xiDecayVertex).position().x(), (*xiDecayVertex).position().y(), (*xiDecayVertex).position().z());
         std::vector<double> xiVtxEVec;
         xiVtxEVec.push_back( phi_vFit_vertex->error().cxx() );
         xiVtxEVec.push_back( phi_vFit_vertex->error().cyx() );
         xiVtxEVec.push_back( phi_vFit_vertex->error().cyy() );
         xiVtxEVec.push_back( phi_vFit_vertex->error().czx() );
         xiVtxEVec.push_back( phi_vFit_vertex->error().czy() );
         xiVtxEVec.push_back( phi_vFit_vertex->error().czz() );
         SMatrixSym3D xiVtxCovMatrix(xiVtxEVec.begin(), xiVtxEVec.end());
         const Vertex::CovarianceMatrix xiVtxCov(xiVtxCovMatrix);
         double xiVtxChi2(xiDecayVertex->chiSquared());
         double xiVtxNdof(xiDecayVertex->degreesOfFreedom());
         double xiNormalizedChi2 = xiVtxChi2/xiVtxNdof;

         double xiRVtxMag = 99999.0;
         double xiLVtxMag = 99999.0;
         double xiSigmaRvtxMag = 999.0;
         double xiSigmaLvtxMag = 999.0;
         double xiV0Angle = -100.0;

         GlobalVector xiV0LineOfFlight = GlobalVector (xiVtx.x() - xVtx,
                                                       xiVtx.y() - yVtx,
                                                       xiVtx.z() - zVtx);

         SMatrixSym3D xiTotalCov;
         if(isVtxPV) xiTotalCov = xiVtxCovMatrix + vtxPrimary->covariance();
         else xiTotalCov = xiVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

         SVector3 distanceVector3D(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), xiV0LineOfFlight.z());
         SVector3 distanceVector2D(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), 0.0);

         xiV0Angle = angle(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), xiV0LineOfFlight.z(),
                           xiTotalP.x(), xiTotalP.y(), xiTotalP.z());
         xiLVtxMag = xiV0LineOfFlight.mag();
         xiRVtxMag = xiV0LineOfFlight.perp();
         xiSigmaLvtxMag = sqrt(ROOT::Math::Similarity(xiTotalCov, distanceVector3D)) / xiLVtxMag;
         xiSigmaRvtxMag = sqrt(ROOT::Math::Similarity(xiTotalCov, distanceVector2D)) / xiRVtxMag;

         if( xiNormalizedChi2 > xiChi2Cut ||
             xiRVtxMag < xiRVtxCut ||
             xiRVtxMag / xiSigmaRvtxMag < xiRVtxSigCut ||
             xiLVtxMag < xiLVtxCut ||
             xiLVtxMag / xiSigmaLvtxMag < xiLVtxSigCut ||
             cos(xiV0Angle) < xiCollinCut
         ) continue;

         RecoChargedCandidate PionCand(theTrackRefs[trdx]->charge(), Particle::LorentzVector(batPionTotalP.x(),
                                                 batPionTotalP.y(), batPionTotalP.z(),
                                                 batPionTotalE), xiVtx);
         PionCand.setTrack(theTrackRefs[trdx]);
         AddFourMomenta addp4;

         if(doDSToPhiPis)
         {
           VertexCompositeCandidate* theDSToPhiPi = 0;
           theDSToPhiPi = new VertexCompositeCandidate(0, xiP4, xiVtx, xiVtxCov, xiVtxChi2, xiVtxNdof);
           theDSToPhiPi->addDaughter(thePhi);
           theDSToPhiPi->addDaughter(PionCand);
           if(theTrackRefs[trdx]->charge()>0) { theDSToPhiPi->setPdgId(431); theDSToPhiPi->setCharge(1); }
           else { theDSToPhiPi->setPdgId(-431); theDSToPhiPi->setCharge(-1); }
           addp4.set( *theDSToPhiPi );
           if( theDSToPhiPi->mass() < dsMass + dsMassCut &&
               theDSToPhiPi->mass() > dsMass - dsMassCut ) theDSToPhiPis.push_back( *theDSToPhiPi );
           if(theDSToPhiPi) delete theDSToPhiPi;
           theDSToPhiPi = 0;
         }
      }
    }
  };

  // Xi reconstruction
  if((doXis || doOmegas || doLambdaCToLamPis) && theLambdas.size() > 0)
  {

    for(unsigned it=0; it<theLambdas.size(); ++it){
        
      const reco::VertexCompositeCandidate & theLambda = theLambdas[it];

      // check lambda mass to be within +- 20MeV...
      float massWindow = 0.020;
      if(theLambda.mass() > lambdaMass + massWindow || theLambda.mass() < lambdaMass - massWindow) continue;
      // after finding values also cut on xi vertex probability > 0.0001 and batPion Dxyz/sigma > 0.5

      int lamPDG = theLambda.pdgId();
      bool lamIsParticle = (lamPDG > 0);

      //get daughter tracks from V0 candidate
      vector<RecoChargedCandidate> v0daughters;
      vector<TrackRef> theDaughterTracks;

      //check momentum and add pion first, proton second
      if (theLambda.daughter(0)->momentum().mag2() < theLambda.daughter(1)->momentum().mag2()){
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theLambda.daughter(0))) );
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theLambda.daughter(1))) );
       } else {
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theLambda.daughter(1))) );
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theLambda.daughter(0))) );
       }
       for(unsigned int j = 0; j < v0daughters.size(); ++j) {
         theDaughterTracks.push_back(v0daughters[j].track());
       }

       vector<TransientTrack> tracksForKalmanFit;
       for (unsigned int ndx = 0; ndx < theDaughterTracks.size(); ndx++) {
           tracksForKalmanFit.push_back(TransientTrack(theDaughterTracks[ndx], &(*bFieldHandle)));
       }

       vector<double> lamDauMasses;
       if (lamIsParticle && tracksForKalmanFit[0].charge() > 0) {
           lamDauMasses.push_back(protonMass);
           lamDauMasses.push_back(piMass);
       } else if (!lamIsParticle && tracksForKalmanFit[0].charge() < 0) {
           lamDauMasses.push_back(protonMass);
           lamDauMasses.push_back(piMass);
       } else {
           lamDauMasses.push_back(piMass);
           lamDauMasses.push_back(protonMass);
       }

      // Loop over tracks and vertex good charged track pairs
      for(unsigned int trdx = 0; trdx < theTrackRefs.size(); trdx++) {

         if ( theTrackRefs[trdx].isNull() ) continue;

         math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
         double dzvtx = theTrackRefs[trdx]->dz(bestvtx);
         double dxyvtx = theTrackRefs[trdx]->dxy(bestvtx);
         if(fabs(dzvtx)>50 || fabs(dxyvtx)>50) continue;

         double dzerror = sqrt(theTrackRefs[trdx]->dzError()*theTrackRefs[trdx]->dzError()+zVtxError*zVtxError);
         double dxyerror = sqrt(theTrackRefs[trdx]->d0Error()*theTrackRefs[trdx]->d0Error()+xVtxError*yVtxError);
         double dauLongImpactSig = dzvtx/dzerror;
         double dauTransImpactSig = dxyvtx/dxyerror;
         if( fabs(dauTransImpactSig) < batDauTransImpactSigCut || fabs(dauLongImpactSig) < batDauLongImpactSigCut ) continue;

         bool match = false;
         // check if the pion is used in any v0 candidate in the collection and flag it
         for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
           if (theTrackRefs[trdx]->charge() == theDaughterTracks[j]->charge() &&
               theTrackRefs[trdx]->momentum() == theDaughterTracks[j]->momentum() ) match = true;
           if (match) break;
         }
         if (match) continue; // Track is already used in making the V0

         // check if pion is in *any* good V0
         // Placeholder
       
         TransientTrack pionTT(theDaughterTracks[0], &(*bFieldHandle) );
         TransientTrack protonTT(theDaughterTracks[1], &(*bFieldHandle) );
         TransientTrack batPionTT(theTrackRefs[trdx], &(*bFieldHandle) );

         if (!batPionTT.isValid()) continue;
         if (!pionTT.isValid()) continue;
         if (!protonTT.isValid()) continue;

         //Creating a KinematicParticleFactory
         KinematicParticleFactoryFromTransientTrack pFactory;

         //initial chi2 and ndf before kinematic fits.
         float chi = 0.;
         float ndf = 0.;
         vector<RefCountedKinematicParticle> lambdaParticles;
         lambdaParticles.push_back(pFactory.particle(pionTT,piMass,chi,ndf,piMass_sigma));
         lambdaParticles.push_back(pFactory.particle(protonTT,protonMass,chi,ndf,protonMass_sigma));

         KinematicParticleVertexFitter fitter;
         RefCountedKinematicTree lambdaVertexFitTree;
         lambdaVertexFitTree = fitter.fit(lambdaParticles);
         if (!lambdaVertexFitTree->isValid()) continue;
  
         lambdaVertexFitTree->movePointerToTheTop();

         RefCountedKinematicParticle lambda_vFit = lambdaVertexFitTree->currentParticle();
         RefCountedKinematicVertex lambda_vFit_vertex = lambdaVertexFitTree->currentDecayVertex();

         lambdaVertexFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle lambdaPi1 = lambdaVertexFitTree->currentParticle();
         lambdaVertexFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle lambdaPi2 = lambdaVertexFitTree->currentParticle();

//         KinematicParameters lambdaPiKP = lambdaPi1->currentState().kinematicParameters();
//         KinematicParameters lambdaPKP  = lambdaPi2->currentState().kinematicParameters();

         // now apply Lambda mass constraint
         KinematicParticleFitter csFitterLambda;
         KinematicConstraint * lam_c = new MassKinematicConstraint(lambdaMass,lambdaMass_sigma);
         // add mass constraint to the lambda fit to do a constrained fit:  

         lambdaVertexFitTree->movePointerToTheTop();
         lambdaVertexFitTree = csFitterLambda.fit(lam_c,lambdaVertexFitTree);
         if (!lambdaVertexFitTree->isValid()) continue;
         lambdaVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle lambda_vFit_withMC = lambdaVertexFitTree->currentParticle();

         if (!lambda_vFit_withMC->currentState().isValid()) continue; 

         if(doXis || doLambdaCToLamPis)
         {
           if(doXis && !doLambdaCToLamPis)
           {
             if(lamIsParticle && theTrackRefs[trdx]->charge()>0) continue;
             if(!lamIsParticle && theTrackRefs[trdx]->charge()<0) continue;
           }

           if(!doXis && doLambdaCToLamPis)
           {
             if(lamIsParticle && theTrackRefs[trdx]->charge()<0) continue;
             if(!lamIsParticle && theTrackRefs[trdx]->charge()>0) continue;
           }
 
           vector<RefCountedKinematicParticle> xiFitParticles;

           xiFitParticles.push_back(pFactory.particle(batPionTT,piMass,chi,ndf,piMass_sigma));
           xiFitParticles.push_back(lambda_vFit_withMC);

           //fit Xi 
           RefCountedKinematicTree xiFitTree = fitter.fit(xiFitParticles);
           if (!xiFitTree->isValid()) continue;

           xiFitTree->movePointerToTheTop();
           RefCountedKinematicParticle xiCand = xiFitTree->currentParticle();
           if (!xiCand->currentState().isValid()) continue;

           RefCountedKinematicVertex xiDecayVertex = xiFitTree->currentDecayVertex();
           if (!xiDecayVertex->vertexIsValid()) continue;

/*
           if ( xiDecayVertex->chiSquared()<0 || xiDecayVertex->chiSquared()>1000 ) continue;

           float xiC2Prob =
              ChiSquaredProbability((double)(xiDecayVertex->chiSquared()),(double)(xiDecayVertex->degreesOfFreedom()));
           if (xiC2Prob < 0.0001) continue;
*/
           if ( xiCand->currentState().mass() > 3 ) continue;

           // get children from final Xi fit
           xiFitTree->movePointerToTheFirstChild();
           RefCountedKinematicParticle batPionCand = xiFitTree->currentParticle();
           xiFitTree->movePointerToTheNextChild();
           RefCountedKinematicParticle lambdaCandMC = xiFitTree->currentParticle();
         
           if(!batPionCand->currentState().isValid() || !lambdaCandMC->currentState().isValid()) continue;

           // get batchlor pion and V0 parameters from Xi fit
           KinematicParameters batPionKP = batPionCand->currentState().kinematicParameters();
           KinematicParameters lambdaCandKP = lambdaCandMC->currentState().kinematicParameters();

           GlobalVector xiTotalP = GlobalVector (xiCand->currentState().globalMomentum().x(),
                                                 xiCand->currentState().globalMomentum().y(),
                                                 xiCand->currentState().globalMomentum().z());

           GlobalVector batPionTotalP = GlobalVector(batPionKP.momentum().x(),batPionKP.momentum().y(),batPionKP.momentum().z());
           GlobalVector lambdaTotalP = GlobalVector(lambdaCandKP.momentum().x(),lambdaCandKP.momentum().y(),lambdaCandKP.momentum().z());

           double batPionTotalE = sqrt( batPionTotalP.mag2() + piMassSquared );
           double lambdaTotalE = sqrt( lambdaTotalP.mag2() + lambdaMass*lambdaMass );
           double xiTotalE = batPionTotalE + lambdaTotalE;

           const Particle::LorentzVector xiP4(xiTotalP.x(), xiTotalP.y(), xiTotalP.z(), xiTotalE);

           Particle::Point xiVtx((*xiDecayVertex).position().x(), (*xiDecayVertex).position().y(), (*xiDecayVertex).position().z());
           std::vector<double> xiVtxEVec;
           xiVtxEVec.push_back( lambda_vFit_vertex->error().cxx() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().cyx() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().cyy() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().czx() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().czy() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().czz() );
           SMatrixSym3D xiVtxCovMatrix(xiVtxEVec.begin(), xiVtxEVec.end());
           const Vertex::CovarianceMatrix xiVtxCov(xiVtxCovMatrix);
           double xiVtxChi2(xiDecayVertex->chiSquared());
           double xiVtxNdof(xiDecayVertex->degreesOfFreedom());
           double xiNormalizedChi2 = xiVtxChi2/xiVtxNdof;

           double xiRVtxMag = 99999.0;
           double xiLVtxMag = 99999.0;
           double xiSigmaRvtxMag = 999.0;
           double xiSigmaLvtxMag = 999.0;
           double xiV0Angle = -100.0;

           GlobalVector xiV0LineOfFlight = GlobalVector (xiVtx.x() - xVtx,
                                                         xiVtx.y() - yVtx,
                                                         xiVtx.z() - zVtx);

           SMatrixSym3D xiTotalCov;
           if(isVtxPV) xiTotalCov = xiVtxCovMatrix + vtxPrimary->covariance();
           else xiTotalCov = xiVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

           SVector3 distanceVector3D(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), xiV0LineOfFlight.z());
           SVector3 distanceVector2D(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), 0.0);

           xiV0Angle = angle(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), xiV0LineOfFlight.z(),
                             xiTotalP.x(), xiTotalP.y(), xiTotalP.z());
           xiLVtxMag = xiV0LineOfFlight.mag();
           xiRVtxMag = xiV0LineOfFlight.perp();
           xiSigmaLvtxMag = sqrt(ROOT::Math::Similarity(xiTotalCov, distanceVector3D)) / xiLVtxMag;
           xiSigmaRvtxMag = sqrt(ROOT::Math::Similarity(xiTotalCov, distanceVector2D)) / xiRVtxMag;

           if( xiNormalizedChi2 > xiChi2Cut ||
               xiRVtxMag < xiRVtxCut ||
               xiRVtxMag / xiSigmaRvtxMag < xiRVtxSigCut ||
               xiLVtxMag < xiLVtxCut ||
               xiLVtxMag / xiSigmaLvtxMag < xiLVtxSigCut ||
               cos(xiV0Angle) < xiCollinCut
           ) continue;

           // Cuts finished, Create the VertexCompositeCandidate object for Xi that will be stored in the Event
           RecoChargedCandidate PionCand(theTrackRefs[trdx]->charge(), Particle::LorentzVector(batPionTotalP.x(),
                                                   batPionTotalP.y(), batPionTotalP.z(),
                                                   batPionTotalE), xiVtx);
           PionCand.setTrack(theTrackRefs[trdx]);

           AddFourMomenta addp4;
           // Store the daughter Candidates in the VertexCompositeCandidates 
           // if they pass mass cuts
           if(doXis)
           {
             if(lamIsParticle && PionCand.charge()>0) continue;
             if(!lamIsParticle && PionCand.charge()<0) continue;

             VertexCompositeCandidate* theXi = 0;
             theXi = new VertexCompositeCandidate(0, xiP4, xiVtx, xiVtxCov, xiVtxChi2, xiVtxNdof);
             theXi->addDaughter(theLambda);
             theXi->addDaughter(PionCand);
             if(lamIsParticle) { theXi->setPdgId(3312); theXi->setCharge(-1); }
             else { theXi->setPdgId(-3312); theXi->setCharge(1); }
             addp4.set( *theXi );
             if( theXi->mass() < xiMass + xiMassCut &&
                 theXi->mass() > xiMass - xiMassCut ) theXis.push_back( *theXi );
             if(theXi) delete theXi;
             theXi = 0;
           }
           if(doLambdaCToLamPis)
           {
             if(lamIsParticle && PionCand.charge()<0) continue;
             if(!lamIsParticle && PionCand.charge()>0) continue;

             VertexCompositeCandidate* theLambdaCToLamPi = 0;
             theLambdaCToLamPi = new VertexCompositeCandidate(0, xiP4, xiVtx, xiVtxCov, xiVtxChi2, xiVtxNdof);
             theLambdaCToLamPi->addDaughter(theLambda);
             theLambdaCToLamPi->addDaughter(PionCand);
             if(lamIsParticle) { theLambdaCToLamPi->setPdgId(4122); theLambdaCToLamPi->setCharge(1); }
             else { theLambdaCToLamPi->setPdgId(-4122); theLambdaCToLamPi->setCharge(-1); } 
             addp4.set( *theLambdaCToLamPi );
             if( theLambdaCToLamPi->mass() < lambdaCMass + lambdaCMassCut &&
                 theLambdaCToLamPi->mass() > lambdaCMass - lambdaCMassCut ) theLambdaCToLamPis.push_back( *theLambdaCToLamPi );
             if(theLambdaCToLamPi) delete theLambdaCToLamPi;
             theLambdaCToLamPi = 0;
           }
        }

         if(doOmegas)
         {
           if(lamIsParticle && theTrackRefs[trdx]->charge()>0) continue;
           if(!lamIsParticle && theTrackRefs[trdx]->charge()<0) continue;

           vector<RefCountedKinematicParticle> omegaFitParticles;

           omegaFitParticles.push_back(pFactory.particle(batPionTT,kaonMass,chi,ndf,kaonMass_sigma));
           omegaFitParticles.push_back(lambda_vFit_withMC);

           //fit Omega 
           RefCountedKinematicTree omegaFitTree = fitter.fit(omegaFitParticles);
           if (!omegaFitTree->isValid()) continue;

           omegaFitTree->movePointerToTheTop();
           RefCountedKinematicParticle omegaCand = omegaFitTree->currentParticle();
           if (!omegaCand->currentState().isValid()) continue;

           RefCountedKinematicVertex omegaDecayVertex = omegaFitTree->currentDecayVertex();
           if (!omegaDecayVertex->vertexIsValid()) continue;

           if ( omegaDecayVertex->chiSquared()<0 || omegaDecayVertex->chiSquared()>1000 ) continue;

           float omegaC2Prob =
              ChiSquaredProbability((double)(omegaDecayVertex->chiSquared()),(double)(omegaDecayVertex->degreesOfFreedom()));
           if (omegaC2Prob < 0.0001) continue;

           if ( omegaCand->currentState().mass() > 3 ) continue;

           // get children from final Omega fit
           omegaFitTree->movePointerToTheFirstChild();
           RefCountedKinematicParticle batPionCand = omegaFitTree->currentParticle();
           omegaFitTree->movePointerToTheNextChild();
           RefCountedKinematicParticle lambdaCandMC = omegaFitTree->currentParticle();
         
           if(!batPionCand->currentState().isValid() || !lambdaCandMC->currentState().isValid()) continue;

           // get batchlor pion and V0 parameters from Xi fit
           KinematicParameters batPionKP = batPionCand->currentState().kinematicParameters();
           KinematicParameters lambdaCandKP = lambdaCandMC->currentState().kinematicParameters();

           GlobalVector omegaTotalP = GlobalVector (omegaCand->currentState().globalMomentum().x(),
                                                    omegaCand->currentState().globalMomentum().y(),
                                                    omegaCand->currentState().globalMomentum().z());

           GlobalVector batPionTotalP = GlobalVector(batPionKP.momentum().x(),batPionKP.momentum().y(),batPionKP.momentum().z());
           GlobalVector lambdaTotalP = GlobalVector(lambdaCandKP.momentum().x(),lambdaCandKP.momentum().y(),lambdaCandKP.momentum().z());


           double batPionTotalE = sqrt( batPionTotalP.mag2() + kaonMassSquared );
           double lambdaTotalE = sqrt( lambdaTotalP.mag2() + lambdaMass*lambdaMass );
           double omegaTotalE = batPionTotalE + lambdaTotalE;

           const Particle::LorentzVector omegaP4(omegaTotalP.x(), omegaTotalP.y(), omegaTotalP.z(), omegaTotalE);

           Particle::Point omegaVtx((*omegaDecayVertex).position().x(), (*omegaDecayVertex).position().y(), (*omegaDecayVertex).position().z());
           std::vector<double> omegaVtxEVec;
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().cxx() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().cyx() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().cyy() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().czx() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().czy() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().czz() );
           SMatrixSym3D omegaVtxCovMatrix(omegaVtxEVec.begin(), omegaVtxEVec.end());
           const Vertex::CovarianceMatrix omegaVtxCov(omegaVtxCovMatrix);
           double omegaVtxChi2(omegaDecayVertex->chiSquared());
           double omegaVtxNdof(omegaDecayVertex->degreesOfFreedom());
           double omegaNormalizedChi2 = omegaVtxChi2/omegaVtxNdof;

           double omegaRVtxMag = 99999.0;
           double omegaLVtxMag = 99999.0;
           double omegaSigmaRvtxMag = 999.0;
           double omegaSigmaLvtxMag = 999.0;
           double omegaV0Angle = -100.0;

           GlobalVector omegaV0LineOfFlight = GlobalVector (omegaVtx.x() - xVtx,
                                                            omegaVtx.y() - yVtx,
                                                            omegaVtx.z() - zVtx);

           SMatrixSym3D omegaTotalCov;
           if(isVtxPV) omegaTotalCov = omegaVtxCovMatrix + vtxPrimary->covariance();
           else omegaTotalCov = omegaVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

           SVector3 distanceVector3D(omegaV0LineOfFlight.x(), omegaV0LineOfFlight.y(), omegaV0LineOfFlight.z());
           SVector3 distanceVector2D(omegaV0LineOfFlight.x(), omegaV0LineOfFlight.y(), 0.0);

           omegaV0Angle = angle(omegaV0LineOfFlight.x(), omegaV0LineOfFlight.y(), omegaV0LineOfFlight.z(),
                             omegaTotalP.x(), omegaTotalP.y(), omegaTotalP.z());
           omegaLVtxMag = omegaV0LineOfFlight.mag();
           omegaRVtxMag = omegaV0LineOfFlight.perp();
           omegaSigmaLvtxMag = sqrt(ROOT::Math::Similarity(omegaTotalCov, distanceVector3D)) / omegaLVtxMag;
           omegaSigmaRvtxMag = sqrt(ROOT::Math::Similarity(omegaTotalCov, distanceVector2D)) / omegaRVtxMag;

           if( omegaNormalizedChi2 > xiChi2Cut ||
               omegaRVtxMag < xiRVtxCut ||
               omegaRVtxMag / omegaSigmaRvtxMag < xiRVtxSigCut ||
               omegaLVtxMag < xiLVtxCut ||
               omegaLVtxMag / omegaSigmaLvtxMag < xiLVtxSigCut ||
               cos(omegaV0Angle) < xiCollinCut
           ) continue;

           // Cuts finished, Create the VertexCompositeCandidate object for Omega that will be stored in the Event
           VertexCompositeCandidate* theOmega = 0;
           theOmega = new VertexCompositeCandidate(0, omegaP4, omegaVtx, omegaVtxCov, omegaVtxChi2, omegaVtxNdof);

           RecoChargedCandidate PionCand(theTrackRefs[trdx]->charge(), Particle::LorentzVector(batPionTotalP.x(),
                                                   batPionTotalP.y(), batPionTotalP.z(),
                                                   batPionTotalE), omegaVtx);
           PionCand.setTrack(theTrackRefs[trdx]);

           AddFourMomenta addp4;
           // Store the daughter Candidates in the VertexCompositeCandidates 
           // if they pass mass cuts
           theOmega->addDaughter(theLambda);
           theOmega->addDaughter(PionCand);
           if(lamIsParticle) { theOmega->setPdgId(3334);  theOmega->setCharge(-2); }
           else { theOmega->setPdgId(-3334);  theOmega->setCharge(2); } 
           addp4.set( *theOmega );
           if( theOmega->mass() < omegaMass + omegaMassCut &&
               theOmega->mass() > omegaMass - omegaMassCut ) theOmegas.push_back( *theOmega );
           if(theOmega) delete theOmega;
           theOmega = 0;
	 }
      }
    }
  }
}

// Get methods
const reco::VertexCompositeCandidateCollection& V0Fitter::getKshorts() const {
  return theKshorts;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getPhis() const {
  return thePhis;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getLambdas() const {
  return theLambdas;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getXis() const {
  return theXis;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getOmegas() const {
  return theOmegas;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getD0() const {
  return theD0s;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getDSToKsK() const {
  return theDSToKsKs;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getDSToPhiPi() const {
  return theDSToPhiPis;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getDPM() const {
  return theDPMs;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getLambdaCToLamPi() const {
  return theLambdaCToLamPis;
}

const reco::VertexCompositeCandidateCollection& V0Fitter::getLambdaCToKsP() const {
  return theLambdaCToKsPs;
}

void V0Fitter::resetAll() {
  theKshorts.clear();
  thePhis.clear();
  theLambdas.clear();
  theXis.clear();
  theOmegas.clear();
  theD0s.clear();
  theDSToKsKs.clear();
  theDSToPhiPis.clear();
  theDPMs.clear();
  theLambdaCToLamPis.clear();
  theLambdaCToKsPs.clear();
}

// Experimental
double V0Fitter::findV0MassError(const GlobalPoint &vtxPos, std::vector<reco::TransientTrack> dauTracks) { 
  return -1.;
}

/*
double V0Fitter::findV0MassError(const GlobalPoint &vtxPos, std::vector<reco::TransientTrack> dauTracks) {
  // Returns -99999. if trajectory states fail at vertex position

  // Load positive track trajectory at vertex into vector, then negative track
  std::vector<TrajectoryStateClosestToPoint> sortedTrajStatesAtVtx;
  for( unsigned int ndx = 0; ndx < dauTracks.size(); ndx++ ) {
    if( dauTracks[ndx].trajectoryStateClosestToPoint(vtxPos).isValid() ) {
      std::cout << "From TSCP: " 
		<< dauTracks[ndx].trajectoryStateClosestToPoint(vtxPos).perigeeParameters().transverseCurvature()
		<< "; From Track: " << dauTracks[ndx].track().qoverp() << std::endl;
    }
    if( sortedTrajStatesAtVtx.size() == 0 ) {
      if( dauTracks[ndx].charge() > 0 ) {
	sortedTrajStatesAtVtx.push_back( dauTracks[ndx].trajectoryStateClosestToPoint(vtxPos) );
      }
      else {
	sortedTrajStatesAtVtx.push_back( dauTracks[ndx].trajectoryStateClosestToPoint(vtxPos) );
      }
    }
  }
  std::vector<PerigeeTrajectoryParameters> param;
  std::vector<PerigeeTrajectoryError> paramError;
  std::vector<GlobalVector> momenta;

  for( unsigned int ndx2 = 0; ndx2 < sortedTrajStatesAtVtx.size(); ndx2++ ) {
    if( sortedTrajStatesAtVtx[ndx2].isValid() ) {
      param.push_back( sortedTrajStatesAtVtx[ndx2].perigeeParameters() );
      paramError.push_back( sortedTrajStatesAtVtx[ndx2].perigeeError() );
      momenta.push_back( sortedTrajStatesAtVtx[ndx2].momentum() );
    }
    else return -99999.;
  }
  return 0;
}
*/



