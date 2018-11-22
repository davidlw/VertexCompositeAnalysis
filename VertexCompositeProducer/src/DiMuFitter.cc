// -*- C++ -*-
//
// Package:    DiMuProducer
// Class:      DiMuFitter
// 
/**\class DiMuFitter DiMuFitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/DiMuFitter.cc
//
*/

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiMuFitter.h"
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/PatExamples/plugins/MuonAnalyzer.h"

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

const float electronMass = 0.000510999;
const float muonMass = 0.10565837;
const float DiMuMass = 3.096916;
float DiMuMass_sigma = DiMuMass*1.e-6;

// Constructor and (empty) destructor
DiMuFitter::DiMuFitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) {

  using std::string;

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_pfcands = iC.consumes<reco::PFCandidateCollection>(theParameters.getParameter<edm::InputTag>("pfCandAlgorithm"));
  token_muons = iC.consumes<reco::MuonCollection>(theParameters.getParameter<edm::InputTag>("muonRecoAlgorithm"));

  // Second, initialize post-fit cuts
  mllCutMin = theParameters.getParameter<double>(string("mllCutMin"));
  mllCutMax = theParameters.getParameter<double>(string("mllCutMax"));
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
  DiMuMassCut = theParameters.getParameter<double>(string("DiMuMassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  VtxChiProbCut = theParameters.getParameter<double>(string("VtxChiProbCut"));
  dPtCut = theParameters.getParameter<double>(string("dPtCut"));
  alphaCut = theParameters.getParameter<double>(string("alphaCut"));
  muonId = theParameters.getParameter<std::string>(string("muonId"));
  isMuonId = theParameters.getParameter<bool>(string("isMuonId")); 
  isPFMuon = theParameters.getParameter<bool>(string("isPFMuon")); 
  isGlobalMuon = theParameters.getParameter<bool>(string("isGlobalMuon"));
  isWrongSign = theParameters.getParameter<bool>(string("isWrongSign"));

  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }
}

DiMuFitter::~DiMuFitter() {
}

// Method containing the algorithm for vertex reconstruction
void DiMuFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  float dauMass = muonMass;
  float dauMass_sigma = dauMass * 1.e-6;
  float dauMassSquared = dauMass * dauMass;

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<reco::PFCandidateCollection> thePfCandHandle;
  Handle<reco::MuonCollection> theMuonHandle;

  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks, theTrackHandle); 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_pfcands, thePfCandHandle);
  iEvent.getByToken(token_muons, theMuonHandle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);  

  if( !theTrackHandle->size() ) return;
  if( !theMuonHandle->size() ) return;

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
  math::XYZPoint bestvtx(xVtx,yVtx,zVtx);

   for( unsigned ic = 0; ic < theMuonHandle->size(); ic++ ) {

     const reco::Muon& cand1 = (*theMuonHandle)[ic];
     if(isMuonId && !muon::isGoodMuon(cand1, muon::selectionTypeFromString(muonId))) continue;  //DataFormats/MuonReco/interface/MuonSelectors.h, TMOneStationTight = 12
     if(isPFMuon && !cand1.isPFMuon()) continue; 
     if(isGlobalMuon && !cand1.isGlobalMuon()) continue;

//   recoMu.numberOfMatchedStations() > 1
/*
     const reco::PFCandidate& cand1 = (*thePfCandHandle)[ic];
     int type1 = cand1.particleId();
     if(isEE && type1 != reco::PFCandidate::e) continue;
     if(isMuMu && type1 != reco::PFCandidate::mu) continue;

     reco::TrackRef trackRef1 = cand1.trackRef();
*/

     reco::TrackRef trackRef1 = cand1.track();
     if(trackRef1.isNull()) continue;

     bool quality_ok = true;
     if (qualities.size()!=0) {
       quality_ok = false;
       for (unsigned int ndx_ = 0; ndx_ < qualities.size(); ndx_++) {
         if (trackRef1->quality(qualities[ndx_])){
           quality_ok = true;
           break;
         }
       }
     }
     if( !quality_ok ) continue;

     TransientTrack tmpTk1( *trackRef1, magField );
     if( trackRef1->normalizedChi2() < tkChi2Cut &&
         trackRef1->numberOfValidHits() >= tkNhitsCut &&
         trackRef1->hitPattern().pixelLayersWithMeasurement() > 0 &&
         trackRef1->pt() > tkPtCut && fabs(trackRef1->eta()) < tkEtaCut) {

       double dzvtx = trackRef1->dz(bestvtx);
       double dxyvtx = trackRef1->dxy(bestvtx);
       double dzerror = sqrt(trackRef1->dzError()*trackRef1->dzError()+zVtxError*zVtxError);
       double dxyerror = sqrt(trackRef1->d0Error()*trackRef1->d0Error()+xVtxError*yVtxError);

       double dauLongImpactSig = dzvtx/dzerror;
       double dauTransImpactSig = dxyvtx/dxyerror;

       if( fabs(dzvtx)>20. || fabs(dxyvtx)>0.3 ) continue;
       if( fabs(dauTransImpactSig) < dauTransImpactSigCut || fabs(dauLongImpactSig) < dauLongImpactSigCut ) continue;
     }

     for( unsigned fc = ic+1; fc < theMuonHandle->size(); fc++ ) {

       const reco::Muon& cand2 = (*theMuonHandle)[fc];
       if(isMuonId && !muon::isGoodMuon(cand2, muon::selectionTypeFromString(muonId))) continue;  
       if(isPFMuon && !cand2.isPFMuon()) continue;
       if(isGlobalMuon && !cand2.isGlobalMuon()) continue;

/*
       const reco::PFCandidate& cand2 = (*thePfCandHandle)[fc];
       int type2 = cand2.particleId();
       if(isEE && type2 != reco::PFCandidate::e) continue;
       if(isMuMu && type2 != reco::PFCandidate::mu) continue;
*/
       double totalE = sqrt( cand1.p() + dauMassSquared ) +
                       sqrt( cand2.p() + dauMassSquared );
       double totalESq = totalE*totalE;

       double totalPSq =
         ( cand1.momentum() + cand2.momentum() ).mag2();

       double mass = sqrt( totalESq - totalPSq);

       if( (mass > mllCutMax || mass < mllCutMin) && (mass > mllCutMax || mass < mllCutMin)) continue;

       reco::TrackRef trackRef2 = cand2.track();
       if(trackRef2.isNull()) continue;
       
       if (qualities.size()!=0) {
         quality_ok = false;
         for (unsigned int ndx_ = 0; ndx_ < qualities.size(); ndx_++) {
           if (trackRef2->quality(qualities[ndx_])){
             quality_ok = true;
             break;
           }
         }
       }
       if( !quality_ok ) continue;

       TransientTrack tmpTk2( *trackRef2, magField );
       if( trackRef2->normalizedChi2() < tkChi2Cut &&
           trackRef2->numberOfValidHits() >= tkNhitsCut &&
           trackRef2->hitPattern().pixelLayersWithMeasurement() > 0 &&
           trackRef2->pt() > tkPtCut && fabs(trackRef2->eta()) < tkEtaCut) {

         double dzvtx = trackRef2->dz(bestvtx);
         double dxyvtx = trackRef2->dxy(bestvtx);
         double dzerror = sqrt(trackRef2->dzError()*trackRef2->dzError()+zVtxError*zVtxError);
         double dxyerror = sqrt(trackRef2->d0Error()*trackRef2->d0Error()+xVtxError*yVtxError);

         double dauLongImpactSig = dzvtx/dzerror;
         double dauTransImpactSig = dxyvtx/dxyerror;

         if( fabs(dzvtx)>20. || fabs(dxyvtx)>0.3 ) continue;
         if( fabs(dauTransImpactSig) < dauTransImpactSigCut || fabs(dauLongImpactSig) < dauLongImpactSigCut ) continue;
       }

//       reco::PFCandidate posCand;
//       reco::PFCandidate negCand;
       reco::Muon posCand;
       reco::Muon negCand;
       TrackRef positiveTrackRef;
       TrackRef negativeTrackRef;
       TransientTrack* posTransTkPtr = 0;
       TransientTrack* negTransTkPtr = 0;

       if(!isWrongSign && trackRef1->charge() < 0. &&
          trackRef2->charge() > 0.) {
         posCand = cand2;
         negCand = cand1;
         negativeTrackRef = trackRef1;
         positiveTrackRef = trackRef2;
         negTransTkPtr = &tmpTk1;
         posTransTkPtr = &tmpTk2;
       }
       else if(!isWrongSign && trackRef1->charge() > 0. &&
               trackRef2->charge() < 0.) {
         posCand = cand1;
         negCand = cand2;
         negativeTrackRef = trackRef2;
         positiveTrackRef = trackRef1;
         negTransTkPtr = &tmpTk2;
         posTransTkPtr = &tmpTk1;
       }
       else if(isWrongSign &&  trackRef1->charge()* trackRef2->charge() > 0.0) {
         posCand = cand1;
         negCand = cand2;
         negativeTrackRef = trackRef2;
         positiveTrackRef = trackRef1;
         negTransTkPtr = &tmpTk2;
         posTransTkPtr = &tmpTk1;
       }
       else continue;

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

       // Get trajectory states for the tracks at POCA for later cuts
       TrajectoryStateClosestToPoint posTSCP =
         posTransTkPtr->trajectoryStateClosestToPoint( cxPt );
       TrajectoryStateClosestToPoint negTSCP =
         negTransTkPtr->trajectoryStateClosestToPoint( cxPt );

       if( !posTSCP.isValid() || !negTSCP.isValid() ) continue;

       // Create the vertex fitter object and vertex the tracks
       float posCandTotalE=0.0;
       float negCandTotalE=0.0;
       float DiMuTotalE=0.0;

       KinematicParticleFactoryFromTransientTrack pFactory;

       float chi = 0.0;
       float ndf = 0.0;

       vector<RefCountedKinematicParticle> DiMuParticles;
       DiMuParticles.push_back(pFactory.particle(*posTransTkPtr,dauMass,chi,ndf,dauMass_sigma));
       DiMuParticles.push_back(pFactory.particle(*negTransTkPtr,dauMass,chi,ndf,dauMass_sigma));

       KinematicParticleVertexFitter DiMuFitter;
       RefCountedKinematicTree DiMuVertex;
       DiMuVertex = DiMuFitter.fit(DiMuParticles);

       if( !DiMuVertex->isValid() ) continue;

       DiMuVertex->movePointerToTheTop();
       RefCountedKinematicParticle DiMuCand = DiMuVertex->currentParticle();

       if (!DiMuCand->currentState().isValid()) continue;

       RefCountedKinematicVertex DiMuDecayVertex = DiMuVertex->currentDecayVertex();
       if (!DiMuDecayVertex->vertexIsValid()) continue;

       float DiMuC2Prob = TMath::Prob(DiMuDecayVertex->chiSquared(),DiMuDecayVertex->degreesOfFreedom());
       if (DiMuC2Prob < VtxChiProbCut) continue;

       GlobalVector DiMuTotalP = GlobalVector ((posCand.momentum() + negCand.momentum()).x(),
                                             (posCand.momentum() + negCand.momentum()).y(),
                                             (posCand.momentum() + negCand.momentum()).z());

       GlobalVector posCandTotalP = GlobalVector(posCand.momentum().x(),posCand.momentum().y(),posCand.momentum().z());
       GlobalVector negCandTotalP = GlobalVector(negCand.momentum().x(),negCand.momentum().y(),negCand.momentum().z());

       posCandTotalE = sqrt( posCandTotalP.mag2() + dauMassSquared );
       negCandTotalE = sqrt( negCandTotalP.mag2() + dauMassSquared );
       DiMuTotalE = posCandTotalE + negCandTotalE;

       const Particle::LorentzVector DiMuP4(DiMuTotalP.x(), DiMuTotalP.y(), DiMuTotalP.z(), DiMuTotalE);

       Particle::Point DiMuVtx((*DiMuDecayVertex).position().x(), (*DiMuDecayVertex).position().y(), (*DiMuDecayVertex).position().z());
       std::vector<double> DiMuVtxEVec;
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().cxx() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().cyx() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().cyy() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().czx() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().czy() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().czz() );
       SMatrixSym3D DiMuVtxCovMatrix(DiMuVtxEVec.begin(), DiMuVtxEVec.end());
       const Vertex::CovarianceMatrix DiMuVtxCov(DiMuVtxCovMatrix);
       double DiMuVtxChi2(DiMuDecayVertex->chiSquared());
       double DiMuVtxNdof(DiMuDecayVertex->degreesOfFreedom());
       double DiMuNormalizedChi2 = DiMuVtxChi2/DiMuVtxNdof;

       double rVtxMag = 99999.0; 
       double lVtxMag = 99999.0;
       double sigmaRvtxMag = 999.0;
       double sigmaLvtxMag = 999.0;
       double DiMuAngle = -100.0;

       GlobalVector DiMuLineOfFlight = GlobalVector (DiMuVtx.x() - xVtx,
                                                   DiMuVtx.y() - yVtx,
                                                   DiMuVtx.z() - zVtx);

       SMatrixSym3D DiMuTotalCov;
       if(isVtxPV) DiMuTotalCov = DiMuVtxCovMatrix + vtxPrimary->covariance();
       else DiMuTotalCov = DiMuVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

       SVector3 distanceVector3D(DiMuLineOfFlight.x(), DiMuLineOfFlight.y(), DiMuLineOfFlight.z());
       SVector3 distanceVector2D(DiMuLineOfFlight.x(), DiMuLineOfFlight.y(), 0.0);

       DiMuAngle = angle(DiMuLineOfFlight.x(), DiMuLineOfFlight.y(), DiMuLineOfFlight.z(),
                        DiMuTotalP.x(), DiMuTotalP.y(), DiMuTotalP.z());
       lVtxMag = DiMuLineOfFlight.mag();
       rVtxMag = DiMuLineOfFlight.perp();
       sigmaLvtxMag = sqrt(ROOT::Math::Similarity(DiMuTotalCov, distanceVector3D)) / lVtxMag;
       sigmaRvtxMag = sqrt(ROOT::Math::Similarity(DiMuTotalCov, distanceVector2D)) / rVtxMag;

       TVector3 svpvVec;
       svpvVec.SetXYZ(DiMuVtx.x() - xVtx,
                      DiMuVtx.y() - yVtx,
                      DiMuVtx.z() - zVtx);
       TVector3 dVec;
       dVec.SetXYZ(DiMuTotalP.x(), DiMuTotalP.y(), DiMuTotalP.z());
       double alpha = svpvVec.Angle(dVec);
          
       if( DiMuNormalizedChi2 > chi2Cut ||
           rVtxMag < rVtxCut ||
           rVtxMag / sigmaRvtxMag < rVtxSigCut ||
           lVtxMag < lVtxCut ||
           lVtxMag / sigmaLvtxMag < lVtxSigCut ||
           cos(DiMuAngle) < collinCut || alpha > alphaCut
       ) continue;

       VertexCompositeCandidate* theDiMu = 0;
       theDiMu = new VertexCompositeCandidate(0, DiMuP4, DiMuVtx, DiMuVtxCov, DiMuVtxChi2, DiMuVtxNdof);

       RecoChargedCandidate
         thePosCand(1, Particle::LorentzVector(posCandTotalP.x(),
                                               posCandTotalP.y(), posCandTotalP.z(),
                                               posCandTotalE), DiMuVtx);
       thePosCand.setTrack(positiveTrackRef);

       RecoChargedCandidate
         theNegCand(-1, Particle::LorentzVector(negCandTotalP.x(),
                                                negCandTotalP.y(), negCandTotalP.z(),
                                                negCandTotalE), DiMuVtx);
       theNegCand.setTrack(negativeTrackRef);

       if(isWrongSign)
       {
         thePosCand.setCharge(trackRef1->charge());
         theNegCand.setCharge(trackRef2->charge());
       }

       AddFourMomenta addp4;
       theDiMu->addDaughter(thePosCand);
       theDiMu->addDaughter(theNegCand);
       theDiMu->setPdgId(443);
       addp4.set( *theDiMu );
       if( theDiMu->mass() < DiMuMass + DiMuMassCut &&
           theDiMu->mass() > DiMuMass - DiMuMassCut &&
	   theDiMu->pt() > dPtCut ) {
         theDiMus.push_back( *theDiMu );
       }

       if(theDiMu) delete theDiMu;
     }
   }
}
// Get methods

const reco::VertexCompositeCandidateCollection& DiMuFitter::getDiMu() const {
  return theDiMus;
}

void DiMuFitter::resetAll() {
    theDiMus.clear();
}
