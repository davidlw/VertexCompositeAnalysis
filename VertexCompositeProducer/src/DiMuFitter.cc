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
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

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

const double muonMass = 0.10565837;

// Constructor and (empty) destructor
DiMuFitter::DiMuFitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) :
  muonSelection(theParameters.getParameter<std::string>("muonSelection")),
  candidateSelection(theParameters.getParameter<std::string>("candidateSelection"))
{

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_muons = iC.consumes<pat::MuonCollection>(theParameters.getParameter<edm::InputTag>("muonRecoAlgorithm"));

  // Second, initialize post-fit cuts
  mllCutMin = theParameters.getParameter<double>("mllCutMin");
  mllCutMax = theParameters.getParameter<double>("mllCutMax");
  tkdXYCut = theParameters.getParameter<double>("tkdXYCut");
  tkdZCut = theParameters.getParameter<double>("tkdZCut");
  tkDCACut = theParameters.getParameter<double>("tkDCACut");
  chi2Cut = theParameters.getParameter<double>("vtxChi2Cut");
  rVtxCut = theParameters.getParameter<double>("rVtxCut");
  rVtxSigCut = theParameters.getParameter<double>("vtxSignificance2DCut");
  lVtxCut = theParameters.getParameter<double>("lVtxCut");
  lVtxSigCut = theParameters.getParameter<double>("vtxSignificance3DCut");
  collinCut = theParameters.getParameter<double>("collinearityCut");
  dauTransImpactSigCut = theParameters.getParameter<double>("dauTransImpactSigCut");
  dauLongImpactSigCut = theParameters.getParameter<double>("dauLongImpactSigCut");
  VtxChiProbCut = theParameters.getParameter<double>("VtxChiProbCut");
  alphaCut = theParameters.getParameter<double>("alphaCut");
  isWrongSign = theParameters.getParameter<bool>("isWrongSign");
}

DiMuFitter::~DiMuFitter() {
}

// Method containing the algorithm for vertex reconstruction
void DiMuFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  const double& dauMass = muonMass;
  const double& dauMassSquared = dauMass * dauMass;
  float dauMass_sigma = dauMass * 1.e-6;

  // Handles for tracks, B-field, and tracker geometry
  edm::Handle<reco::VertexCollection> theVertexHandle;
  edm::Handle<reco::PFCandidateCollection> thePfCandHandle;
  edm::Handle<pat::MuonCollection> theMuonHandle;

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  edm::ESHandle<MagneticField> bFieldHandle;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_muons, theMuonHandle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);

  if( !theMuonHandle->size() ) return;

  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  magField = bFieldHandle.product();

  // Find the primary vertex
  bool isVtxPV = 0;
  double xVtx = theBeamSpotHandle->position().x();
  double yVtx = theBeamSpotHandle->position().y();
  double zVtx = 0.0;
  double xVtxError = theBeamSpotHandle->BeamWidthX();
  double yVtxError = theBeamSpotHandle->BeamWidthY();
  double zVtxError = 0.0;
  const reco::Vertex& vtxPrimary = (theVertexHandle->size()>0 ? *(theVertexHandle->begin()) : reco::Vertex());
  if (!vtxPrimary.isFake() && vtxPrimary.tracksSize()>=2) {
    isVtxPV = 1;
    xVtx = vtxPrimary.x();
    yVtx = vtxPrimary.y();
    zVtx = vtxPrimary.z();
    xVtxError = vtxPrimary.xError();
    yVtxError = vtxPrimary.yError();
    zVtxError = vtxPrimary.zError();
  }
  const math::XYZPoint bestvtx(xVtx, yVtx, zVtx);

  // Select the muon candidates
  pat::MuonCollection muonColl;
  for (const auto& muon : *theMuonHandle) {

    if (!muonSelection(muon)) continue;

    const reco::TrackRef& trackRef = muon.track();
    if(trackRef.isNull()) continue;

    const double& dzvtx  = trackRef->dz(bestvtx);
    const double& dxyvtx = trackRef->dxy(bestvtx);
    if( fabs(dzvtx)>tkdZCut || fabs(dxyvtx)>tkdXYCut ) continue;

    const double& dzerror  = std::sqrt(trackRef->dzError()*trackRef->dzError()+zVtxError*zVtxError);
    const double& dxyerror = std::sqrt(trackRef->d0Error()*trackRef->d0Error()+xVtxError*yVtxError);
    const double& dauLongImpactSig  = dzvtx/dzerror;
    const double& dauTransImpactSig = dxyvtx/dxyerror;
    if( fabs(dauTransImpactSig) < dauTransImpactSigCut || fabs(dauLongImpactSig) < dauLongImpactSigCut ) continue;

    muonColl.push_back(muon);
  }

  // Loop over the selected muon candidates
  for (uint ic = 0; ic < muonColl.size(); ic++) {

    const pat::Muon& cand1 = muonColl[ic];
    const reco::TrackRef& trackRef1 = cand1.track();
    const reco::TransientTrack tmpTk1( *trackRef1, magField );

    for( uint fc = ic+1; fc < muonColl.size(); fc++ ) {

       const pat::Muon& cand2 = muonColl[fc];

       const double& totalE = std::sqrt( cand1.p() + dauMassSquared ) + std::sqrt( cand2.p() + dauMassSquared );
       const double& totalESq = totalE*totalE;
       const double& totalPSq = ( cand1.momentum() + cand2.momentum() ).mag2();
       const double& mass = std::sqrt( totalESq - totalPSq);
       if( (mass > mllCutMax || mass < mllCutMin) && (mass > mllCutMax || mass < mllCutMin)) continue;

       const reco::TrackRef& trackRef2 = cand2.track();
       const reco::TransientTrack tmpTk2( *trackRef2, magField );

       if(isWrongSign==false && (cand1.charge()*cand2.charge()) > 0.0) continue;
       if(isWrongSign==true  && (cand1.charge()*cand2.charge()) < 0.0) continue;
       const auto& posCand = (cand1.charge()>0 ? cand1 : cand2);
       const auto& negCand = (cand1.charge()<0 ? cand1 : cand2);
       const auto& posTransTkPtr = (cand1.charge()>0 ? &tmpTk1 : &tmpTk2);
       const auto& negTransTkPtr = (cand1.charge()<0 ? &tmpTk1 : &tmpTk2);

       // Trajectory states to calculate DCA for the 2 tracks
       if( !posTransTkPtr->impactPointTSCP().isValid() || !negTransTkPtr->impactPointTSCP().isValid() ) continue;
       const FreeTrajectoryState& posState = posTransTkPtr->impactPointTSCP().theState();
       const FreeTrajectoryState& negState = negTransTkPtr->impactPointTSCP().theState();

       // Measure distance between tracks at their closest approach
       ClosestApproachInRPhi cApp;
       cApp.calculate(posState, negState);
       if( !cApp.status() ) continue;

       const double& dca = fabs( cApp.distance() );
       if (dca < 0. || dca > tkDCACut) continue;
       const GlobalPoint& cxPt = cApp.crossingPoint();

       // Get trajectory states for the tracks at POCA for later cuts
       const TrajectoryStateClosestToPoint& posTSCP = posTransTkPtr->trajectoryStateClosestToPoint( cxPt );
       const TrajectoryStateClosestToPoint& negTSCP = negTransTkPtr->trajectoryStateClosestToPoint( cxPt );
       if( !posTSCP.isValid() || !negTSCP.isValid() ) continue;

       // Create the vertex fitter object and vertex the tracks
       float chi = 0.0 , ndf = 0.0;
       std::vector<RefCountedKinematicParticle> DiMuParticles;
       KinematicParticleFactoryFromTransientTrack pFactory;
       DiMuParticles.push_back(pFactory.particle(*posTransTkPtr, dauMass, chi, ndf, dauMass_sigma));
       DiMuParticles.push_back(pFactory.particle(*negTransTkPtr, dauMass, chi, ndf, dauMass_sigma));

       KinematicParticleVertexFitter DiMuFitter;
       const RefCountedKinematicTree& DiMuVertex = DiMuFitter.fit(DiMuParticles);
       if( !DiMuVertex->isValid() ) continue;

       DiMuVertex->movePointerToTheTop();
       const RefCountedKinematicParticle& DiMuCand = DiMuVertex->currentParticle();
       if (!DiMuCand->currentState().isValid()) continue;

       const RefCountedKinematicVertex& DiMuDecayVertex = DiMuVertex->currentDecayVertex();
       if (!DiMuDecayVertex->vertexIsValid()) continue;

       const double& DiMuVtxChi2 = DiMuDecayVertex->chiSquared();
       const double& DiMuVtxNdof = DiMuDecayVertex->degreesOfFreedom();
       const double& DiMuC2Prob  = TMath::Prob(DiMuVtxChi2, DiMuVtxNdof);
       if (DiMuC2Prob < VtxChiProbCut) continue;

       const double& DiMuNormalizedChi2 = DiMuVtxChi2/DiMuVtxNdof;
       if(DiMuNormalizedChi2 > chi2Cut) continue;

       const auto candMom = posCand.momentum() + negCand.momentum();
       const GlobalVector& DiMuTotalP    = GlobalVector(candMom.x(), candMom.y(), candMom.z());
       const GlobalVector& posCandTotalP = GlobalVector(posCand.momentum().x(), posCand.momentum().y(),posCand.momentum().z());
       const GlobalVector& negCandTotalP = GlobalVector(negCand.momentum().x(), negCand.momentum().y(),negCand.momentum().z());

       const float posCandTotalE = std::sqrt( posCandTotalP.mag2() + dauMassSquared );
       const float negCandTotalE = std::sqrt( negCandTotalP.mag2() + dauMassSquared );
       const float DiMuTotalE = posCandTotalE + negCandTotalE;

       const reco::Particle::LorentzVector DiMuP4(DiMuTotalP.x(), DiMuTotalP.y(), DiMuTotalP.z(), DiMuTotalE);

       const reco::Particle::Point DiMuVtx(DiMuDecayVertex->position().x(), DiMuDecayVertex->position().y(), DiMuDecayVertex->position().z());
       std::vector<double> DiMuVtxEVec;
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().cxx() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().cyx() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().cyy() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().czx() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().czy() );
       DiMuVtxEVec.push_back( DiMuDecayVertex->error().czz() );
       const SMatrixSym3D DiMuVtxCovMatrix(DiMuVtxEVec.begin(), DiMuVtxEVec.end());
       const reco::Vertex::CovarianceMatrix DiMuVtxCov(DiMuVtxCovMatrix);

       const GlobalVector& DiMuLineOfFlight = GlobalVector(DiMuVtx.x() - xVtx, DiMuVtx.y() - yVtx, DiMuVtx.z() - zVtx);

       const double& lVtxMag = DiMuLineOfFlight.mag();
       const double& rVtxMag = DiMuLineOfFlight.perp();
       if (rVtxMag < rVtxCut || lVtxMag < lVtxCut) continue;

       const double& DiMuAngle = angle(DiMuLineOfFlight.x(), DiMuLineOfFlight.y(), DiMuLineOfFlight.z(), DiMuTotalP.x(), DiMuTotalP.y(), DiMuTotalP.z());
       if (std::cos(DiMuAngle) < collinCut) continue;

       const SVector3 distanceVector3D(DiMuLineOfFlight.x(), DiMuLineOfFlight.y(), DiMuLineOfFlight.z());
       const SVector3 distanceVector2D(DiMuLineOfFlight.x(), DiMuLineOfFlight.y(), 0.0);

       const SMatrixSym3D& DiMuTotalCov = DiMuVtxCovMatrix + (isVtxPV ? vtxPrimary.covariance() : theBeamSpotHandle->rotatedCovariance3D());

       const double& sigmaLvtxMag = std::sqrt(ROOT::Math::Similarity(DiMuTotalCov, distanceVector3D)) / lVtxMag;
       const double& sigmaRvtxMag = std::sqrt(ROOT::Math::Similarity(DiMuTotalCov, distanceVector2D)) / rVtxMag;
       if ((rVtxMag / sigmaRvtxMag < rVtxSigCut) || (lVtxMag / sigmaLvtxMag < lVtxSigCut)) continue;

       const TVector3 svpvVec(DiMuVtx.x() - xVtx, DiMuVtx.y() - yVtx, DiMuVtx.z() - zVtx);
       const TVector3 dVec(DiMuTotalP.x(), DiMuTotalP.y(), DiMuTotalP.z());
       const double& alpha = svpvVec.Angle(dVec);
       if (alpha > alphaCut) continue;

       const auto& diMu = reco::CompositeCandidate(0, DiMuP4, DiMuVtx);
       pat::CompositeCandidate theDiMu(diMu);
       theDiMu.addUserFloat("vertexNdof", DiMuVtxNdof);
       theDiMu.addUserFloat("vertexChi2", DiMuVtxChi2);
       theDiMu.addUserData<reco::Vertex::CovarianceMatrix>("vertexCovariance", DiMuVtxCov);

       theDiMu.addDaughter(posCand, "posDau");
       theDiMu.addDaughter(negCand, "negDau");
       AddFourMomenta addP4;
       addP4.set(theDiMu);

       // Store the dimuon candidate
       if (candidateSelection(theDiMu)) theDiMus.push_back(theDiMu);
     }
   }
}
// Get methods

const pat::CompositeCandidateCollection& DiMuFitter::getDiMu() const {
  return theDiMus;
}

void DiMuFitter::resetAll() {
    theDiMus.clear();
}
