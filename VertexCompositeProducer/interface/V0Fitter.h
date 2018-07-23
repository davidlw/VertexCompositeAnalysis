// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      V0Fitter
// 
/**\class V0Fitter V0Fitter.h VertexCompositeAnalysis/VertexCompositeProducer/interface/V0Fitter.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#ifndef VertexCompositeAnalysis__V0_FITTER_H
#define VertexCompositeAnalysis__V0_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Math/interface/angle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
//#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <fstream>

class V0Fitter {
 public:
  V0Fitter(const edm::ParameterSet& theParams, edm::ConsumesCollector && iC);
//	   const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::ConsumesCollector && iC);
  ~V0Fitter();

  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  // Switching to L. Lista's reco::Candidate infrastructure for V0 storage
  const reco::VertexCompositeCandidateCollection& getKshorts() const;
  const reco::VertexCompositeCandidateCollection& getPhis() const;
  const reco::VertexCompositeCandidateCollection& getLambdas() const;
  const reco::VertexCompositeCandidateCollection& getXis() const;
  const reco::VertexCompositeCandidateCollection& getOmegas() const;
  const reco::VertexCompositeCandidateCollection& getD0() const;
  const reco::VertexCompositeCandidateCollection& getDSToKsK() const;
  const reco::VertexCompositeCandidateCollection& getDSToPhiPi() const;
  const reco::VertexCompositeCandidateCollection& getDPM() const;
  const reco::VertexCompositeCandidateCollection& getLambdaCToLamPi() const;
  const reco::VertexCompositeCandidateCollection& getLambdaCToKsP() const;
  void resetAll();

 private:
  // STL vector of VertexCompositeCandidate that will be filled with VertexCompositeCandidates by fitAll()
  reco::VertexCompositeCandidateCollection theKshorts;
  reco::VertexCompositeCandidateCollection thePhis;
  reco::VertexCompositeCandidateCollection theLambdas;
  reco::VertexCompositeCandidateCollection theXis;
  reco::VertexCompositeCandidateCollection theOmegas;
  reco::VertexCompositeCandidateCollection theD0s;
  reco::VertexCompositeCandidateCollection theDSToKsKs;
  reco::VertexCompositeCandidateCollection theDSToPhiPis;
  reco::VertexCompositeCandidateCollection theDPMs;
  reco::VertexCompositeCandidateCollection theLambdaCToLamPis;
  reco::VertexCompositeCandidateCollection theLambdaCToKsPs;

  // Tracker geometry for discerning hit positions
  const TrackerGeometry* trackerGeom;

  const MagneticField* magField;

  edm::InputTag recoAlg;
  edm::InputTag vtxAlg;
  edm::EDGetTokenT<reco::TrackCollection> token_tracks;
  edm::EDGetTokenT<reco::VertexCollection> token_vertices;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;

  bool useRefTrax;
  bool storeRefTrax;
  bool doKshorts;
  bool doPhis;
  bool doLambdas;
  bool doXis;
  bool doOmegas;
  bool doD0s;
  bool doDSToKsKs;
  bool doDSToPhiPis;
  bool doDPMs;
  bool doLambdaCToLamPis;
  bool doLambdaCToKsPs;

  bool doVertexFit;

  /*bool doPostFitCuts;
    bool doTkQualCuts;*/

  // Cuts
  double tkDCACut;
  double tkChi2Cut;
  int    tkNhitsCut;
  double tkPtCut;
  double chi2Cut;
  double rVtxCut;
  double rVtxSigCut;
  double lVtxCut;
  double lVtxSigCut;
  double collinCut;
  double xiChi2Cut;
  double xiRVtxCut;
  double xiRVtxSigCut;
  double xiLVtxCut;
  double xiLVtxSigCut;
  double xiCollinCut;
  double kShortMassCut;
  double phiMassCut;
  double lambdaMassCut;
  double d0MassCut;
  double dsMassCut;
  double dpmMassCut;
  double lambdaCMassCut;
  double xiMassCut;
  double omegaMassCut;
  double dauTransImpactSigCut;
  double dauLongImpactSigCut;
  double batDauTransImpactSigCut;
  double batDauLongImpactSigCut;
  double mPiPiCutMin;
  double mPiPiCutMax;
  double mKKCutMin;
  double mKKCutMax;
  double innerHitPosCut;

  std::vector<reco::TrackBase::TrackQuality> qualities;

  edm::InputTag vtxFitter;

  // Helper method that does the actual fitting using the KalmanVertexFitter
  double findV0MassError(const GlobalPoint &vtxPos, std::vector<reco::TransientTrack> dauTracks);

  // Applies cuts to the VertexCompositeCandidates after they are fitted/created.
  //void applyPostFitCuts();

  // Stuff for debug file output.
  std::ofstream mPiPiMassOut;

  inline void initFileOutput() {
    mPiPiMassOut.open("mPiPi.txt", std::ios::app);
  }
  inline void cleanupFileOutput() {
    mPiPiMassOut.close();
  }
};

#endif
