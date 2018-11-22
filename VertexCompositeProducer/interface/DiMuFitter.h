// -*- C++ -*-
//
// Package:    DiMuProducer
// Class:      DiMuFitter
// 
/**\class DiMuFitter DiMuFitter.h VertexCompositeAnalysis/VertexCompositeProducer/interface/DiMuFitter.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li
//
//

#ifndef VertexCompositeAnalysis__DiMu_FITTER_H
#define VertexCompositeAnalysis__DiMu_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/angle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
//include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <fstream>

class DiMuFitter {
 public:
  DiMuFitter(const edm::ParameterSet& theParams, edm::ConsumesCollector && iC);
  ~DiMuFitter();

  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  const reco::VertexCompositeCandidateCollection& getDiMu() const;
  void resetAll();

 private:
  // STL vector of VertexCompositeCandidate that will be filled with VertexCompositeCandidates by fitAll()
  reco::VertexCompositeCandidateCollection theDiMus;

  // Tracker geometry for discerning hit positions
  const TrackerGeometry* trackerGeom;

  const MagneticField* magField;

  edm::InputTag recoAlg;
  edm::InputTag vtxAlg;
  edm::EDGetTokenT<reco::TrackCollection> token_tracks;
  edm::EDGetTokenT<reco::VertexCollection> token_vertices;
  edm::EDGetTokenT<reco::PFCandidateCollection> token_pfcands;
  edm::EDGetTokenT<reco::MuonCollection> token_muons;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;

  // Cuts
  double mllCutMin;
  double mllCutMax;
  double tkDCACut;
  double tkChi2Cut;
  int    tkNhitsCut;
  double tkPtCut;
  double tkEtaCut;
  double chi2Cut;
  double rVtxCut;
  double rVtxSigCut;
  double lVtxCut;
  double lVtxSigCut;
  double collinCut;
  double DiMuMassCut;
  double dauTransImpactSigCut;
  double dauLongImpactSigCut;
  double VtxChiProbCut;
  double dPtCut;
  double alphaCut;
  bool   isEE;
  bool   isMuMu;
  std::string muonId;
  bool   isMuonId;
  bool   isPFMuon;
  bool   isGlobalMuon;
  bool   isWrongSign;

  std::vector<reco::TrackBase::TrackQuality> qualities;
};

#endif
