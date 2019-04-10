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

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
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

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Math/interface/angle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <fstream>

class DiMuFitter {
 public:
  DiMuFitter(const edm::ParameterSet& theParams, edm::ConsumesCollector && iC);
  ~DiMuFitter();
  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  const pat::CompositeCandidateCollection& getDiMu() const;
  void resetAll();

 private:
  // STL vector of VertexCompositeCandidate that will be filled with VertexCompositeCandidates by fitAll()
  pat::CompositeCandidateCollection theDiMus;

  // Tracker geometry for discerning hit positions
  const MagneticField* magField;

  edm::EDGetTokenT<reco::VertexCollection> token_vertices;
  edm::EDGetTokenT<pat::MuonCollection> token_muons;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;

  // Cuts
  double mllCutMin;
  double mllCutMax;
  double tkdXYCut;
  double tkdZCut;
  double tkDCACut;
  double chi2Cut;
  double rVtxCut;
  double rVtxSigCut;
  double lVtxCut;
  double lVtxSigCut;
  double collinCut;
  double dauTransImpactSigCut;
  double dauLongImpactSigCut;
  double VtxChiProbCut;
  double alphaCut;
  bool   isWrongSign;
  StringCutObjectSelector<pat::Muon> muonSelection;
  StringCutObjectSelector<reco::Candidate, true> candidateSelection;
};

#endif
