// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      BFitter
// 
/**\class BFitter BFitter.h VertexCompositeAnalysis/VertexCompositeProducer/interface/BFitter.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li
//
//

#ifndef VertexCompositeAnalysis__B_FITTER_H
#define VertexCompositeAnalysis__B_FITTER_H

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
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/EgammaObjects/interface/GBRForest.h"

#include <string>
#include <fstream>
#include <typeinfo>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>

class BFitter {
 public:
  BFitter(const edm::ParameterSet& theParams, edm::ConsumesCollector && iC);
  ~BFitter();

  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  const reco::VertexCompositeCandidateCollection& getB() const;
//  const std::vector<float>& getMVAVals() const; 

  void resetAll();

 private:
  reco::VertexCompositeCandidateCollection theBs;

  // Tracker geometry for discerning hit positions
  const TrackerGeometry* trackerGeom;

  const MagneticField* magField;

  edm::InputTag recoAlg;
  edm::InputTag vtxAlg;
  edm::InputTag d0Alg;
  edm::EDGetTokenT<reco::TrackCollection> token_tracks;
  edm::EDGetTokenT<reco::VertexCollection> token_vertices;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> token_d0s;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > token_dedx;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;

  // Cuts
  double mPiDCutMin;
  double mPiDCutMax;
  double batTkDCACut;
  double batTkChi2Cut;
  int    batTkNhitsCut;
  double batTkPtErrCut;
  double batTkPtCut;
  double batTkEtaCut;
  double batTkPtSumCut;
  double batTkEtaDiffCut;
  double bVtxChi2Cut;
  double bRVtxCut;
  double bRVtxSigCut;
  double bLVtxCut;
  double bLVtxSigCut;
  double bCollinCut2D;
  double bCollinCut3D;
  double bMassCut;
  double batDauTransImpactSigCut;
  double batDauLongImpactSigCut;
  double bVtxChiProbCut;
  double bPtCut;
  double bAlphaCut;
  double bAlpha2DCut;
  bool   isWrongSignB;

  std::vector<reco::TrackBase::TrackQuality> qualities;

  //setup mva selector
/*
  bool useAnyMVA_;
  std::vector<bool> useMVA_;
  std::vector<double> min_MVA_;
  std::string mvaType_;
  std::string forestLabel_;
  GBRForest * forest_;
  bool useForestFromDB_;

  std::vector<float> mvaVals_;
*/
//  auto_ptr<edm::ValueMap<float> >mvaValValueMap;
//  MVACollection mvas; 

//  std::string dbFileName_;

};

#endif
