// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      D0Fitter
// 
/**\class D0Fitter D0Fitter.h VertexCompositeAnalysis/VertexCompositeProducer/interface/D0Fitter.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li
//
//

#ifndef VertexCompositeAnalysis__D0_FITTER_H
#define VertexCompositeAnalysis__D0_FITTER_H

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
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
//#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

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

class D0Fitter {
 public:
  D0Fitter(const edm::ParameterSet& theParams, edm::ConsumesCollector && iC);
  ~D0Fitter();

  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  // Switching to L. Lista's reco::Candidate infrastructure for D0 storage
  const reco::VertexCompositeCandidateCollection& getD0() const;
  const std::vector<float>& getMVAVals() const; 

//  auto_ptr<edm::ValueMap<float> > getMVAMap() const;
  void resetAll();

 private:
  // STL vector of VertexCompositeCandidate that will be filled with VertexCompositeCandidates by fitAll()
  reco::VertexCompositeCandidateCollection theD0s;

  // Tracker geometry for discerning hit positions
  const TrackerGeometry* trackerGeom;

  const MagneticField* magField;

  edm::InputTag recoAlg;
  edm::InputTag vtxAlg;
  edm::EDGetTokenT<reco::TrackCollection> token_tracks;
  edm::EDGetTokenT<reco::VertexCollection> token_vertices;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > token_dedx;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
  ///cesar -begin
  edm::InputTag mvaTrackRecoSrcLabel_;
  edm::EDGetTokenT<std::vector<float>> mvaTrackRecoSrc_;
  ///cesar -end

  // Cuts
  double mPiKCutMin;
  double mPiKCutMax;
  double tkDCACut;
  double tkChi2Cut;
  int    tkNhitsCut;
  double tkPtErrCut;
  double tkPtCut;
  double tkEtaCut;
  double tkPtSumCut;
  double tkEtaDiffCut;
  double chi2Cut;
  double rVtxCut;
  double rVtxSigCut;
  double lVtxCut;
  double lVtxSigCut;
  double collinCut2D;
  double collinCut3D;
  double d0MassCut;
  double dauTransImpactSigCut;
  double dauLongImpactSigCut;
  double VtxChiProbCut;
  double dPtCut;
  double alphaCut;
  double alpha2DCut;
  bool   isWrongSign;

  std::vector<reco::TrackBase::TrackQuality> qualities;

  //setup mva selector
  bool useAnyMVA_;
  std::vector<bool> useMVA_;
  std::vector<double> min_MVA_;
  std::string mvaType_;
  std::string forestLabel_;
  GBRForest * forest_;
  bool useForestFromDB_;

  std::vector<float> mvaVals_;

//  auto_ptr<edm::ValueMap<float> >mvaValValueMap;
//  MVACollection mvas; 

  std::string dbFileName_;

};

#endif
