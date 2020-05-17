// -*- C++ -*-
//
// Package:    ParticleProducer
// Class:      ParticleFitter
// 
/**\class ParticleFitter ParticleFitter.h VertexCompositeAnalysis/VertexCompositeProducer/interface/ParticleFitter.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andre Stahl
//
//

#ifndef VertexCompositeAnalysis__Particle_FITTER_H
#define VertexCompositeAnalysis__Particle_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Math/interface/angle.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/VirtualKinematicParticleFactory.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "combinations.hpp"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>
#include <TVector3.h>


struct ParticleComparator {
  template<class T>
  const bool isEqual(const T& x, const T& y) const
  {
    return std::fabs(x-y) <= 1.E-6;
  }
  template<class T>
  const bool isLess(const T& x, const T& y) const
  {
    return !isEqual(x, y) && x < y;
  }
  const bool isParticleEqual(const pat::GenericParticle& x, const pat::GenericParticle& y) const
  {
    return (isEqual(x.pt(), y.pt()) && isEqual(x.eta(), y.eta()) && isEqual(x.phi(), y.phi()) && isEqual(x.charge(), y.charge()));
  }
  const bool isParticleLess(const pat::GenericParticle& x, const pat::GenericParticle& y) const
  {
    return (isLess(x.pt(), y.pt()) ||
            (isEqual(x.pt(), y.pt()) && isLess(x.eta(), y.eta())) ||
            (isEqual(x.pt(), y.pt()) && isEqual(x.eta(), y.eta()) && isLess(x.phi(), y.phi())) ||
            (isEqual(x.pt(), y.pt()) && isEqual(x.eta(), y.eta()) && isEqual(x.phi(), y.phi()) && isLess(x.charge(), y.charge())));
  }
  inline bool operator()(const pat::GenericParticle& x, const pat::GenericParticle& y) const
  {
    return isParticleLess(x, y);
  }
};


struct ParticleMassComparator : ParticleComparator {
  inline bool operator()(const pat::GenericParticle& x, const pat::GenericParticle& y) const
  {
    return (isParticleLess(x, y) ||
            (isEqual(x.pt(), y.pt()) && isEqual(x.eta(), y.eta()) && isEqual(x.phi(), y.phi()) && isEqual(x.charge(), y.charge()) && isLess(x.mass(), y.mass())));
  }
};


const std::map<uint, double> MASS_ = {
  {11, 0.000511}, {13, 0.10565837}, {15, 1.77686}, // leptons
  {211, 0.13957018}, {310, 0.497614}, {321, 0.493677}, {333, 1.019445}, // light and strange mesons
  {411, 1.86962}, {421, 1.86484}, {431, 1.96847}, {511, 5.27929}, // charmed and bottom mesons
  {2212, 0.938272013}, {3122, 1.115683}, {3312, 1.32171}, {3334, 1.67245}, {4122, 2.28646} // baryons
};
const std::map<uint, float> WIDTH_ = {{211, 3.5E-7f}, {321, 1.6E-5f}, {2212, 1.6E-5f}};


typedef std::set<pat::GenericParticle, ParticleComparator> ParticleSet;
typedef std::set<pat::GenericParticle, ParticleMassComparator> ParticleMassSet;
typedef ROOT::Math::SVector<double, 3> SVector3;


class ParticleDaughter {
 public:
  ParticleDaughter();
  ParticleDaughter(const edm::ParameterSet& pSet, const edm::ParameterSet& config, edm::ConsumesCollector&& iC);
  ~ParticleDaughter();

  const int& pdgId() const { return pdgId_; }
  const int& charge() const { return charge_; }
  const double& mass() const { return mass_; }
  const float& width() const { return width_; }
  const pat::GenericParticleCollection& particles() const { return particles_; }
  const bool useSource() const { return !token_source_.isUninitialized(); }
  
  template <class T>
  void addParticles(const edm::Event& event, const edm::EDGetTokenT<std::vector<T> >& token, const reco::Vertex& vertex, const bool embedInfo=true);
  void addParticles(const edm::Event& event);
  void fillInfo(const edm::ParameterSet& pSet, const edm::ParameterSet& config, edm::ConsumesCollector& iC);
  void clear();

 private:
  template <class T>
  void addData(pat::GenericParticle& c, const edm::Ref<std::vector<T> >& p, const bool& embedInfo);
  void addData(pat::GenericParticle& c, const reco::TrackRef& p, const bool& embedInfo);
  void addData(pat::GenericParticle& c, const reco::PFCandidateRef& p, const bool& embedInfo);
  void addData(pat::GenericParticle& c, const pat::MuonRef& p, const bool& embedInfo);
  void addData(pat::GenericParticle& c, const pat::ElectronRef& p, const bool& embedInfo);
  void setMVA (pat::GenericParticle& c, const size_t& i, const edm::Handle<std::vector<float> >& m);
  void setDeDx(pat::GenericParticle& c, const edm::Handle<edm::ValueMap<reco::DeDxData> >& m);

  int pdgId_;
  int charge_;
  double mass_;
  float width_;
  std::string selection_;
  std::string finalSelection_;
  pat::GenericParticleCollection particles_;

  edm::EDGetTokenT<pat::GenericParticleCollection> token_source_;
  edm::EDGetTokenT<std::vector<float> >            token_mva_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > token_dedx_;
};


class ParticleFitter {
 public:
  ParticleFitter(const edm::ParameterSet& pSet, edm::ConsumesCollector&& iC);
  ~ParticleFitter();

  typedef std::map<double, std::vector<double> > DoubleMap;

  const pat::GenericParticleCollection& particles() const { return candidates_; }
  const bool hasNoDaughters() const { return daughters_.empty(); }

  void setVertex(const edm::Event& iEvent);
  void addParticles(ParticleDaughter& d, const edm::Event& iEvent);
  void fillDaughters(const edm::Event& iEvent);
  bool isUniqueDaughter(ParticleSet& set, const pat::GenericParticle& dau);
  void makeCandidates();
  void swapDaughters(DoubleMap& swapDauColls, const pat::GenericParticle& cand);
  void setBestMass(pat::GenericParticle& cand, const DoubleMap& swapDauColls);
  void addSwapCandidates(pat::GenericParticleCollection& swapCandColl, const pat::GenericParticle& cand, const DoubleMap& swapDauColls);
  bool fitCandidate(pat::GenericParticle& cand);
  void addExtraInfo(pat::GenericParticle& cand);
  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void clear();

 private:
  int pdgId_;
  bool doSwap_;
  double mass_, width_;
  reco::Vertex vertex_;
  reco::Vertex beamSpot2D_;
  std::vector<ParticleDaughter> daughters_;
  pat::GenericParticleCollection candidates_;

  edm::ESHandle<MagneticField> bFieldHandle_;

  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot_;
  edm::EDGetTokenT<reco::VertexCollection> token_vertices_;
  edm::EDGetTokenT<pat::ElectronCollection> token_electrons_;
  edm::EDGetTokenT<pat::MuonCollection> token_muons_;
  edm::EDGetTokenT<pat::TauCollection> token_taus_;
  edm::EDGetTokenT<pat::PhotonCollection> token_photons_;
  edm::EDGetTokenT<reco::TrackCollection> token_tracks_;
  edm::EDGetTokenT<reco::PFCandidateCollection> token_pfParticles_;
  edm::EDGetTokenT<pat::JetCollection> token_jets_;

  StringCutObjectSelector<pat::GenericParticle, true> preSelection_;
  StringCutObjectSelector<pat::GenericParticle, true> pocaSelection_;
  StringCutObjectSelector<pat::GenericParticle, true> postSelection_;
  StringCutObjectSelector<pat::GenericParticle, true> finalSelection_;
};


#endif
