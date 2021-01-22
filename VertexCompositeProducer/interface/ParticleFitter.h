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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Math/interface/angle.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/GaussianSumVertexFit/interface/GsfVertexFitter.h"
#include "RecoVertex/GaussianSumVertexFit/interface/AdaptiveGsfVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/VirtualKinematicParticleFactory.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicPerigeeConversions.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/FinalTreeBuilder.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "HepPDT/ParticleID.hh"

#include "combinations.hpp"


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
            (isParticleEqual(x, y) && isLess(x.mass(), y.mass())));
  }
};
struct ParticleTreeComparator : ParticleComparator {
  inline bool operator()(const pat::GenericParticle& x, const pat::GenericParticle& y) const
  {
    return (isLess(x.mass(), y.mass()) ||
            (isEqual(x.mass(), y.mass()) && isLess(x.charge(), y.charge())) ||
            (isEqual(x.mass(), y.mass()) && isEqual(x.charge(), y.charge()) && isParticleLess(x, y)));
  }
};


const std::map<uint, double> MASS_ = {
  {0, -1.}, // default
  {11, 0.000511}, {13, 0.10565837}, {15, 1.77686}, // leptons
  {23, 91.188}, {24, 80.38}, // bosons
  {211, 0.13957018}, {310, 0.497614}, {321, 0.493677}, {333, 1.019445}, // light and strange mesons
  {411, 1.86962}, {421, 1.86484}, {431, 1.96847}, {511, 5.27929}, // charmed and bottom mesons
  {443, 3.096900}, {100443, 3.68609}, {553, 9.4603}, {100553, 10.0233}, {200553, 10.3552}, // quarkonia
  {2212, 0.938272013}, {3122, 1.115683}, {3312, 1.32171}, {3334, 1.67245}, {4122, 2.28646} // baryons
};
const std::map<uint, float> WIDTH_ = {{211, 3.5E-7f}, {321, 1.6E-5f}, {2212, 1.6E-5f}};


typedef std::set<pat::GenericParticle, ParticleComparator> ParticleSet;
typedef std::set<pat::GenericParticle, ParticleMassComparator> ParticleMassSet;


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
  void init(const edm::EventSetup& iSetup);
  void clear();
  template<class T>
  void clear(std::vector<T>& v) { std::vector<T>().swap(v); };

  pat::GenericParticleCollection particles_;

 private:
  template <class T>
  void addInfo(pat::GenericParticle& c, const T& p);
  void addInfo(pat::GenericParticle& c, const reco::Conversion& p);
  template <class T>
  void addData(pat::GenericParticle& c, const edm::Ref<std::vector<T> >& p, const bool& embedInfo);
  void addData(pat::GenericParticle& c, const reco::TrackRef& p, const bool& embedInfo);
  void addData(pat::GenericParticle& c, const reco::PFCandidateRef& p, const bool& embedInfo);
  void addData(pat::GenericParticle& c, const pat::MuonRef& p, const bool& embedInfo);
  void addData(pat::GenericParticle& c, const pat::ElectronRef& p, const bool& embedInfo);
  void setMVA (pat::GenericParticle& c, const size_t& i, const edm::Handle<std::vector<float> >& m);
  void setDeDx(pat::GenericParticle& c, const std::map<std::string, edm::Handle<edm::ValueMap<reco::DeDxData> > >& m);
  void addMuonL1Info(pat::GenericParticle& c, const edm::Handle<pat::TriggerObjectStandAloneMatch>& m);

  int pdgId_;
  int charge_;
  double mass_;
  float width_;
  std::string selection_;
  std::string finalSelection_;

  edm::ParameterSet conf_;
  PropagateToMuon* propToMuon_;

  edm::EDGetTokenT<pat::GenericParticleCollection> token_source_;
  edm::EDGetTokenT<std::vector<float> >            token_mva_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneMatch> token_muonL1Info_;

  std::map<std::string, edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > > tokens_dedx_;
};


class ParticleFitter {
 public:
  ParticleFitter(const edm::ParameterSet& pSet, edm::ConsumesCollector&& iC);
  ~ParticleFitter();

  typedef std::map<double, std::vector<double> > DoubleMap;
  typedef edm::RefProd<pat::GenericParticleCollection> GenericParticleRefProd;
  typedef std::vector<math::XYZTLorentzVector> LorentzVectorColl;
  typedef std::vector<reco::TransientTrack> TransTrackColl;
  typedef std::vector<RefCountedKinematicParticle> KinParColl;
  typedef std::tuple<float, float, float, size_t, bool> VertexTuple;
  typedef std::tuple<float, float, float, float, signed char> ParticleTuple;
  typedef std::tuple<KinParColl, TransTrackColl, std::map<ParticleTuple, size_t> > ParticleInfo;

  const reco::VertexCollection& vertices() const { return vertices_; };
  const pat::GenericParticleCollection& particles() const { return candidates_; };
  const pat::GenericParticleCollection& daughters() const { return particles_; };
  const bool hasNoDaughters() const { return daughters_.empty(); };

  reco::VertexRef getVertexRef(const reco::Vertex& vertex);
  math::XYZTLorentzVector getP4(const GlobalVector& p, const double& m);
  pat::GenericParticleRef addParticle(const pat::GenericParticle& particle);
  void matchPrimaryVertex(pat::GenericParticle& cand, const TransTrackColl& tracks={}, FreeTrajectoryState fts={}, const double& thr=1.E-9);
  RefCountedKinematicTree fitVertex(const ParticleInfo& parInfo, const int& fitAlgo, GlobalPoint decP, const reco::Vertex& priVtx={});
  RefCountedKinematicTree fitVertex(const ParticleInfo& parInfo, const int& fitAlgo);

  void setVtxProd(const reco::VertexRefProd& prod) { vtxProd_ = prod; };
  void setDauProd(const GenericParticleRefProd& prod) { dauProd_ = prod; };
  void setVertex(const edm::Event& iEvent);
  void getNTracks(const edm::Event& iEvent);
  void addParticles(ParticleDaughter& d, const edm::Event& iEvent);
  void fillDaughters(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool isUniqueDaughter(ParticleSet& set, const pat::GenericParticle& dau);
  void makeCandidates();
  void swapDaughters(DoubleMap& swapDauColls, const pat::GenericParticle& cand, const pat::GenericParticleCollection& dauColl);
  void setBestMass(pat::GenericParticle& cand, const DoubleMap& swapDauColls);
  void addSwapCandidates(pat::GenericParticleCollection& swapCandColl, const pat::GenericParticle& cand, const DoubleMap& swapDauColls);
  bool fitCandidate(pat::GenericParticle& cand, const pat::GenericParticleCollection& dauColl);
  void addExtraInfo(pat::GenericParticle& cand);
  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void clear();
  template<class T>
  void clear(std::vector<T>& v) { std::vector<T>().swap(v); };

 private:
  int pdgId_;
  bool doSwap_, matchVertex_, vtxSortByTrkSize_;
  double mass_, width_;
  std::vector<UInt_t> fitAlgoV_;
  std::vector<double> puMap_;
  reco::BeamSpot beamSpot_;
  reco::Vertex beamSpot2D_, vertex_;
  reco::VertexCollection vertices_, priVertices_;
  std::vector<ParticleDaughter> daughters_;
  pat::GenericParticleCollection candidates_, particles_;
  std::map<VertexTuple, reco::VertexRef> vertexRefMap_;
  std::map<ParticleTuple, pat::GenericParticleRef> particleRefMap_;

  reco::VertexRefProd vtxProd_;
  GenericParticleRefProd dauProd_;

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
  edm::EDGetTokenT<reco::ConversionCollection> token_convPhotons_;

  StringCutObjectSelector<pat::GenericParticle, true> preSelection_;
  StringCutObjectSelector<pat::GenericParticle, true> preMassSelection_;
  StringCutObjectSelector<pat::GenericParticle, true> pocaSelection_;
  StringCutObjectSelector<pat::GenericParticle, true> postSelection_;
  StringCutObjectSelector<pat::GenericParticle, true> finalSelection_;
};


#endif
