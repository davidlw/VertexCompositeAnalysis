// -*- C++ -*-
//
// Package:    ParticleProducer
// Class:      ParticleFitter
// 
/**\class ParticleFitter ParticleFitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/ParticleFitter.cc
//
*/

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ParticleFitter.h"


// ParticleFitter: constructor and (empty) destructor
ParticleFitter::ParticleFitter(const edm::ParameterSet& theParameters, edm::ConsumesCollector && iC) :
  pdgId_(theParameters.getParameter<int>("pdgId")),
  doSwap_(theParameters.getParameter<bool>("doSwap")),
  preSelection_(theParameters.getParameter<std::string>("preSelection")),
  preMassSelection_(theParameters.getParameter<std::string>("preMassSelection")),
  pocaSelection_(theParameters.getParameter<std::string>("pocaSelection")),
  postSelection_(theParameters.getParameter<std::string>("postSelection")),
  finalSelection_(theParameters.getParameter<std::string>("finalSelection"))
{
  // get candidate information
  if (theParameters.existsAs<double>("mass")) {
    mass_ = theParameters.getParameter<double>("mass");
  }
  else if (MASS_.find(pdgId_)!=MASS_.end()) {
    mass_ = MASS_.at(pdgId_);
  }
  width_ = mass_*0.15;
  if (theParameters.existsAs<double>("width")) {
    width_ = theParameters.getParameter<double>("width");
  }

  // get daughter information
  const auto daughterVPset = theParameters.getParameter<std::vector<edm::ParameterSet> >("daughterInfo");
  for (const auto& pSet : daughterVPset) {
    ParticleDaughter daughter;
    daughter.fillInfo(pSet, theParameters, iC);
    daughters_.push_back(daughter);
  }

  // get vertex fitter information
  fitAlgoV_  = theParameters.existsAs<std::vector<UInt_t> >("fitAlgo" ) ? theParameters.getParameter<std::vector<UInt_t>>("fitAlgo" ) : std::vector<UInt_t>({0});

  // get input tags
  token_beamSpot_ = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_vertices_ = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("primaryVertices"));
  token_electrons_ = iC.consumes<pat::ElectronCollection>(theParameters.getParameter<edm::InputTag>("electrons"));
  token_muons_ = iC.consumes<pat::MuonCollection>(theParameters.getParameter<edm::InputTag>("muons"));
  token_taus_ = iC.consumes<pat::TauCollection>(theParameters.getParameter<edm::InputTag>("taus"));
  token_photons_ = iC.consumes<pat::PhotonCollection>(theParameters.getParameter<edm::InputTag>("photons"));
  token_tracks_ = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("tracks"));
  token_pfParticles_ = iC.consumes<reco::PFCandidateCollection>(theParameters.getParameter<edm::InputTag>("pfParticles"));
  token_jets_ = iC.consumes<pat::JetCollection>(theParameters.getParameter<edm::InputTag>("jets"));
  token_convPhotons_ = iC.consumes<reco::ConversionCollection>(theParameters.getParameter<edm::InputTag>("conversions"));

  // initialize attributes
  vertex_ = reco::VertexRef();
  beamSpot2D_ = reco::Vertex();
  candidates_ = {};
  particles_ = {};
  vertices_ = {};
};


ParticleFitter::~ParticleFitter()
{
};

// Method containing the algorithm for vertex reconstruction
void ParticleFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get the magnetic field and tracking geometry
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
  iSetup.get<GlobalTrackingGeometryRecord>().get(trkGeometryHandle_);
  // fill daughters particles
  fillDaughters(iEvent);
  // create candidates
  makeCandidates();
};


void ParticleFitter::clear()
{
  vertex_ = reco::VertexRef();
  beamSpot2D_ = reco::Vertex();
  candidates_.clear();
  particles_.clear();
  vertices_.clear();
  std::for_each(daughters_.begin(), daughters_.end(), [](ParticleDaughter &d){ d.clear(); });
};


reco::VertexRef ParticleFitter::getVertexRef(const reco::Vertex& vertex)
{
  const auto& tuple = std::make_tuple(vertex.x(), vertex.y(), vertex.z(), vertex.tracksSize());
  auto& ref = vertexRefMap_[tuple];
  if (ref.isNull()) {
    vertices_.push_back(vertex);
    ref = reco::VertexRef(vtxProd_, vertices_.size()-1);
  }
  return ref;
};


void ParticleFitter::setVertex(const edm::Event& iEvent)
{
  // extract the input collections
  edm::Handle<reco::VertexCollection> vertexHandle;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(token_vertices_, vertexHandle);
  iEvent.getByToken(token_beamSpot_, beamSpotHandle);
  // select the vertices (sort based on track multiplicity)
  reco::VertexCollection vertices;
  for (size_t iVtx=0; iVtx<vertexHandle->size(); iVtx++) {
    const auto& pv = vertexHandle->at(iVtx);
    if (!pv.isFake() && pv.tracksSize() >= 2) {
      vertices.push_back(pv);
      const auto& tuple = std::make_tuple(pv.x(), pv.y(), pv.z(), pv.tracksSize());
      vertexRefMap_[tuple] = reco::VertexRef(vertexHandle, iVtx);
    }
  }
  auto byTracksSize = [] (const reco::Vertex& v1, const reco::Vertex& v2) -> bool { return v1.tracksSize() > v2.tracksSize(); };
  std::sort(vertices.begin(), vertices.end(), byTracksSize);
  // set the primary vertex
  const auto beamSpotVertex = reco::Vertex(beamSpotHandle->position(), beamSpotHandle->rotatedCovariance3D());
  const auto primaryVertex = (vertices.empty() ? beamSpotVertex : vertices[0]);
  vertex_ = getVertexRef(primaryVertex);
  // set the beam spot
  beamSpot_ = *beamSpotHandle;
  // set the 2D beam spot
  const auto beamSpot2DValue = reco::Vertex::Point(beamSpotHandle->x0(), beamSpotHandle->y0(), 0.);
  reco::Vertex::Error beamSpot2DError;
  beamSpot2DError(0,0) = std::pow(beamSpotHandle->BeamWidthX(), 2.);
  beamSpot2DError(1,1) = std::pow(beamSpotHandle->BeamWidthY(), 2.);
  beamSpot2D_ = reco::Vertex(beamSpot2DValue, beamSpot2DError);
};


void ParticleFitter::addParticles(ParticleDaughter& d, const edm::Event& iEvent)
{
  const auto pdgId = std::abs(d.pdgId());
  const auto charge = (d.charge()!=-99 ? d.charge() : HepPDT::ParticleID(pdgId).threeCharge());
  const auto& vertex = (vertex_->isFake() ? beamSpot2D_ : *vertex_);
  if (d.useSource()) { d.addParticles(iEvent); }
  else if (pdgId==0 ) { d.addParticles(iEvent, token_pfParticles_, vertex); }
  else if (pdgId<=6 ) { d.addParticles(iEvent, token_jets_, vertex);        }
  else if (pdgId==11) { d.addParticles(iEvent, token_electrons_, vertex);   }
  else if (pdgId==13) { d.addParticles(iEvent, token_muons_, vertex);       }
  else if (pdgId==15) { d.addParticles(iEvent, token_taus_, vertex);        }
  else if (d.pdgId()== 22) { d.addParticles(iEvent, token_photons_, vertex);     }
  else if (d.pdgId()==-22) { d.addParticles(iEvent, token_convPhotons_, vertex); }
  else if (charge!=0) { d.addParticles(iEvent, token_tracks_, vertex);      }
  else                { d.addParticles(iEvent, token_pfParticles_, vertex); }
};


void ParticleFitter::fillDaughters(const edm::Event& iEvent)
{
  for (auto& daughter : daughters_) {
    addParticles(daughter, iEvent);
  }
  // insert daughter particles
  for (const auto& daughter : daughters_) {
    for (const auto& particle : daughter.particles_) {
      addParticle(particle);
    }
  }
};


pat::GenericParticleRef ParticleFitter::addParticle(const pat::GenericParticle& particle)
{
  const auto& tuple = std::make_tuple(particle.pt(), particle.eta(), particle.phi(), particle.mass(), particle.charge());
  auto& ref = particleRefMap_[tuple];
  if (ref.isNonnull()) { return ref; }
  if (particle.hasUserData("daughters")) {
    auto& daughters = *const_cast<pat::GenericParticleRefVector*>(particle.userData<pat::GenericParticleRefVector>("daughters"));
    pat::GenericParticleRefVector dauColl(dauProd_.id());
    for (const auto& dau : daughters) { dauColl.push_back(addParticle(*dau)); }
    daughters = dauColl;
  }
  if (particle.hasUserData("primaryVertex")) {
    auto& priVtx = *const_cast<reco::VertexRef*>(particle.userData<reco::VertexRef>("primaryVertex"));
    priVtx = getVertexRef(*priVtx);
  }
  particles_.push_back(particle);
  ref = pat::GenericParticleRef(dauProd_, particles_.size()-1);
  return ref;
};


bool ParticleFitter::isUniqueDaughter(ParticleSet& set, const pat::GenericParticle& dau)
{
  if (dau.hasUserData("daughters")) {
    const auto& gDauColl = *dau.userData<pat::GenericParticleRefVector>("daughters");
    for (const auto& gDau : gDauColl) {
      if (!isUniqueDaughter(set, particles_.at(gDau.key()))) { return false; }
    }
  }
  else if (!set.insert(dau).second) { return false; }
  return true;
};


void ParticleFitter::makeCandidates()
{
  // define daughter container
  const auto nDaughters = daughters_.size();
  std::vector<pat::GenericParticleCollection> dauColls;
  for (const auto& daughter : daughters_) {
    if (daughter.particles().empty()) break;
    dauColls.push_back(daughter.particles());
  }
  if (dauColls.size()!=nDaughters) return;
  // make daughter combinations
  const auto combinations = make_combinations(dauColls);
  if (combinations.empty()) return;
  // loop over all combinations of daughters
  std::set<ParticleTuple> candidates;
  for (const auto& combination : combinations) {
    // get unique set of daughters
    ParticleSet daughters;
    for (const auto& daughter : combination) {
      if(!daughters.insert(*daughter).second) break;
    }
    // check if all daughters are unique (are in set)
    if (daughters.size()!=nDaughters) continue;
    // create candidate
    pat::GenericParticle cand;
    int charge = 0;
    reco::Candidate::LorentzVector p4(0,0,0,0);
    for (const auto& daughter : daughters) {
      charge += daughter.charge();
      p4 += daughter.p4();
    }
    cand.setCharge(charge);
    cand.setP4(p4);
    cand.setPdgId(pdgId_);
    cand.setStatus(3);
    // check if candidate has not been found
    const auto& candTuple = std::make_tuple(cand.pt(), cand.eta(), cand.phi(), (doSwap_ ? 0.0 : cand.mass()), cand.charge());
    if (!candidates.insert(candTuple).second) continue;
    // add daughter related information
    float dauPtSum = 0.;
    for (const auto& daughter : daughters) {
      dauPtSum += daughter.pt();
    }
    cand.addUserFloat("dauPtSum", dauPtSum);
    if (nDaughters==2) {
      auto dau1 = daughters.cbegin();
      cand.addUserFloat("dauEtaDiff", std::abs(std::next(dau1)->eta()-dau1->eta()));
    }
    // make preselection
    if (!preSelection_(cand)) continue;
    pat::GenericParticleRefVector dauRefColl(dauProd_.id());
    LorentzVectorColl daughtersP4;
    for (const auto& daughter : daughters) {
      const auto& tuple = std::make_tuple(daughter.pt(), daughter.eta(), daughter.phi(), daughter.mass(), daughter.charge());
      dauRefColl.push_back(particleRefMap_.at(tuple));
      daughtersP4.push_back(daughter.p4());
    }
    cand.addUserData<pat::GenericParticleRefVector>("daughters", dauRefColl);
    cand.addUserData<LorentzVectorColl>("daughtersP4", daughtersP4);
    auto iniCand = cand;
    // check if grandaughters are unique
    bool hasDuplicate = false;
    for (const auto& dau : daughters) {
      if (dau.hasUserData("daughters") && !isUniqueDaughter(daughters, dau)) { hasDuplicate = true; break; }
    }
    if (hasDuplicate) continue;
    // make swapped daughters
    DoubleMap swapDauColls;
    swapDaughters(swapDauColls, iniCand);
    setBestMass(cand, swapDauColls);
    // make mass preselection
    if (!preMassSelection_(cand)) continue;
    // fit candidate and make postselection
    if (!fitCandidate(cand)) continue;
    // add extra information
    addExtraInfo(cand);
    // make final selection
    if (!finalSelection_(cand)) continue;
    // add candidate
    if (std::abs(cand.mass()-mass_) < width_) {
      candidates_.push_back(cand);
    }
    // add swapped candidates
    pat::GenericParticleCollection swapCandColl;
    addSwapCandidates(swapCandColl, cand, swapDauColls);
    for (const auto& swapCand : swapCandColl) {
      if (std::abs(swapCand.mass()-mass_) < width_) {
        candidates_.push_back(swapCand);
      }
    }
  }
  // sort candidates
  std::sort(candidates_.begin(), candidates_.end(), ParticleMassComparator());
};


void ParticleFitter::swapDaughters(DoubleMap& swapDauColls, const pat::GenericParticle& cand)
{
  if (!doSwap_) return;
  // extract daughter collection
  const auto& dauRefColl = *cand.userData<pat::GenericParticleRefVector>("daughters");
  pat::GenericParticleCollection dauColl(dauRefColl.size());
  for (const auto& dau : dauRefColl) { dauColl.push_back(particles_.at(dau.key())); }
  // loop over permutations of daughters
  auto perColl = dauColl;
  while (std::next_permutation(perColl.begin(), perColl.end(), ParticleComparator())) {
    std::vector<double> swapDauColl;
    reco::Candidate::LorentzVector p4(0,0,0,0);
    for (size_t i=0; i<dauColl.size(); i++) {
      const auto& dau = dauColl[i];
      const auto& per = perColl[i];
      const auto& sourceID1 = (dau.hasUserInt("sourceID") ? dau.userInt("sourceID") : 0);
      const auto& sourceID2 = (per.hasUserInt("sourceID") ? per.userInt("sourceID") : 0);
      if (sourceID1!=sourceID2) break;
      if (sourceID1==0 && !ParticleComparator().isParticleEqual(dau, per)) break;
      if (cand.charge()!=0 && dau.charge()!=per.charge()) break;
      swapDauColl.push_back(per.mass());
      p4 += reco::Candidate::LorentzVector(dau.px(), dau.py(), dau.pz(), std::sqrt(dau.p4().P2()+per.massSqr()));
    }
    if (swapDauColl.size()==dauColl.size() && p4.mass()!=cand.mass()) {
      swapDauColls.insert({p4.mass(), swapDauColl});
    }
  }
};


void ParticleFitter::setBestMass(pat::GenericParticle& cand, const DoubleMap& swapDauColls)
{
  if (!doSwap_) return;
  // find mass closest to expected value
  auto mass = cand.mass();
  for (const auto& swapM : swapDauColls) {
    mass = ((std::abs(swapM.first - mass_) < std::abs(mass - mass_)) ? swapM.first : mass);
  }
  cand.addUserFloat("bestMass", mass);
};


void ParticleFitter::addSwapCandidates(pat::GenericParticleCollection& swapCandColl, const pat::GenericParticle& cand, const DoubleMap& swapDauColls)
{
  if (!doSwap_) return;
  // add swapped candidates
  for (const auto& swapM : swapDauColls) {
    pat::GenericParticle swapCand = cand;
    const auto& dauColl = *swapCand.userData<pat::GenericParticleRefVector>("daughters");
    auto& daughtersP4 = *const_cast<LorentzVectorColl*>(swapCand.userData<LorentzVectorColl>("daughtersP4"));
    reco::Candidate::LorentzVector p4(0,0,0,0);
    for (size_t i=0; i<dauColl.size(); i++) {
      const auto& dau = particles_.at(dauColl[i].key());
      const auto energy = std::sqrt(GlobalVector(dau.px(), dau.py(), dau.pz()).mag2() + swapM.second[i]*swapM.second[i]);
      daughtersP4[i] = reco::Candidate::LorentzVector(dau.px(), dau.py(), dau.pz(), energy);
      p4 += daughtersP4[i];
    }
    swapCand.setP4(p4);
    swapCandColl.push_back(swapCand);
  }
};


FitResult ParticleFitter::fitVertex(const KinParColl& particles, const int& fitAlgo, const reco::Vertex& priVtx, const reco::Vertex& decVtx)
{
  if (fitAlgo>2) return FitResult();
  // perform unconstrained vertex fit
  RefCountedKinematicTree result;
  GlobalPoint decP(decVtx.x(), decVtx.y(), decVtx.z());
  if (fitAlgo==0 || (fitAlgo==1 && decP==GlobalPoint())) {
    edm::ParameterSet algoConf;
    algoConf.addParameter("maxDistance", 0.01); // default: 0.01
    algoConf.addParameter("maxNbrOfIterations", 100); // default: 100
    KinematicParticleVertexFitter fitter(algoConf);
    result = fitter.fit(particles);
    if (!result->isValid()) return FitResult();
    decP = result->currentParticle()->currentState().globalPosition();
  }
  // perform constrained vertex fit
  if (fitAlgo>0) {
    edm::ParameterSet algoConf;
    algoConf.addParameter("maxDelta", 0.01); // default: 0.01
    algoConf.addParameter("maxNbrOfIterations", 1000); // default: 1000
    algoConf.addParameter("maxReducedChiSq", 225.); // default: 225.
    algoConf.addParameter("minChiSqImprovement", 50.); // default: 50.
    KinematicConstrainedVertexFitter fitter;
    fitter.setParameters(algoConf);
    if (fitAlgo==1) {
      GlobalPoint pv = GlobalPoint(priVtx.x(), priVtx.y(), priVtx.z());
      MultiTrackPointingKinematicConstraint constrain(pv);
      result = fitter.fit(particles, &constrain, &decP);
    }
    else if (fitAlgo==2) {
      MultiTrackMassKinematicConstraint constrain(mass_, particles.size());
      result = fitter.fit(particles, &constrain, &decP);
    }
  }
  // check result
  if (!result->isValid()) return FitResult();
  const auto& fitVertex = result->currentDecayVertex();
  if (!fitVertex->vertexIsValid()) return FitResult();
  const auto& candState = result->currentParticle()->currentState();
  if (!candState.isValid()) return FitResult();
  // extract information
  std::vector<GlobalVector> fitMom;
  for (const auto& particle : result->daughterParticles()) {
    if (!particle->currentState().isValid()) return FitResult();
    fitMom.push_back(particle->currentState().kinematicParameters().momentum());
  }
  const auto decVtxP = reco::Particle::Point(fitVertex->position().x(), fitVertex->position().y(), fitVertex->position().z());
  const auto fitVtx = reco::Vertex(decVtxP, fitVertex->error().matrix(), fitVertex->chiSquared(), fitVertex->degreesOfFreedom(), particles.size());
  return std::make_tuple(fitMom, fitVtx, candState);
};


FitResult ParticleFitter::fitVertex(const TransTrackColl& particles, const int& fitAlgo)
{
  if (fitAlgo<3) return FitResult();
  // perform vertex fit
  CachingVertex<5> fitVertex;
  if (fitAlgo==3) {
    edm::ParameterSet algoConf;
    algoConf.addParameter("maxDistance", 0.01); // default: 0.01
    algoConf.addParameter("maxNbrOfIterations", 10); // default: 10
    KalmanVertexFitter fitter(algoConf, true);
    fitVertex = fitter.vertex(particles);
  }
  else if (fitAlgo==4) {
    edm::ParameterSet algoConf;
    algoConf.addParameter("maxDistance", 0.01); // default: 0.01
    algoConf.addParameter("maxNbrOfIterations", 10); // default: 10
    algoConf.addParameter("limitComponents", false);
    algoConf.addParameter("smoothTracks", true);
    GsfVertexFitter fitter(algoConf);
    fitVertex = fitter.vertex(particles);
  }
  else if (fitAlgo==5) {
    edm::ParameterSet algoConf;
    algoConf.addParameter("maxshift", 0.0001); // default: 0.0001
    algoConf.addParameter("maxlpshift", 0.1); // default: 0.1
    algoConf.addParameter("maxstep", 30); // default: 30
    algoConf.addParameter("weightthreshold", 0.001); // default 0.001
    auto fitter = AdaptiveVertexFitter(GeometricAnnealing(), DefaultLinearizationPointFinder(), KalmanVertexUpdator<5>(), KalmanVertexTrackCompatibilityEstimator<5>(), KalmanVertexSmoother());
    fitter.setParameters(algoConf);
    fitVertex = fitter.vertex(particles);
  }
  else if (fitAlgo==6) {
    edm::ParameterSet algoConf;
    algoConf.addParameter("maxshift", 0.0001); // default: 0.0001
    algoConf.addParameter("maxlpshift", 0.1); // default: 0.1
    algoConf.addParameter("maxstep", 30); // default: 30
    algoConf.addParameter("weightthreshold", 0.001); // default 0.001
    algoConf.addParameter("limitComponents", false);
    AdaptiveGsfVertexFitter fitter(algoConf);
    fitVertex = fitter.vertex(particles);
  }
  // check result
  if (!fitVertex.isValid() || !fitVertex.vertexState().isValid()) return FitResult();
  // extract information
  std::vector<GlobalVector> fitMom;
  for (const auto& track : fitVertex.tracks()) {
    if (!track->refittedStateAvailable()) return FitResult();
    fitMom.push_back(track->refittedState()->freeTrajectoryState().momentum());
  }
  const auto nDoF = std::max(std::round(fitVertex.degreesOfFreedom()), 0.f);
  const auto decVtxP = reco::Particle::Point(fitVertex.position().x(), fitVertex.position().y(), fitVertex.position().z());
  const auto fitVtx = reco::Vertex(decVtxP, fitVertex.error().matrix(), fitVertex.totalChiSquared(), nDoF, particles.size());
  return std::make_tuple(fitMom, fitVtx, KinematicState());
};
  

bool ParticleFitter::fitCandidate(pat::GenericParticle& cand)
{
  if (fitAlgoV_.empty() || !cand.hasUserData("daughters") || daughters_.size()<2) return true;
  const auto& magField = bFieldHandle_.product();
  // get the daughters transient tracks
  const auto& dauColl = *cand.userData<pat::GenericParticleRefVector>("daughters");
  std::vector<std::pair<std::vector<std::pair<reco::TransientTrack, double> >, size_t> > daughters;
  for (size_t iDau=0; iDau<dauColl.size(); iDau++) {
    const auto& daughter = particles_.at(dauColl[iDau].key());
    std::vector<std::pair<reco::TransientTrack, double> > tracks;
    if (daughter.numberOfTracks()>0) {
      for (size_t iTrk=0; iTrk<daughter.numberOfTracks(); iTrk++) {
        const auto& trackMass = daughter.userFloat(Form("trackMass%lu", iTrk));
        tracks.push_back({reco::TransientTrack(*daughter.track(iTrk), magField, trkGeometryHandle_), trackMass});
      }
    }
    else if (daughter.track().isNonnull()) {
      tracks.push_back({reco::TransientTrack(*daughter.track(), magField, trkGeometryHandle_), daughter.mass()});
    }
    else if (daughter.hasUserData("kinematicParameters")) {
      tracks.push_back({reco::TransientTrack(), daughter.mass()});
    }
    if (!tracks.empty()) {
      daughters.push_back({tracks, iDau});
    }
  }
  if (daughters.size()<2) return true;
  // measure distance between daughter tracks at their point of closest approach
  if (daughters.size()==2 && daughters[0].first.size()==1 && daughters[1].first.size()==1 &&
      daughters[0].first[0].first.isValid() && daughters[0].first[0].first.impactPointTSCP().isValid() &&
      daughters[1].first[0].first.isValid() && daughters[1].first[0].first.impactPointTSCP().isValid()) {
    ClosestApproachInRPhi cApp;
    const auto& stateDau1 = daughters[0].first[0].first.impactPointTSCP().theState();
    const auto& stateDau2 = daughters[1].first[0].first.impactPointTSCP().theState();
    cApp.calculate(stateDau1, stateDau2);
    if (!cApp.status()) return false;
    const auto& cxPt = cApp.crossingPoint();
    reco::Candidate::LorentzVector p4(0,0,0,0), swapP4(0,0,0,0);
    for (size_t i=0; i<2; i++) {
      const auto& tscpDau = daughters[i].first[0].first.trajectoryStateClosestToPoint(cxPt);
      if (!tscpDau.isValid()) return false;
      for (size_t j=0; j<2; j++) {
        if (!doSwap_ && i!=j) continue;
        const auto& daughter = particles_.at(dauColl[daughters[j].second].key());
        const auto dauEnergy = std::sqrt(tscpDau.momentum().mag2() + daughter.massSqr());
        const auto tscpDauP4 = reco::Candidate::LorentzVector(tscpDau.momentum().x(), tscpDau.momentum().y(), tscpDau.momentum().z(), dauEnergy);
        if (i==j) { p4 += tscpDauP4; }
        else { swapP4 += tscpDauP4; }
      }
    }
    cand.setP4(p4);
    if (doSwap_) {
      const auto& mass = reco::Candidate::PolarLorentzVector(swapP4).M();
      cand.addUserFloat("bestMass", ((std::abs(cand.mass() - mass_) < std::abs(mass - mass_)) ? cand.mass() : mass), true);
    }
    cand.addUserFloat("dca", cApp.distance());
    if (!pocaSelection_(cand)) return false;
  }
  // prepare particles for decay vertex fit
  std::map<size_t, std::vector<size_t> > dauParIdx;
  std::vector<double> masses;
  TransTrackColl tracks;
  KinParColl particles;
  for (const auto& d : daughters) {
    const auto& daughter = particles_.at(dauColl[d.second].key());
    std::vector<size_t> parIdx;
    if (d.first[0].first.isValid()) {
      for (const auto& t : d.first) {
        if (!t.first.isValid() || !t.first.impactPointTSCP().isValid()) return false;
        float chi = 0., ndf = 0., width = daughter.userFloat("width");
        KinematicParticleFactoryFromTransientTrack pFactory;
        tracks.push_back(t.first);
        masses.push_back(t.second);
        particles.push_back(pFactory.particle(t.first, t.second, chi, ndf, width));
        parIdx.push_back(particles.size()-1);
      }
    }
    else {
      float chi = 0., ndf = 0.;
      const auto& kinPar = *daughter.userData<KinematicParameters>("kinematicParameters");
      const auto& kinParError = *daughter.userData<KinematicParametersError>("kinematicParametersError");
      const auto state = KinematicState(kinPar, kinParError, daughter.charge(), magField);
      VirtualKinematicParticleFactory pFactory;
      masses.push_back(d.first[0].second);
      particles.push_back(pFactory.particle(state, chi, ndf, NULL));
      parIdx.push_back(particles.size()-1);
    }
    if (!parIdx.empty()) {
      dauParIdx[d.second] = parIdx;
    }
  }
  if (dauParIdx.size()<2) return false;
  // fit decay vertex
  FitResult fit0;
  for (const auto& a : fitAlgoV_) {
    if (cand.hasUserFloat(Form("vertexProb_%d",a))) continue;
    if (a>2 && tracks.size()!=particles.size()) throw std::logic_error(Form("Invalid fit algorithm (%d) for virtual daughters", a));
    FitResult fitResult;
    if (a<3) {
      if (std::get<0>(fit0).empty()) { fit0 = fitVertex(particles, 0); }
      fitResult = (a==0 ? fit0 : fitVertex(particles, a, *vertex_, std::get<1>(fit0)));
    }
    else { fitResult = fitVertex(tracks, a); }
    // extract fit result
    const auto& dauMom = std::get<0>(fitResult);
    const auto& decVtx = std::get<1>(fitResult);
    const auto& candSt = std::get<2>(fitResult);
    // add basic fit information
    if (!dauMom.empty()) {
      const auto lbl = (a==*fitAlgoV_.begin() ? "" : Form("_%d",a));
      cand.addUserFloat(Form("normChi2%s",lbl), decVtx.chi2()/decVtx.ndof());
      cand.addUserFloat(Form("vertexProb%s",lbl), TMath::Prob(decVtx.chi2(), decVtx.ndof()));
      if (candSt.isValid()) {
        cand.addUserData<KinematicParameters>(Form("kinematicParameters%s",lbl), candSt.kinematicParameters());
      }
    }
    if (a!=fitAlgoV_[0]) continue;
    if (dauMom.empty()) return false;
    // update particle kinematics from fit
    reco::Candidate::LorentzVector candP4(0, 0, 0, 0);
    auto& daughtersP4 = *const_cast<LorentzVectorColl*>(cand.userData<LorentzVectorColl>("daughtersP4"));
    for (size_t iDau=0; iDau<dauColl.size(); iDau++) {
      if (dauParIdx.find(iDau)!=dauParIdx.end()) {
        reco::Candidate::LorentzVector dauP4(0, 0, 0, 0);
        for (const auto& iPar : dauParIdx.at(iDau)) {
          const auto& fitMom = dauMom[iPar];
          const auto fitEnergy = std::sqrt(fitMom.mag2() + masses[iPar]*masses[iPar]);
          dauP4 += reco::Candidate::LorentzVector(fitMom.x(), fitMom.y(), fitMom.z(), fitEnergy);
        }
        daughtersP4[iDau] = dauP4;
      }
      candP4 += daughtersP4[iDau];
    }
    cand.setP4(candP4);
    // apply post-fit selection
    cand.setVertex(decVtx.position());
    if (!postSelection_(cand)) return false;
    // add information to fitted candidate
    cand.addUserData<reco::Vertex>("decayVertex", decVtx);
    if (candSt.isValid()) {
      cand.addUserData<KinematicParametersError>("kinematicParametersError", candSt.kinematicParametersError());
    }
  }
  cand.addUserData<reco::VertexRef>("primaryVertex", vertex_);
  return true;
};


void ParticleFitter::addExtraInfo(pat::GenericParticle& cand)
{
  if (!cand.hasUserData("decayVertex") || !cand.hasUserData("primaryVertex")) return;
  const auto& decayVertex = *cand.userData<reco::Vertex>("decayVertex");
  const auto& priVtx = *cand.userData<reco::VertexRef>("primaryVertex");
  const auto& primaryVertex = (priVtx.isAvailable() ? *priVtx : vertices_.at(priVtx.key()));
  // compute lifetime information
  const auto lineOfFlight = decayVertex.position() - primaryVertex.position();
  const auto lVtxMag = lineOfFlight.r();
  const auto rVtxMag = lineOfFlight.rho();
  const auto angle3D = angle(lineOfFlight.x(), lineOfFlight.y(), lineOfFlight.z(), cand.px(), cand.py(), cand.pz());
  const auto angle2D = angle(lineOfFlight.x(), lineOfFlight.y(), 0.0, cand.px(), cand.py(), 0.0);
  const auto distanceVector3D = SVector3(lineOfFlight.x(), lineOfFlight.y(), lineOfFlight.z());
  const auto distanceVector2D = SVector3(lineOfFlight.x(), lineOfFlight.y(), 0.0);
  const auto totalCov = decayVertex.covariance() + primaryVertex.covariance();
  const auto sigmaLvtxMag = std::sqrt(ROOT::Math::Similarity(totalCov, distanceVector3D)) / lVtxMag;
  const auto sigmaRvtxMag = std::sqrt(ROOT::Math::Similarity(totalCov, distanceVector2D)) / rVtxMag;
  // set lifetime infomation
  cand.addUserFloat("lVtxMag", lVtxMag);
  cand.addUserFloat("rVtxMag", rVtxMag);
  cand.addUserFloat("lVtxSig", (lVtxMag/sigmaLvtxMag));
  cand.addUserFloat("rVtxSig", (rVtxMag/sigmaRvtxMag));
  cand.addUserFloat("angle3D", angle3D);
  cand.addUserFloat("angle2D", angle2D);
};


// ParticleDaughter: constructor and (empty) destructor
ParticleDaughter::ParticleDaughter()
{
  pdgId_ = 0;
  charge_ = -99;
  mass_ = 0.;
  width_ = 0.;
  selection_ = "";
  finalSelection_ = "";
  particles_ = {}; 
};


ParticleDaughter::ParticleDaughter(const edm::ParameterSet& pSet, const edm::ParameterSet& config, edm::ConsumesCollector&& iC) :
  ParticleDaughter()
{
  fillInfo(pSet, config, iC);
};


ParticleDaughter::~ParticleDaughter()
{
};

void ParticleDaughter::clear() {
  particles_.clear();
};


void ParticleDaughter::fillInfo(const edm::ParameterSet& pSet, const edm::ParameterSet& config, edm::ConsumesCollector& iC)
{
  if (pSet.existsAs<int>("pdgId")) {
    pdgId_ = pSet.getParameter<int>("pdgId");
  }
  if (pSet.existsAs<int>("charge")) {
    charge_ = pSet.getParameter<int>("charge");
  }
  if (pSet.existsAs<double>("mass")) {
    mass_ = pSet.getParameter<double>("mass");
  }
  else if (MASS_.find(pdgId_)!=MASS_.end()) {
    mass_ = MASS_.at(pdgId_);
  }
  if (pSet.existsAs<double>("width")) {
    width_ = pSet.getParameter<double>("width");
  }
  else if (WIDTH_.find(pdgId_)!=WIDTH_.end()) {
    width_ = WIDTH_.at(pdgId_);
  }
  else { width_ = mass_*1e-6; }
  if (pSet.existsAs<std::string>("selection")) {
    selection_ = pSet.getParameter<std::string>("selection");
  }
  if (pSet.existsAs<std::string>("finalSelection")) {
    finalSelection_ = pSet.getParameter<std::string>("finalSelection");
  }
  if (pSet.existsAs<edm::InputTag>("source")) {
    token_source_ = iC.consumes<pat::GenericParticleCollection>(pSet.getParameter<edm::InputTag>("source"));
  }
  const auto& mvaSet = (pSet.existsAs<edm::InputTag>("mva") ? pSet : config);
  if (mvaSet.existsAs<edm::InputTag>("mva")) {
    token_mva_ = iC.consumes<std::vector<float> >(mvaSet.getParameter<edm::InputTag>("mva"));
  }
  if (config.existsAs<edm::InputTag>("dedxHarmonic2")) {
    token_dedx_ = iC.consumes<edm::ValueMap<reco::DeDxData> >(config.getParameter<edm::InputTag>("dedxHarmonic2"));
  }
  if (pSet.existsAs<edm::InputTag>("muonL1Info")) {
    token_muonL1Info_ = iC.consumes<pat::TriggerObjectStandAloneMatch>(pSet.getParameter<edm::InputTag>("muonL1Info"));
  }
};


template <class T>
void ParticleDaughter::addParticles(const edm::Event& event, const edm::EDGetTokenT<std::vector<T> >& token, const reco::Vertex& vertex, const bool embedInfo)
{
  // extract input collections
  edm::Handle<std::vector<T> > handle;
  if (!token.isUninitialized()) event.getByToken(token, handle);
  edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxMap;
  if (!token_dedx_.isUninitialized()) event.getByToken(token_dedx_, dEdxMap);
  edm::Handle<std::vector<float> > mvaColl;
  if (!token_mva_.isUninitialized()) event.getByToken(token_mva_, mvaColl);
  edm::Handle<pat::TriggerObjectStandAloneMatch> muonL1Info;
  if (!token_muonL1Info_.isUninitialized()) event.getByToken(token_muonL1Info_, muonL1Info);
  // set selections
  StringCutObjectSelector<T, true> selection(selection_);
  StringCutObjectSelector<pat::GenericParticle, true> finalSelection(finalSelection_);
  // add particles
  if (handle.isValid()) {
    ParticleMassSet particles;
    for (size_t i=0; i<handle->size(); i++) {
      const auto& p = edm::Ref<std::vector<T> >(handle, i);
      if (!selection(*p)) continue;
      pat::GenericParticle cand;
      addInfo(cand, *p);
      if (charge_!=-99 && cand.charge()!=charge_) continue;
      cand.setPdgId(pdgId_);
      cand.setStatus(1);
      addData(cand, p, embedInfo);
      if (cand.track().isNonnull()) {
        const auto& trk = *cand.track();
        const float dz = trk.dz(vertex.position());
        const float dxy = trk.dxy(vertex.position());
        const float dzerror  = std::sqrt(trk.dzError()*trk.dzError()+vertex.zError()*vertex.zError());
        const float dxyerror = std::sqrt(trk.d0Error()*trk.d0Error()+vertex.xError()*vertex.yError());
        cand.addUserFloat("dz", dz);
        cand.addUserFloat("dxy", dxy);
        cand.addUserFloat("dzSig", dz/dzerror);
        cand.addUserFloat("dxySig", dxy/dxyerror);
      }
      else if (cand.charge()!=0) continue; // ignore if charged particle has no track
      cand.addUserFloat("width", width_);
      setDeDx(cand, dEdxMap);
      setMVA(cand, i, mvaColl);
      addMuonL1Info(cand, muonL1Info);
      if (finalSelection(cand)) {
        particles.insert(cand);
      }
    }
    particles_.assign(particles.begin(), particles.end());
  }
};


void ParticleDaughter::addParticles(const edm::Event& event)
{
  edm::Handle<pat::GenericParticleCollection> handle;
  event.getByToken(token_source_, handle);
  StringCutObjectSelector<pat::GenericParticle, true> selection(selection_);
  StringCutObjectSelector<pat::GenericParticle, true> finalSelection(finalSelection_);
  if (handle.isValid()) {
    ParticleMassSet particles;
    for (const auto& p : *handle) {
      if (!selection(p)) continue;
      pat::GenericParticle cand(p);
      if (charge_!=-99 && cand.charge()!=charge_) continue;
      if (cand.hasUserData("daughters")) {
        cand.setStatus(2);
      }
      if (finalSelection(cand)) {
        particles.insert(cand);
      }
    }
    particles_.assign(particles.begin(), particles.end());
  }
};


template <class T>
void ParticleDaughter::addInfo(pat::GenericParticle& c, const T& p)
{
  const auto p4 = reco::Candidate::PolarLorentzVector(p.pt(), p.eta(), p.phi(), mass_);
  c.setP4(p4);
  c.setCharge(p.charge());
  c.setVertex(p.vertex());
};


void ParticleDaughter::addInfo(pat::GenericParticle& c, const reco::Conversion& p)
{
  const auto& vtx = p.conversionVertex();
  reco::TrackCollection tracks;
  if (vtx.isValid() && vtx.nTracks(0.5)==2) {
    tracks = vtx.refittedTracks();
  }
  else {
    for (const auto& t: p.tracks()) { tracks.push_back(*t); }
  }
  int charge = 0;
  reco::Particle::PolarLorentzVector p4(0, 0, 0, 0);
  for (const auto& t: tracks) {
    charge += t.charge();
    p4 += reco::Particle::PolarLorentzVector(t.pt(), t.eta(), t.phi(), 0.000511);
  }
  reco::TrackRefVector trackRefs;
  for (size_t i=0; i<tracks.size(); i++) {
    trackRefs.push_back(reco::TrackRef(&tracks, i));
    c.addUserFloat(Form("trackMass%lu", i), 0.000511);
  }
  c.setP4(p4);
  c.setCharge(charge);
  c.setVertex(vtx.position());
  c.setTracks(trackRefs, true);
};


template <class T>
void ParticleDaughter::addData(pat::GenericParticle& c, const edm::Ref<std::vector<T> >& p, const bool& embedInfo)
{
  c.addUserData<T>("src", *p);
};


void ParticleDaughter::addData(pat::GenericParticle& c, const reco::TrackRef& p, const bool& embedInfo)
{
  c.setTrack(p, embedInfo);
  if (embedInfo) c.addUserData<reco::TrackRef>("trackRef", p);
  c.addUserInt("sourceID", 1);
};


void ParticleDaughter::addData(pat::GenericParticle& c, const reco::PFCandidateRef& p, const bool& embedInfo)
{
  if (c.pdgId()==0) { c.setPdgId(p->pdgId()); }
  c.setTrack(p->trackRef(), embedInfo);
  if (embedInfo) c.addUserData<reco::TrackRef>("trackRef", p->trackRef());
  c.addUserData<reco::PFCandidate>("src", *p);
};


void ParticleDaughter::addData(pat::GenericParticle& c, const pat::MuonRef& p, const bool& embedInfo)
{
  auto track = p->track();
  if (!track.id().isValid()) { track = dynamic_cast<const reco::Muon*>(p->originalObject())->track(); }
  c.setTrack(track, embedInfo);
  if (embedInfo) c.addUserData<reco::TrackRef>("trackRef", track);
  c.addUserData<pat::Muon>("src", *p);
};


void ParticleDaughter::addData(pat::GenericParticle& c, const pat::ElectronRef& p, const bool& embedInfo)
{
  const auto& trkC = reco::TrackCollection({*dynamic_cast<const reco::Track*>(p->gsfTrack().get())});
  c.setTrack(reco::TrackRef(&trkC, 0), true);
  c.addUserData<pat::Electron>("src", *p);
};


void ParticleDaughter::setMVA(pat::GenericParticle& c, const size_t& i, const edm::Handle<std::vector<float> >& mvaColl)
{
  if (mvaColl.isValid() && i<mvaColl->size()) {
    c.addUserFloat("mva", mvaColl->at(i));
  }
};


void ParticleDaughter::setDeDx(pat::GenericParticle& c, const edm::Handle<edm::ValueMap<reco::DeDxData> >& dEdxMap)
{
  const auto& track = (c.hasUserData("trackRef") ? *c.userData<reco::TrackRef>("trackRef") : c.track());
  if (track.isNull()) return;
  if (dEdxMap.isValid() && dEdxMap->contains(track.id())) {
    const auto& dEdx = (*dEdxMap)[track].dEdx();
    c.addUserFloat("dEdx", dEdx);
  }
};


void ParticleDaughter::addMuonL1Info(pat::GenericParticle& c, const edm::Handle<pat::TriggerObjectStandAloneMatch>& muonL1Info)
{
  if (std::abs(c.pdgId())!=13 || !c.userData<pat::Muon>("src")) return;
  const auto& p = c.userData<pat::Muon>("src")->originalObjectRef();
  if (muonL1Info.isValid() && muonL1Info->contains(p.id())) {
    const auto& muonP4 = (*muonL1Info)[p]->p4();
    c.addUserData<reco::Particle::LorentzVector>("muonL1", muonP4);
  }
};
