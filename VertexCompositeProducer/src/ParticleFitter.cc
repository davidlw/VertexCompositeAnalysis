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

  // initialize attributes
  vertex_ = reco::Vertex();
  beamSpot2D_ = reco::Vertex();
  candidates_ = {};
};


ParticleFitter::~ParticleFitter() {
};

// Method containing the algorithm for vertex reconstruction
void ParticleFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // set the magnetic field
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
  // fill daughters particles
  fillDaughters(iEvent);
  // create candidates
  makeCandidates();
};


void ParticleFitter::clear() {
  vertex_ = reco::Vertex();
  beamSpot2D_ = reco::Vertex();
  candidates_.clear();
  std::for_each(daughters_.begin(), daughters_.end(), [](ParticleDaughter &d){ d.clear(); });
};


void ParticleFitter::setVertex(const edm::Event& iEvent) {
  // extract the input collections
  edm::Handle<reco::VertexCollection> vertexHandle;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(token_vertices_, vertexHandle);
  iEvent.getByToken(token_beamSpot_, beamSpotHandle);
  // set the primary vertex
  const auto beamSpotVertex = reco::Vertex(beamSpotHandle->position(), beamSpotHandle->rotatedCovariance3D());
  const auto primaryVertex = (vertexHandle->size()>0 ? vertexHandle->at(0) : reco::Vertex());
  const bool isGoodVtx = !primaryVertex.isFake() && primaryVertex.tracksSize()>=2;
  vertex_ = (isGoodVtx ? primaryVertex : beamSpotVertex);
  // set the 2D beam spot
  const auto beamSpot2DValue = reco::Vertex::Point(beamSpotHandle->position().x(), beamSpotHandle->position().y(), 0.);
  reco::Vertex::Error beamSpot2DError;
  beamSpot2DError(0,0) = std::pow(beamSpotHandle->BeamWidthX(), 2.);
  beamSpot2DError(1,1) = std::pow(beamSpotHandle->BeamWidthY(), 2.);
  beamSpot2D_ = reco::Vertex(beamSpot2DValue, beamSpot2DError);
};


void ParticleFitter::addParticles(ParticleDaughter& d, const edm::Event& iEvent) {
  const auto pdgId = std::abs(d.pdgId());
  const auto charge = d.charge();
  const auto& vertex = (vertex_.isFake() ? beamSpot2D_ : vertex_);
  if (d.useSource()) { d.addParticles(iEvent); }
  else if (pdgId==0 ) { d.addParticles(iEvent, token_pfParticles_, vertex); }
  else if (pdgId<=6 ) { d.addParticles(iEvent, token_jets_, vertex);        }
  else if (pdgId==11) { d.addParticles(iEvent, token_electrons_, vertex);   }
  else if (pdgId==13) { d.addParticles(iEvent, token_muons_, vertex);       }
  else if (pdgId==15) { d.addParticles(iEvent, token_taus_, vertex);        }
  else if (pdgId==22) { d.addParticles(iEvent, token_photons_, vertex);     }
  else if (charge!=0) { d.addParticles(iEvent, token_tracks_, vertex);      }
  else                { d.addParticles(iEvent, token_pfParticles_, vertex); }
};


void ParticleFitter::fillDaughters(const edm::Event& iEvent) {
  for (auto& daughter : daughters_) {
    addParticles(daughter, iEvent);
  }
};


bool ParticleFitter::isUniqueDaughter(ParticleSet& set, const pat::GenericParticle& dau) {
  if (dau.hasUserData("daughters")) {
    const auto& gDauColl = *dau.userData<pat::GenericParticleCollection>("daughters");
    for (const auto& gDau : gDauColl) {
      if (!isUniqueDaughter(set, gDau)) { return false; }
    }
  }
  else if (!set.insert(dau).second) { return false; }
  return true;
};


void ParticleFitter::makeCandidates() {
  // define daughter container
  const auto nDaughters = daughters_.size();
  std::vector<pat::GenericParticleCollection> daughterColls;
  for (const auto& daughter : daughters_) {
    if (daughter.particles().empty()) break;
    daughterColls.push_back(daughter.particles());
  }
  if (daughterColls.size()!=nDaughters) return;
  // make daughter combinations
  const auto combinations = make_combinations(daughterColls);
  if (combinations.empty()) return;
  // loop over all combinations of daughters
  std::set<std::tuple<float, float, float, float, signed char> > candidates;
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
    pat::GenericParticleCollection daughterColl(daughters.begin(), daughters.end());
    cand.addUserData<pat::GenericParticleCollection>("daughters", daughterColl);
    auto iniCand = cand;
    // check if grandaughters are unique
    bool hasDuplicate = false;
    for (const auto& dau : daughterColl) {
      if (dau.hasUserData("daughters") && !isUniqueDaughter(daughters, dau)) { hasDuplicate = true; break; }
    }
    if (hasDuplicate) continue;
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
    // make swapped daughters
    DoubleMap swapDauColls;
    swapDaughters(swapDauColls, iniCand);
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
  const auto& dauColl = *cand.userData<pat::GenericParticleCollection>("daughters");
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
  cand.addUserFloat("mass", mass);
};


void ParticleFitter::addSwapCandidates(pat::GenericParticleCollection& swapCandColl, const pat::GenericParticle& cand, const DoubleMap& swapDauColls)
{
  if (!doSwap_) return;
  // add swapped candidates
  for (const auto& swapM : swapDauColls) {
    pat::GenericParticle swapCand = cand;
    auto& dauColl = *const_cast<pat::GenericParticleCollection*>(swapCand.userData<pat::GenericParticleCollection>("daughters"));
    reco::Candidate::LorentzVector p4(0,0,0,0);
    for (size_t i=0; i<dauColl.size(); i++) {
      auto& dau = dauColl[i];
      const auto energy = std::sqrt(GlobalVector(dau.px(), dau.py(), dau.pz()).mag2() + swapM.second[i]*swapM.second[i]);
      dau.setP4(reco::Candidate::LorentzVector(dau.px(), dau.py(), dau.pz(), energy));
      p4 += dau.p4();
    }
    swapCand.setP4(p4);
    swapCandColl.push_back(swapCand);
  }
};


bool ParticleFitter::fitCandidate(pat::GenericParticle& cand) {
  if (!cand.hasUserData("daughters") || daughters_.size()<2) return true;
  const auto iniP4 = cand.p4();
  const auto& magField = bFieldHandle_.product();
  // get the daughters transient tracks
  auto& daughterColl = *const_cast<pat::GenericParticleCollection*>(cand.userData<pat::GenericParticleCollection>("daughters"));
  std::vector<std::pair<reco::TransientTrack, size_t> > daughters;
  for (size_t iDau=0; iDau<daughterColl.size(); iDau++) {
    const auto& trk = daughterColl[iDau].track();
    if (trk.isNonnull()) {
      daughters.push_back({reco::TransientTrack(*trk, magField), iDau});
    }
    else if (daughterColl[iDau].hasUserData("kinematicParameters")) {
      daughters.push_back({reco::TransientTrack(), iDau});
    }
  }
  if (daughters.size()<2) return true;
  // measure distance between daughter tracks at their point of closest approach
  if (daughters.size()==2 &&
      daughters[0].first.isValid() && daughters[0].first.impactPointTSCP().isValid() &&
      daughters[1].first.isValid() && daughters[1].first.impactPointTSCP().isValid()) {
    ClosestApproachInRPhi cApp;
    const auto& stateDau1 = daughters[0].first.impactPointTSCP().theState();
    const auto& stateDau2 = daughters[1].first.impactPointTSCP().theState();
    cApp.calculate(stateDau1, stateDau2);
    if (!cApp.status()) return false;
    const auto& cxPt = cApp.crossingPoint();
    reco::Candidate::LorentzVector p4(0,0,0,0), swapP4(0,0,0,0);
    for (size_t i=0; i<2; i++) {
      const auto& tscpDau = daughters[i].first.trajectoryStateClosestToPoint(cxPt);
      if (!tscpDau.isValid()) return false;
      for (size_t j=0; j<2; j++) {
        if (!doSwap_ && i!=j) continue;
        const auto dauEnergy = std::sqrt(tscpDau.momentum().mag2() + daughterColl[daughters[j].second].massSqr());
        const auto tscpDauP4 = reco::Candidate::LorentzVector(tscpDau.momentum().x(), tscpDau.momentum().y(), tscpDau.momentum().z(), dauEnergy);
        if (i==j) { p4 += tscpDauP4; }
        else { swapP4 += tscpDauP4; }
      }
    }
    cand.setP4(p4);
    if (doSwap_) {
      const auto& mass = reco::Candidate::PolarLorentzVector(swapP4).M();
      cand.addUserFloat("mass", ((std::abs(cand.mass() - mass_) < std::abs(mass - mass_)) ? cand.mass() : mass), true);
    }
    cand.addUserFloat("dca", cApp.distance());
    if (!pocaSelection_(cand)) return false;
  }
  // prepare particles for decay vertex fit
  std::vector<RefCountedKinematicParticle> particles;
  for (const auto& d : daughters) {
    const auto& daughter = daughterColl[d.second];
    if (d.first.isValid()) {
      if (!d.first.impactPointTSCP().isValid()) return false;
      float chi = 0., ndf = 0., width = daughter.userFloat("width");
      KinematicParticleFactoryFromTransientTrack pFactory;
      particles.push_back(pFactory.particle(d.first, daughter.mass(), chi, ndf, width));
    }
    else {
      float chi = 0., ndf = 0.;
      const auto& kinPar = *daughter.userData<KinematicParameters>("kinematicParameters");
      const auto& kinParError = *daughter.userData<KinematicParametersError>("kinematicParametersError");
      const auto state = KinematicState(kinPar, kinParError, daughter.charge(), magField);
      VirtualKinematicParticleFactory pFactory;
      particles.push_back(pFactory.particle(state, chi, ndf, NULL));
    }
  }
  if (particles.size()<2) return false;
  // fit decay vertex
  KinematicParticleVertexFitter fitter;
  const auto& result = fitter.fit(particles);
  if (!result->isValid()) return false;
  const auto& fitVertex = result->currentDecayVertex();
  if (!fitVertex->vertexIsValid()) return false;
  // check fitted final state particle
  const auto& candState = result->currentParticle()->currentState();
  if (!candState.isValid()) return false;
  // check fitted daughter particles
  std::vector<RefCountedKinematicParticle> daughterParticles;
  for (const auto& p : result->daughterParticles()) {
    if (p->currentState().isValid()) {
      daughterParticles.push_back(p);
    }
  }
  if (daughterParticles.size()!=daughters.size()) return false;
  // update particle kinematics from fit
  cand.addUserData<reco::Candidate::LorentzVector>("initialP4", iniP4);
  reco::Candidate::LorentzVector candP4(0, 0, 0, 0);
  for (size_t iDau=0; iDau<daughters.size(); iDau++) {
    auto& daughter = daughterColl[daughters[iDau].second];
    const auto iniP4 = daughter.p4();
    const auto fitKP = daughterParticles[iDau]->currentState().kinematicParameters();
    const auto fitEnergy = std::sqrt(fitKP.momentum().mag2() + daughter.massSqr());
    const auto fitP4 = reco::Candidate::LorentzVector(fitKP.momentum().x(), fitKP.momentum().y(), fitKP.momentum().z(), fitEnergy);
    daughter.addUserData<reco::Candidate::LorentzVector>("initialP4", iniP4);
    daughter.setP4(fitP4);
    candP4 += fitP4;
  }
  cand.setP4(candP4);
  // set data
  const auto decayVertexPos = reco::Particle::Point(fitVertex->position().x(), fitVertex->position().y(), fitVertex->position().z());
  reco::Vertex decayVertex(decayVertexPos, fitVertex->error().matrix(), fitVertex->chiSquared(), fitVertex->degreesOfFreedom(), daughters.size());
  cand.setVertex(decayVertexPos);
  cand.addUserFloat("normChi2", decayVertex.chi2()/decayVertex.ndof());
  cand.addUserFloat("vertexProb", TMath::Prob(decayVertex.chi2(), decayVertex.ndof()));
  // add information to fitted candidate
  if (!postSelection_(cand)) return false;
  cand.addUserData<reco::Vertex>("primaryVertex", vertex_);
  cand.addUserData<reco::Vertex>("decayVertex", decayVertex);
  cand.addUserData<KinematicParameters>("kinematicParameters", candState.kinematicParameters());
  cand.addUserData<KinematicParametersError>("kinematicParametersError", candState.kinematicParametersError());
  return true;
};


void ParticleFitter::addExtraInfo(pat::GenericParticle& cand) {
  if (!cand.hasUserData("decayVertex") || !cand.hasUserData("primaryVertex")) return;
  const auto& decayVertex = *cand.userData<reco::Vertex>("decayVertex");
  const auto& primaryVertex = *cand.userData<reco::Vertex>("primaryVertex");
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
  const auto pVector3D = SVector3(cand.px(), cand.py(), cand.pz());
  const auto pVector2D = SVector3(cand.px(), cand.py(), 0);;
  const auto sigmaPseudoLvtxMag = std::sqrt(ROOT::Math::Similarity(totalCov, pVector3D)) / cand.p();
  const auto sigmaPseudoRvtxMag = std::sqrt(ROOT::Math::Similarity(totalCov, pVector2D)) / cand.pt();
  // set lifetime infomation
  cand.addUserFloat("lVtxMag", lVtxMag);
  cand.addUserFloat("rVtxMag", rVtxMag);
  cand.addUserFloat("lVtxSig", (lVtxMag/sigmaLvtxMag));
  cand.addUserFloat("rVtxSig", (rVtxMag/sigmaRvtxMag));
  cand.addUserFloat("angle3D", angle3D);
  cand.addUserFloat("angle2D", angle2D);
  cand.addUserFloat("pdlErr3D", sigmaPseudoLvtxMag);
  cand.addUserFloat("pdlErr2D", sigmaPseudoRvtxMag);
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


ParticleDaughter::~ParticleDaughter() {
};

void ParticleDaughter::clear() {
  particles_.clear();
};


void ParticleDaughter::fillInfo(const edm::ParameterSet& pSet, const edm::ParameterSet& config, edm::ConsumesCollector& iC) {
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
    token_dedx_ = iC.consumes<edm::ValueMap<reco::DeDxData>>(config.getParameter<edm::InputTag>("dedxHarmonic2"));
  }
};


template <class T>
void ParticleDaughter::addParticles(const edm::Event& event, const edm::EDGetTokenT<std::vector<T> >& token, const reco::Vertex& vertex, const bool embedInfo) {
  // extract input collections
  edm::Handle<std::vector<T> > handle;
  if (!token.isUninitialized()) event.getByToken(token, handle);
  edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxMap;
  if (!token_dedx_.isUninitialized()) event.getByToken(token_dedx_, dEdxMap);
  edm::Handle<std::vector<float> > mvaColl;
  if (!token_mva_.isUninitialized()) event.getByToken(token_mva_, mvaColl);
  // set selections
  StringCutObjectSelector<T, true> selection(selection_);
  StringCutObjectSelector<pat::GenericParticle, true> finalSelection(finalSelection_);
  // add particles
  if (handle.isValid()) {
    ParticleMassSet particles;
    for (size_t i=0; i<handle->size(); i++) {
      const auto& p = edm::Ref<std::vector<T> >(&(*handle), i);
      if (!selection(*p)) continue;
      if (charge_!=-99 && p->charge()!=charge_) continue;
      const auto p4 = reco::Candidate::PolarLorentzVector(p->pt(), p->eta(), p->phi(), mass_);
      pat::GenericParticle cand;
      cand.setPdgId(pdgId_);
      cand.setP4(p4);
      cand.setCharge(p->charge());
      cand.setVertex(p->vertex());
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
      if (finalSelection(cand)) {
        particles.insert(cand);
      }
    }
    particles_.assign(particles.begin(), particles.end());
  }
};


void ParticleDaughter::addParticles(const edm::Event& event) {
  edm::Handle<pat::GenericParticleCollection> handle;
  event.getByToken(token_source_, handle);
  StringCutObjectSelector<pat::GenericParticle, true> selection(selection_);
  StringCutObjectSelector<pat::GenericParticle, true> finalSelection(finalSelection_);
  if (handle.isValid()) {
    ParticleMassSet particles;
    for (const auto& p : *handle) {
      if (!selection(p)) continue;
      if (charge_!=-99 && p.charge()!=charge_) continue;
      if (p.hasUserData("daughters")) {
        const_cast<pat::GenericParticle*>(&p)->setStatus(2);
      }
      if (finalSelection(p)) {
        particles.insert(p);
      }
    }
    particles_.assign(particles.begin(), particles.end());
  }
};


template <class T>
void ParticleDaughter::addData(pat::GenericParticle& c, const edm::Ref<std::vector<T> >& p, const bool& embedInfo) {
  c.addUserData<T>("src", *p);
};


void ParticleDaughter::addData(pat::GenericParticle& c, const reco::TrackRef& p, const bool& embedInfo) {
  c.setTrack(p, embedInfo);
  c.addUserInt("sourceID", 1);
};


void ParticleDaughter::addData(pat::GenericParticle& c, const reco::PFCandidateRef& p, const bool& embedInfo) {
  if (c.pdgId()==0) { c.setPdgId(p->pdgId()); }
  c.setTrack(p->trackRef(), embedInfo);
  c.addUserData<reco::PFCandidate>("src", *p);
};


void ParticleDaughter::addData(pat::GenericParticle& c, const pat::MuonRef& p, const bool& embedInfo) {
  c.setTrack(p->track(), embedInfo);
  c.setCombinedMuon(p->combinedMuon(), embedInfo);
  c.addUserData<pat::Muon>("src", *p);
};


void ParticleDaughter::addData(pat::GenericParticle& c, const pat::ElectronRef& p, const bool& embedInfo) {
  c.setTrack(p->track(), embedInfo);
  c.setGsfTrack(p->gsfTrack(), embedInfo);
  c.setSuperCluster(p->superCluster(), embedInfo);
  c.setEcalIso(p->ecalIso());
  c.setHcalIso(p->hcalIso());
  c.addUserData<pat::Electron>("src", *p);
};


void ParticleDaughter::setMVA(pat::GenericParticle& cand, const size_t& i, const edm::Handle<std::vector<float> >& mvaColl) {
  if (mvaColl.isValid()) {
    cand.addUserFloat("mva", mvaColl->at(i));
  }
};


void ParticleDaughter::setDeDx(pat::GenericParticle& cand, const edm::Handle<edm::ValueMap<reco::DeDxData> >& dEdxMap) {
  if (cand.track().isNull()) return;
  if (dEdxMap.isValid() && dEdxMap->contains(cand.track().id())) {
    const auto& dEdx = (*dEdxMap)[cand.track()].dEdx();
    cand.addUserFloat("dEdx", dEdx);
  }
};
