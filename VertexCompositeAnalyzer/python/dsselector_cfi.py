import FWCore.ParameterSet.Config as cms

dsselector = cms.EDProducer('VertexCompositeSelector',
  doGenMatching = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(431),
  PID_dau1 = cms.untracked.int32(310),
  PID_dau2 = cms.untracked.int32(321),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDsToKsKCandidatesNew:DSToKsK"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),

  selectGenMatch = cms.untracked.bool(False),
  selectGenUnMatch = cms.untracked.bool(False),
  selectGenMatchSwap = cms.untracked.bool(False),
  selectGenMatchUnSwap = cms.untracked.bool(False),

  useAnyMVA = cms.bool(False),
  useExistingMVA = cms.bool(False),
  mvaType = cms.string('BDT'),
  GBRForestLabel = cms.string('DsInpPb'),
  GBRForestFileName = cms.string(''),
  MVACollection = cms.InputTag("generalDsToKsKCandidatesNew:MVAValues"),
  mvaMax = cms.untracked.double(999.9),
  mvaMin = cms.untracked.double(-999.9),
  mvaCuts = cms.vdouble(-1.,0,0,0,0),

  trkPMin = cms.untracked.double(0.),
  trkPtMin = cms.untracked.double(0.),
  trkEtaMax = cms.untracked.double(999.),
  trkPSumMin = cms.untracked.double(0.),
  trkPtSumMin = cms.untracked.double(0.),
  trkPtAsymMin = cms.untracked.double(0.),
  trkEtaDiffMax = cms.untracked.double(999.),
  trkPtErrMax = cms.untracked.double(999.),
  trkNHitMin = cms.untracked.int32(0),
  candpTMin = cms.untracked.double(-999.),
  candpTMax = cms.untracked.double(999.),
  candYMin = cms.untracked.double(-999.),
  candYMax = cms.untracked.double(999.),
  cand3DDecayLengthSigMin = cms.untracked.double(0.),
  cand3DPointingAngleMax = cms.untracked.double(999.),
  candVtxProbMin = cms.untracked.double(0.)
                              )

dsselectorMC = cms.EDProducer('VertexCompositeSelector',

  doGenMatching = cms.untracked.bool(True),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(True),
  twoLayerDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(431),
  PID_dau1 = cms.untracked.int32(310),
  PID_dau2 = cms.untracked.int32(321),
  deltaR = cms.untracked.double(0.03),

  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDsToKsKCandidatesNew:DSToKsK"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),

  selectGenMatch = cms.untracked.bool(False),
  selectGenUnMatch = cms.untracked.bool(False),
  selectGenMatchSwap = cms.untracked.bool(False),
  selectGenMatchUnSwap = cms.untracked.bool(False),

  useAnyMVA = cms.bool(False),
  useExistingMVA = cms.bool(False),
  mvaType = cms.string('BDT'),
  GBRForestLabel = cms.string('DsInpPb'),
  GBRForestFileName = cms.string(''),
  MVACollection = cms.InputTag("generalDsToKsKCandidatesNew:MVAValues"),
  mvaMax = cms.untracked.double(999.9),
  mvaMin = cms.untracked.double(-999.9),
  mvaCuts = cms.vdouble(-1.,0,0,0,0),

  trkPMin = cms.untracked.double(0.),
  trkPtMin = cms.untracked.double(0.),
  trkEtaMax = cms.untracked.double(999.),
  trkPSumMin = cms.untracked.double(0.),
  trkPtSumMin = cms.untracked.double(0.),
  trkPtAsymMin = cms.untracked.double(0.),
  trkEtaDiffMax = cms.untracked.double(999.),
  trkPtErrMax = cms.untracked.double(999.),
  trkNHitMin = cms.untracked.int32(0),
  candpTMin = cms.untracked.double(-999.),
  candpTMax = cms.untracked.double(999.),
  candYMin = cms.untracked.double(-999.),
  candYMax = cms.untracked.double(999.),
  cand3DDecayLengthSigMin = cms.untracked.double(0.),
  cand3DPointingAngleMax = cms.untracked.double(999.),
  cand3DDCAMin = cms.untracked.double(-999.),
  cand3DDCAMax = cms.untracked.double(999.),
  cand2DDCAMin = cms.untracked.double(-999.),
  cand2DDCAMax = cms.untracked.double(999.),
  candVtxProbMin = cms.untracked.double(0.)
                              )
