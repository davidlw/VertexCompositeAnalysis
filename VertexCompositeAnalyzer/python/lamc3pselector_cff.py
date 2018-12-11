import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3pselector_cfi import *

lamc3pselectorBDTPreCut = lamc3pselector.clone(
  useAnyMVA = cms.bool(True),

  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkNHitMin = cms.untracked.int32(11),

  cand3DPointingAngleMax = cms.untracked.double(1.0),
  cand2DPointingAngleMax = cms.untracked.double(1.0),
)

lamc3pselectorBDTNonPrompt = lamc3pselectorBDTPreCut.clone(
  useAnyMVA = cms.bool(True),
  BDTCutFileName = cms.string('BDTCuts_NonPrompLamC3P_HM185.root')
#  mvaCuts = cms.vdouble(4.428e-01,-2.639e-04,-2.610e-02,0.0,2.997e-01,-8.948e-02,1.145e-02,-5.321e-04)
#  mvaCuts = cms.vdouble(0.45,-0.0047,-0.023,0.23,-0.087,0.011,-0.0005)
#  mvaCuts = cms.vdouble(0.44,0.0,-0.026,0.3,-0.09,0.011,0.0)
)

lamc3pselectorBDTPrompt = lamc3pselectorBDTPreCut.clone(
  useAnyMVA = cms.bool(True),
  BDTCutFileName = cms.string('BDTCuts_PrompLamC3P_HM185.root')
#  mvaCuts = cms.vdouble(0.4,0.19,-0.088,0.14,-0.0054,-0.0016,0.0001)
#  mvaCuts = cms.vdouble(3.999e-01,1.872e-01,-8.781e-02,1.397e-01,-5.427e-03,-1.607e-03,1.252e-04) v11
#  mvaCuts = cms.vdouble(4.185e-01,-4.441e-02,1.322e-01,-5.939e-02,1.577e-01,-1.177e-02,-6.253e-04,7.098e-05) #v12
)

lamc3pselectorMCBDTPreCut = lamc3pselectorMC.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkNHitMin = cms.untracked.int32(11),
)

lamc3pselectorMCBDTNonPrompt = lamc3pselectorMCBDTPreCut.clone(
  useAnyMVA = cms.bool(True),
  BDTCutFileName = cms.string('BDTCuts_NonPrompLamC3P_HM185.root')
#  mvaCuts = cms.vdouble(4.428e-01,-2.639e-04,-2.610e-02,0.0,2.997e-01,-8.948e-02,1.145e-02,-5.321e-04)
)

lamc3pselectorMCBDTPrompt = lamc3pselectorMCBDTPreCut.clone(
  useAnyMVA = cms.bool(True),
  BDTCutFileName = cms.string('BDTCuts_PrompLamC3P_HM185.root')
#  mvaCuts = cms.vdouble(4.185e-01,-4.441e-02,1.322e-01,-5.939e-02,1.577e-01,-1.177e-02,-6.253e-04,7.098e-05) #v12
)     

lamc3pselectorPID = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),
  userPID = cms.bool(True),

  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorPID2 = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),
  userPID = cms.bool(True),

  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.2),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorCut = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkEtaMax = cms.untracked.double(1.5),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorCutNew = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorCutNew2 = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),

  trkPtSumMin = cms.untracked.double(1.6),
  trkPtMin = cms.untracked.double(0.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)


lamc3pselectorCutNew3 = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),

  trkPtSumMin = cms.untracked.double(1.6),
  trkPtMin = cms.untracked.double(0.7),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorCutNew4 = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),

  trkPtSumMin = cms.untracked.double(1.6),
  trkPtMin = cms.untracked.double(0.8),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorCutNew5 = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),
  
  trkPtSumMin = cms.untracked.double(1.6),
  trkPtMin = cms.untracked.double(0.7),
  trkEtaMax = cms.untracked.double(1.5),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorCutMC = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),
  
  trkPtMin = cms.untracked.double(0.7),
  trkEtaMax = cms.untracked.double(1.5),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorCutNewMC = lamc3pselector.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

lamc3pselectorWS = lamc3pselector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalLamC3PCandidatesNewWrongSign:LamC3P"),
  MVACollection = cms.InputTag("generalLamC3PCandidatesNewWrongSign:MVAValues")
)

lamc3pselectorWSMC = lamc3pselectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalLamC3PCandidatesNewWrongSign:LamC3P"),
  MVACollection = cms.InputTag("generalLamC3PCandidatesNewWrongSign:MVAValues")
)

lamc3pselectorMCGenMatch = lamc3pselectorMC.clone(
  selectGenMatch = cms.untracked.bool(True)
)

lamc3pselectorMCGenUnMatch = lamc3pselectorMC.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)

lamc3pselectorMCGenMatchSwap = lamc3pselectorMC.clone(
  selectGenMatchSwap = cms.untracked.bool(True)
)

lamc3pselectorMCGenMatchUnSwap = lamc3pselectorMC.clone(
  selectGenMatchUnSwap = cms.untracked.bool(True)
)
