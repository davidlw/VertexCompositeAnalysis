import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.dsselector_cfi import *

dsselectorMCGenMatch = dsselectorMC.clone(
  selectGenMatch = cms.untracked.bool(True)
)

dsselectorMCGenUnMatch = dsselectorMC.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)

dsselectorBDTPreCut = dsselector.clone(
  useAnyMVA = cms.bool(False),

#  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.2),
  trkEtaDiffMax = cms.untracked.double(1.),
#  trkNHitMin = cms.untracked.int32(11),
)

