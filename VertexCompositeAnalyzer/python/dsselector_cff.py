import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.dsselector_cfi import *

dsselectorMCGenMatch = dsselectorMC.clone(
  selectGenMatch = cms.untracked.bool(True)
)

dsselectorMCGenUnMatch = dsselectorMC.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)
