import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalV0Candidates_cfi import *

generalV0CandidatesNew = generalV0Candidates.clone (
    selectD0s = cms.bool(False),
    selectLambdaCs = cms.bool(False),
    selectXis = cms.bool(False),
    selectOmegas = cms.bool(False),

    tkNhitsCut = cms.int32(3),
    tkChi2Cut = cms.double(7.0),
    dauTransImpactSigCut = cms.double(1.0),
    dauLongImpactSigCut = cms.double(1.0),
    xiVtxSignificance3DCut = cms.double(0.0),
    xiVtxSignificance2DCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(2.5),

    collinearityCut = cms.double(0.997),

    innerHitPosCut = cms.double(-1)
)

generalCascadeCandidatesNew = generalV0Candidates.clone(
    selectD0s = cms.bool(False),
    selectLambdaCs = cms.bool(False),
    selectXis = cms.bool(True),
    selectOmegas = cms.bool(True),

    tkNhitsCut = cms.int32(3),
    tkChi2Cut = cms.double(7.0),
    dauTransImpactSigCut = cms.double(2.5),
    dauLongImpactSigCut = cms.double(2.5),
    batDauTransImpactSigCut = cms.double(4.0),
    batDauLongImpactSigCut = cms.double(4.0),
    xiVtxSignificance3DCut = cms.double(2.0),
    xiVtxSignificance2DCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(8),

    collinearityCut = cms.double(0.5),
    xiCollinearityCut = cms.double(-2.0),

    innerHitPosCut = cms.double(-1)
)

generalDsToKsKCandidatesNew = generalV0Candidates.clone (
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(False),
    selectDSToKsKs = cms.bool(True),
    selectDSToPhiPis = cms.bool(False),
    selectXis = cms.bool(False),
    selectOmegas = cms.bool(False),

    tkNhitsCut = cms.int32(3),
    tkChi2Cut = cms.double(7.0),
    dauTransImpactSigCut = cms.double(1.0),
    dauLongImpactSigCut = cms.double(1.0),
    batDauTransImpactSigCut = cms.double(0.0),
    batDauLongImpactSigCut = cms.double(0.0),
    xiVtxSignificance3DCut = cms.double(0.0),
    xiVtxSignificance2DCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(4.5),

    collinearityCut = cms.double(0.995),
    xiCollinearityCut = cms.double(-2.0),

    innerHitPosCut = cms.double(-1)
)

generalDsToPhiPiCandidatesNew = generalV0Candidates.clone (
    selectKshorts = cms.bool(False),
    selectLambdas = cms.bool(False),
    selectPhis = cms.bool(True),
    selectDSToKsKs = cms.bool(False),
    selectDSToPhiPis = cms.bool(True),
    selectXis = cms.bool(False),
    selectOmegas = cms.bool(False),

#    doVertexFit = cms.bool(False),
#    vertexFitter = cms.InputTag('None'),

    mKKCutMin = cms.double(1.0),
    mKKCutMax = cms.double(1.04),

    trackQualities = cms.vstring('highPurity'),
    collinearityCut = cms.double(-2.0),
    dauTransImpactSigCut = cms.double(-1.),
    dauLongImpactSigCut = cms.double(-1.),
    vtxSignificance3DCut = cms.double(-999.),

    batDauTransImpactSigCut = cms.double(-1.),
    batDauLongImpactSigCut = cms.double(-1.),
    xiVtxSignificance3DCut = cms.double(-999.),
    xiVtxSignificance2DCut = cms.double(-999.),
    xiCollinearityCut = cms.double(-2.0),

    innerHitPosCut = cms.double(-1)
)

generalDpmCandidatesNew = generalV0Candidates.clone (
    selectDPMs = cms.bool(True),
    selectLambdas = cms.bool(False),
    selectXis = cms.bool(False),
    selectOmegas = cms.bool(False),

    tkNhitsCut = cms.int32(3),
    tkChi2Cut = cms.double(7.0),
    dauTransImpactSigCut = cms.double(1.0),
    dauLongImpactSigCut = cms.double(1.0),
    batDauTransImpactSigCut = cms.double(0.0),
    batDauLongImpactSigCut = cms.double(0.0),
    xiVtxSignificance3DCut = cms.double(0.0),
    xiVtxSignificance2DCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(4.5),

    collinearityCut = cms.double(0.995),
    xiCollinearityCut = cms.double(-2.0),

    innerHitPosCut = cms.double(-1)
)

generalLambdaCCandidatesNew = generalV0Candidates.clone (
    selectXis = cms.bool(False),
    selectOmegas = cms.bool(False),
    selectLambdaCToLamPis = cms.bool(True),
    selectLambdaCToKsPs = cms.bool(True),

    tkNhitsCut = cms.int32(3),
    tkChi2Cut = cms.double(7.0),
    dauTransImpactSigCut = cms.double(1.0),
    dauLongImpactSigCut = cms.double(1.0),
    batDauTransImpactSigCut = cms.double(0.0),
    batDauLongImpactSigCut = cms.double(0.0),
    xiVtxSignificance3DCut = cms.double(0.0),
    xiVtxSignificance2DCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(4.5),

    collinearityCut = cms.double(0.995),
    xiCollinearityCut = cms.double(-2.0),

    innerHitPosCut = cms.double(-1)
)

generalPhiCandidatesNew = generalV0Candidates.clone(

    selectKshorts = cms.bool(False),
    selectLambdas = cms.bool(False),
    selectPhis = cms.bool(True),
    selectXis = cms.bool(False),
    selectOmegas = cms.bool(False),

    doVertexFit = cms.bool(False),
    vertexFitter = cms.InputTag('None'),

    mKKCutMin = cms.double(0.98),
    mKKCutMax = cms.double(1.06),

    trackQualities = cms.vstring('highPurity'),
    collinearityCut = cms.double(-2.0),

    # Track impact parameter significance 
    dauTransImpactSigCut = cms.double(-1.),
    dauLongImpactSigCut = cms.double(-1.),

    vtxSignificance3DCut = cms.double(-999.)
)
