import FWCore.ParameterSet.Config as cms
process = cms.Process("writeGBRForests")
# CV: needs to be set to 1 so that GBRForestWriter::analyze method gets called exactly once
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")
process.load('Configuration/StandardSequences/Services_cff')

process.gbrForestWriter = cms.EDAnalyzer("GBRForestWriter",
    jobs = cms.VPSet(
        cms.PSet(
            inputFileName = cms.string('weights/TMVAClassification_BDT.weights.xml'),
            inputFileType = cms.string("XML"),
            inputVariables = cms.vstring(['pT','y','VtxProb','3DDecayLengthSignificance','2DDecayLengthSignificance','3DDecayLength','3DPointingAngle','2DPointingAngle','zDCASignificanceDaugther1','zDCASignificanceDaugther2','xyDCASignificanceDaugther1','xyDCASignificanceDaugther2','pTD1','pTD2','EtaD1','EtaD2','NHitD1','NHitD2','pTerrD1','pTerrD2']),
#            inputVariables = cms.vstring(['pT','y','VtxProb','3DDecayLengthSignificance','2DDecayLengthSignificance','3DDecayLength','3DPointingAngle','2DPointingAngle','zDCASignificanceDaugther1','zDCASignificanceDaugther2','xyDCASignificanceDaugther1','xyDCASignificanceDaugther2','pTD1','pTD2','EtaD1','EtaD2','NHitD1','NHitD2','pTerrD1','pTerrD2','bestvtxZ']),
#            inputVariables = cms.vstring(['pT','y','VtxProb','3DDecayLengthSignificance','2DDecayLengthSignificance','3DDecayLength','3DPointingAngle','2DPointingAngle','zDCASignificanceDaugther1','zDCASignificanceDaugther2','xyDCASignificanceDaugther1','xyDCASignificanceDaugther2','pTD1','pTD2','EtaD1','EtaD2']),
            spectatorVariables = cms.vstring(),
            methodName = cms.string("BDT"),
            gbrForestName = cms.string("D0InpPb"),
            outputFileType = cms.string("GBRForest"),
            outputFileName = cms.string("GBRForestfile_BDT_PromptD0InpPb_default_HLT185_v1.root")
        )
    )
)
process.p = cms.Path(process.gbrForestWriter)
