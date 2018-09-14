import FWCore.ParameterSet.Config as cms
process = cms.Process("writeGBRForests")
# CV: needs to be set to 1 so that GBRForestWriter::analyze method gets called exactly once
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")
process.load('Configuration/StandardSequences/Services_cff')

process.gbrForestWriter = cms.EDAnalyzer("GBRForestWriter",
    jobs = cms.VPSet(
        cms.PSet(
            inputFileName = cms.string('weights_prompt_pPb_default2_jpsi/TMVAClassification_BDT.weights.xml'),
            inputFileType = cms.string("XML"),
            inputVariables = cms.vstring(['pT','y','VtxProb','3DDecayLengthSignificance','2DDecayLengthSignificance','3DDecayLength','zDCASignificanceDaugther1','zDCASignificanceDaugther2','xyDCASignificanceDaugther1','xyDCASignificanceDaugther2','NHitD1','NHitD2','nMatchedChamberD1','nMatchedStationD1','EnergyDepositionD1','nMatchedChamberD2','nMatchedStationD2','EnergyDepositionD2','dxSig1_seg','dySig1_seg','ddxdzSig1_seg','ddydzSig1_seg','dxSig2_seg','dySig2_seg','ddxdzSig2_seg','ddydzSig2_seg','pTD1','pTD2','EtaD1','EtaD2']),
            spectatorVariables = cms.vstring(),
            methodName = cms.string("BDT"),
            gbrForestName = cms.string("JPsiInpPb"),
            outputFileType = cms.string("GBRForest"),
            outputFileName = cms.string("GBRForestfile_BDT_PromptJPsiInpPb_pPb_default2_jpsi.root")
        )
    )
)
process.p = cms.Path(process.gbrForestWriter)
