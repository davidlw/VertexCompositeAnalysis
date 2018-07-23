import FWCore.ParameterSet.Config as cms

process = cms.Process("d0ana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/pPb_Skim_D0Both_v7_test/180704_132822/0000/pPb_HM_976.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/98B099A2-2DB0-E611-BE2C-02163E0119D7.root'
)
                            )
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana.root')
                                   )

process.d0ana.useAnyMVA = cms.bool(True)
process.d0ana_wrongsign.useAnyMVA = cms.bool(True)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWS:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
process.d0ana_wrongsign.MVACollection = cms.InputTag("d0selectorWS:MVAValuesNewD0")
#process.d0ana.isSkimMVA = cms.untracked.bool(True)
#process.d0ana_wrongsign.isSkimMVA = cms.untracked.bool(True)
process.d0ana.saveHistogram = cms.untracked.bool(True)
process.d0ana.saveTree = cms.untracked.bool(False)
process.d0ana_wrongsign.saveHistogram = cms.untracked.bool(True)
process.d0ana_wrongsign.saveTree = cms.untracked.bool(False)

process.d0selectorCutNew3.cand3DDecayLengthSigMin = cms.untracked.double(-10000.)
process.d0selectorCutNew3.cand3DPointingAngleMax = cms.untracked.double(10000.)
process.d0selectorCutNew3.candVtxProbMin = cms.untracked.double(-10000.) 

process.d0preselector = process.d0selectorCutNew3.clone()
process.d0preselectorWS = process.d0selectorCutNew3.clone()

#process.d0preselector = process.d0selector.clone(
#  trkPtErrMax = cms.untracked.double(0.1),
#  trkNHitMin = cms.untracked.int32(11)
#)
#process.d0preselectorWS = process.d0selectorWS.clone(
#  trkPtErrMax = cms.untracked.double(0.1),
#  trkNHitMin = cms.untracked.int32(11)
#)

process.d0selector.useAnyMVA = cms.bool(True)
process.d0selectorWS.useAnyMVA = cms.bool(True)
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_pT15_24_scenario2.root')
process.d0selectorWS.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_pT15_24_scenario2.root')
process.d0selector.VertexCompositeCollection = cms.untracked.InputTag("d0preselector:D0")
process.d0selectorWS.VertexCompositeCollection = cms.untracked.InputTag("d0preselectorWS:D0")

process.d0ana_seq = cms.Sequence(process.d0preselector * process.d0selector * process.d0ana)
process.d0ana_wrongsign_seq = cms.Sequence(process.d0preselectorWS * process.d0selectorWS * process.d0ana_wrongsign)

process.p = cms.Path(process.d0ana_seq)
#process.p1 = cms.Path(process.d0ana_wrongsign_seq)
