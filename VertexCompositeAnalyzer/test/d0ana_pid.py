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
process.d0ana.saveHistogram = cms.untracked.bool(True)
process.d0ana.saveTree = cms.untracked.bool(False)

process.d0ana_pos = process.d0ana.clone()
process.d0ana_neg = process.d0ana.clone()
process.d0ana_pos.VertexCompositeCollection = cms.untracked.InputTag("d0selectorPOS:D0")
process.d0ana_pos.MVACollection = cms.InputTag("d0selectorPOS:MVAValuesNewD0")
process.d0ana_neg.VertexCompositeCollection = cms.untracked.InputTag("d0selectorNEG:D0")
process.d0ana_neg.MVACollection = cms.InputTag("d0selectorNEG:MVAValuesNewD0")

process.d0selectorPID2.cand3DDecayLengthSigMin = cms.untracked.double(-10000.)
process.d0selectorPID2.cand3DPointingAngleMax = cms.untracked.double(10000.)
process.d0selectorPID2.candVtxProbMin = cms.untracked.double(-10000.) 

process.d0preselectorPOS = process.d0selectorPID2.clone( selectFlavor = cms.untracked.int32(1) )
process.d0preselectorNEG = process.d0selectorPID2.clone( selectFlavor = cms.untracked.int32(-1) )

process.d0selector.useAnyMVA = cms.bool(True)
process.d0selectorPOS = process.d0selector.clone()
process.d0selectorNEG = process.d0selector.clone()
process.d0selectorPOS.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_pT15_24_pid_flavor1.root')
process.d0selectorPOS.VertexCompositeCollection = cms.untracked.InputTag("d0preselectorPOS:D0")
process.d0selectorNEG.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_pT15_24_pid_flavor-1.root')
process.d0selectorNEG.VertexCompositeCollection = cms.untracked.InputTag("d0preselectorNEG:D0")

process.d0ana_pos_seq = cms.Sequence(process.d0preselectorPOS * process.d0selectorPOS * process.d0ana_pos)
process.d0ana_neg_seq = cms.Sequence(process.d0preselectorNEG * process.d0selectorNEG * process.d0ana_neg)

process.p = cms.Path(process.d0ana_pos_seq)
process.p1 = cms.Path(process.d0ana_neg_seq)
