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
'root://cms-xrd-global.cern.ch//store/user/davidlw/HighMultiplicityEOF/pp_Skim_D0Both_default_v1/180920_103725/0000/pPb_HM_9.root'
#'file:/afs/cern.ch/user/d/davidlw/CMSSW/CMSSW_8_0_28_jpsi/src/VertexCompositeAnalysis/VertexCompositeProducer/test/pPb_HM.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/data/Run2018C/HighMultiplicityEOF/AOD/PromptReco-v2/000/319/463/00000/94837F54-8086-E811-AC90-02163E017EE5.root',
'root://cms-xrd-global.cern.ch//store/data/Run2018C/HighMultiplicityEOF/AOD/PromptReco-v2/000/319/463/00000/02C908C8-7F86-E811-996F-FA163E254550.root',
'root://cms-xrd-global.cern.ch//store/data/Run2018C/HighMultiplicityEOF/AOD/PromptReco-v2/000/319/463/00000/02A952B5-8686-E811-BCA0-02163E010ECC.root',
'root://cms-xrd-global.cern.ch//store/data/Run2018C/HighMultiplicityEOF/AOD/PromptReco-v2/000/319/462/00000/6AAFDF92-8086-E811-BD35-FA163E9925A4.root',
),
          dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
          inputCommands=cms.untracked.vstring(
                  'keep *',
                  'drop *_*totem*_*_*'
          )
                            )

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity185_*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana.root')
                                   )

process.d0ana.useAnyMVA = cms.bool(True)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
#process.d0ana.isSkimMVA = cms.untracked.bool(True)
process.d0ana.saveHistogram = cms.untracked.bool(True)
process.d0ana.saveTree = cms.untracked.bool(False)
process.d0ana.yBins = cms.untracked.vdouble(-2.4,-1.6,-0.8,0.0,0.8,1.6,2.4)

process.d0selector = process.d0selectorBDTPreCut.clone()
process.d0selector.useAnyMVA = cms.bool(True)
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0Inpp_default_N100_WS_v2.root')
process.d0selector.GBRForestLabel = cms.string('D0Inpp')
process.d0selector.multMin = cms.untracked.double(100)
process.d0selector.multMax = cms.untracked.double(1000)

process.npd0selector = process.d0selector.clone()
process.npd0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0Inpp_default_N100_WS_v2.root')
process.npd0ana = process.d0ana.clone()
process.npd0ana.VertexCompositeCollection = cms.untracked.InputTag("npd0selector:D0")
process.npd0ana.MVACollection = cms.InputTag("npd0selector:MVAValuesNewD0")

#process.d0ana_seq = cms.Sequence(process.hltHM * process.d0selector * process.d0ana)
#process.npd0ana_seq = cms.Sequence(process.hltHM * process.npd0selector * process.npd0ana)
process.d0ana_seq = cms.Sequence(process.d0selector * process.d0ana)
process.npd0ana_seq = cms.Sequence(process.npd0selector * process.npd0ana)

process.p = cms.Path(process.d0ana_seq)
process.p1 = cms.Path(process.npd0ana_seq)
