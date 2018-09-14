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
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_SB.root')
process.d0selector.multMin = cms.untracked.double(185)
process.d0selector.multMax = cms.untracked.double(250)

process.d0ana_seq = cms.Sequence(process.hltHM * process.d0selector * process.d0ana)

process.p = cms.Path(process.d0ana_seq)
