import FWCore.ParameterSet.Config as cms

process = cms.Process("jpsiana")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(15000)
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/pPb_Skim_JPsiBoth_default_v1/180827_220926/0000/pPb_HM_JPsi_27.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/480/00000/4C0D189A-1BAF-E611-B5CA-FA163EAF1F45.root',
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/6A19A7DC-92AF-E611-AB81-02163E0142DF.root',
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/84DA1928-90AF-E611-A838-02163E01240B.root'
)
                            )

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity185_*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('jpsiana.root')
                                   )

process.jpsiana.useAnyMVA = cms.bool(True)
process.jpsiana_wrongsign.useAnyMVA = cms.bool(True)
process.jpsiana.VertexCompositeCollection = cms.untracked.InputTag("jpsiselector:DiMu")
process.jpsiana_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("jpsiselectorWS:DiMu")
process.jpsiana.MVACollection = cms.InputTag("jpsiselector:MVAValuesNewDiMu")
process.jpsiana_wrongsign.MVACollection = cms.InputTag("jpsiselectorWS:MVAValuesNewDiMu")
#process.jpsiana.isSkimMVA = cms.untracked.bool(True)
#process.jpsiana_wrongsign.isSkimMVA = cms.untracked.bool(True)
process.jpsiana.saveHistogram = cms.untracked.bool(True)
process.jpsiana.saveTree = cms.untracked.bool(False)
process.jpsiana_wrongsign.saveHistogram = cms.untracked.bool(True)
process.jpsiana_wrongsign.saveTree = cms.untracked.bool(False)

process.jpsiana1 = process.jpsiana.clone()
process.jpsiana2 = process.jpsiana.clone()
process.jpsiana1.VertexCompositeCollection = cms.untracked.InputTag("jpsiselector1:DiMu")
process.jpsiana1.MVACollection = cms.InputTag("jpsiselector1:MVAValuesNewDiMu")
process.jpsiana2.VertexCompositeCollection = cms.untracked.InputTag("jpsiselector2:DiMu")
process.jpsiana2.MVACollection = cms.InputTag("jpsiselector2:MVAValuesNewDiMu")

#process.jpsiselectorCutNew.cand3DDecayLengthSigMin = cms.untracked.double(-10000.)
#process.jpsiselectorCutNew.cand3DPointingAngleMax = cms.untracked.double(10000.)
#process.jpsiselectorCutNew.candVtxProbMin = cms.untracked.double(-10000.) 
#process.jpsipreselector = process.jpsiselectorCutNew.clone()
#process.jpsipreselectorWS = process.jpsiselectorCutNew.clone()

process.jpsiselector.useAnyMVA = cms.bool(True)
process.jpsiselectorWS.useAnyMVA = cms.bool(True)
process.jpsiselector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptJPsiInpPb_scenario1_jpsi.root')
process.jpsiselectorWS.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptJPsiInpPb_scenario1_jpsi.root')
#process.jpsiselector2.VertexCompositeCollection = cms.untracked.InputTag("jpsipreselector:D0")
#process.jpsiselector2WS.VertexCompositeCollection = cms.untracked.InputTag("jpsipreselectorWS:D0")
process.jpsiselector1.useAnyMVA = cms.bool(True)
process.jpsiselector1.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptJPsiInpPb_scenario1_jpsi1.root')
process.jpsiselector2.useAnyMVA = cms.bool(True)
process.jpsiselector2.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptJPsiInpPb_scenario1_jpsi2.root')

process.jpsiselector.multMin = cms.untracked.double(185)
process.jpsiselector.multMax = cms.untracked.double(250)
process.jpsiselectorWS.multMin = cms.untracked.double(185)
process.jpsiselectorWS.multMax = cms.untracked.double(250)
process.jpsiselector1.multMin = cms.untracked.double(185)
process.jpsiselector1.multMax = cms.untracked.double(250)
process.jpsiselector2.multMin = cms.untracked.double(185)
process.jpsiselector2.multMax = cms.untracked.double(250)



#process.jpsiana_seq = cms.Sequence(process.jpsipreselector * process.jpsiselector * process.jpsiana)
#process.jpsiana_wrongsign_seq = cms.Sequence(process.jpsipreselectorWS * process.jpsiselectorWS * process.jpsiana_wrongsign)
process.jpsiana_seq = cms.Sequence(process.hltHM * process.jpsiselector * process.jpsiana)
process.jpsiana_wrongsign_seq = cms.Sequence(process.hltHM * process.jpsiselector2WS * process.jpsiana_wrongsign)
process.jpsiana1_seq = cms.Sequence(process.hltHM * process.jpsiselector1 * process.jpsiana1)
process.jpsiana2_seq = cms.Sequence(process.hltHM * process.jpsiselector2 * process.jpsiana2)


process.p = cms.Path(process.jpsiana_seq)
process.p1 = cms.Path(process.jpsiana1_seq)
process.p2 = cms.Path(process.jpsiana2_seq)
