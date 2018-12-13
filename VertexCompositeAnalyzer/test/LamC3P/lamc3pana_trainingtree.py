import FWCore.ParameterSet.Config as cms

process = cms.Process("lamc3pana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/Skim_LamC3P_v3/181212_205641/0000/PbPb_6.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_260.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_26.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_259.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_258.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_257.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_256.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_255.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_254.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_253.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_252.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_251.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_250.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_25.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_249.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_248.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_247.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_246.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_245.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_244.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_243.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_242.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_241.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_240.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_24.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_239.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_238.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_237.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_236.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_235.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_234.root',
)
                            )

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity185_*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3pselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3panalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('lamc3pana_training.root')
                                   )

process.lamc3pana_seq = cms.Sequence(process.lamc3pana)
#process.lamc3pana_seq = cms.Sequence(process.hltHM * process.lamc3pana)
#process.lamc3pana_wrongsign_seq = cms.Sequence(process.hltHM * process.lamc3pana_wrongsign)

process.p = cms.Path(process.lamc3pana_seq)
#process.p1 = cms.Path(process.lamc3pana_wrongsign_seq)
