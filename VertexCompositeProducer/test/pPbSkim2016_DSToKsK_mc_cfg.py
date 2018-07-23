import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_1.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_2.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_3.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_4.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_5.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_6.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_7.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_8.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_9.root',
'/store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_10.root',
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '80X_mcRun2_pA_v4'

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity280_v*']
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity280_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
#    cut = cms.string("!isFake && abs(z) <= 1 && position.Rho <= 2 && tracksSize >= 5"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

#Reject beam scraping events standard pp configuration
process.NoScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.PAcollisionEventSelection = cms.Sequence(
#                                         process.hfCoincFilter * 
                                         process.PAprimaryVertexFilter *
                                         process.NoScraping
                                         )

process.eventFilter_HM = cms.Sequence( 
#    process.hltHM *
    process.PAcollisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#process.dEdx_step = cms.Path( process.eventFilter_HM * process.produceEnergyLoss )

########## Ds candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalV0Candidates_cff")
process.dsrereco_step = cms.Path( process.eventFilter_HM * process.generalDsToKsKCandidatesNew)

###############################################################################################

#process.load("RiceHIG.Skim2013.ppanalysisSkimContentFull_cff")
#process.load("RiceHIG.Skim2013.ppanalysisSkimContentSlim_cff")
process.load("RiceHIG.Skim2013.ppanalysisSkimContentV0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('pPb_HM_DS.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_HM_step = cms.EndPath(process.output_HM)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.dsrereco_step,
    process.output_HM_step
)
