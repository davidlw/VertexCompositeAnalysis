import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov//store/data/Run2017E/HighMultiplicityEOF5/AOD/17Nov2017-v1/50000/28DB5048-B1EB-E711-A3B2-0CC47AC08BF8.root'
),
          dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
          inputCommands=cms.untracked.vstring(
                  'keep *',
                  'drop *_*totem*_*_*'
          )
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '94X_dataRun2_v11'

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

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNew.trkPtSumCut = cms.double(1.6)
process.generalD0CandidatesNew.trkEtaDiffCut = cms.double(1.0)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.7)
process.generalD0CandidatesNew.alphaCut = cms.double(1.0)
process.generalD0CandidatesNew.alpha2DCut = cms.double(1.0)

#process.generalD0CandidatesNew.useAnyMVA = cms.bool(True)
#process.generalD0CandidatesNew.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_scenario2.root')
process.generalD0CandidatesNewWrongSign = process.generalD0CandidatesNew.clone(isWrongSign = cms.bool(True))

process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew * process.generalD0CandidatesNewWrongSign )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('pp_HM.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_HM_step = cms.EndPath(process.output_HM)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
    process.output_HM_step
)
