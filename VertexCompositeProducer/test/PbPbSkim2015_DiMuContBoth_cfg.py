import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'/store/hidata/HIRun2015/HIOniaL1DoubleMu0/AOD/PromptReco-v1/000/263/757/00000/F2BD4D30-ACBB-E511-9FC2-02163E012647.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '75X_dataRun2_v12'

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltMuon = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMuon.HLTPaths = ['HLT_HIL1DoubleMu0*_v*']
process.hltMuon.andOr = cms.bool(True)
process.hltMuon.throw = cms.bool(False)

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.evtsel_filter = cms.Sequence(process.hltMuon * process.hfCoincFilter3 * process.primaryVertexFilter)

process.eventFilter_step = cms.Path( process.evtsel_filter )

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")

process.generalMuMuContinuimOneStTightCandidatesWrongSign = process.generalMuMuContinuimOneStTightCandidates.clone(isWrongSign = cms.bool(True))
process.generalMuMuContinuimOneStTightPFCandidatesWrongSign = process.generalMuMuContinuimOneStTightPFCandidates.clone(isWrongSign = cms.bool(True))
process.generalMuMuContinuimPFCandidatesWrongSign = process.generalMuMuContinuimPFCandidates.clone(isWrongSign = cms.bool(True))
process.generalMuMuContinuimOneStTightGlobalCandidatesWrongSign = process.generalMuMuContinuimOneStTightGlobalCandidates.clone(isWrongSign = cms.bool(True))
process.generalMuMuContinuimOneStTightPFGlobalCandidatesWrongSign = process.generalMuMuContinuimOneStTightPFGlobalCandidates.clone(isWrongSign = cms.bool(True))
process.generalMuMuContinuimPFGlobalCandidatesWrongSign = process.generalMuMuContinuimPFGlobalCandidates.clone(isWrongSign = cms.bool(True))

#process.dimurereco_step = cms.Path( process.evtsel_filter * process.generalMuMuContinuimOneStTightCandidates * process.generalMuMuContinuimOneStTightPFCandidates * process.generalMuMuContinuimPFCandidates * process.generalMuMuContinuimOneStTightGlobalCandidates * process.generalMuMuContinuimOneStTightPFGlobalCandidates * process.generalMuMuContinuimPFGlobalCandidates)
#process.dimurerecowrongsign_step = cms.Path( process.evtsel_filter * process.generalMuMuContinuimOneStTightCandidatesWrongSign * process.generalMuMuContinuimOneStTightPFCandidatesWrongSign * process.generalMuMuContinuimPFCandidatesWrongSign * process.generalMuMuContinuimOneStTightGlobalCandidates * process.generalMuMuContinuimOneStTightPFGlobalCandidates * process.generalMuMuContinuimPFGlobalCandidates)

process.dimurereco_step = cms.Path( process.evtsel_filter * process.generalMuMuContinuimOneStTightGlobalCandidates)
process.dimurerecowrongsign_step = cms.Path( process.evtsel_filter * process.generalMuMuContinuimOneStTightGlobalCandidatesWrongSign)

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentJPsi_cff")
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('PbPb_DiMuCont.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.eventFilter_step,
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.output_step
)
