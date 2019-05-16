import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run2_2018)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/data/Run2017C/HighMultiplicityEOF1/AOD/12Sep2017-v1/150000/CED96518-90B0-E711-BA06-FA163E51A507.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v10', '')

# Add the VertexComposite producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalMuMuMassMin2CandidatesWrongSign = process.generalMuMuMassMin2Candidates.clone(isWrongSign = cms.bool(True))
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, False)

# Add centrality producer
process.load("RecoHI.HiCentralityAlgos.pACentrality_cfi")
process.pACentrality.srcHFhits = cms.InputTag("reducedHcalRecHits:hfreco")
process.pACentrality.srcEBhits = cms.InputTag("reducedEcalRecHitsEB")
process.pACentrality.srcEEhits = cms.InputTag("reducedEcalRecHitsEE")
process.pACentrality.produceHFtowers = cms.bool(False)
process.pACentrality.produceETmidRapidity = cms.bool(False)
process.pACentrality.producePixelhits = cms.bool(False)
process.pACentrality.producePixelTracks = cms.bool(False)
process.cent_seq = cms.Sequence(process.pACentrality)

# Add muon event selection
process.twoMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("muons"), minNumber = cms.uint32(2))
process.goodMuon = cms.EDFilter("MuonSelector",
            src = cms.InputTag("muons"),
            cut = process.generalMuMuMassMin2Candidates.muonSelection,
            )
process.twoGoodMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodMuon"), minNumber = cms.uint32(2))
process.goodDimuon = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = process.generalMuMuMassMin2Candidates.candidateSelection,
            checkCharge = cms.bool(False),
            decay = cms.string('goodMuon@+ goodMuon@-')
            )
process.oneGoodDimuon = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDimuon"), minNumber = cms.uint32(1))
process.dimuonEvtSel = cms.Sequence(process.twoMuons * process.goodMuon * process.twoGoodMuons * process.goodDimuon * process.oneGoodDimuon)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # Other triggers
    'HLT_FullTrack_Multiplicity85_*', # High multiplicity
    'HLT_FullTrack_Multiplicity100_*', # High multiplicity
    'HLT_FullTrack_Multiplicity130_*', # High multiplicity
    'HLT_FullTrack_Multiplicity155_*', # High multiplicity
    'HLT_L1MinimumBiasHF_OR_*', # Minimum bias
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter * process.NoScraping)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.dimuonEvtSel
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.pcentandep_step = cms.Path(process.eventFilter_HM * process.cent_seq)
process.dimurereco_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalMuMuMassMin2Candidates)
process.dimurerecowrongsign_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalMuMuMassMin2CandidatesWrongSign)

# Add the VertexComposite tree
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")
process.dimucontana.selectEvents = cms.untracked.string("eventFilter_HM_step")
process.dimucontana_wrongsign.selectEvents = cms.untracked.string("eventFilter_HM_step")

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuana.root'))
process.p = cms.EndPath(process.dimucontana * process.dimucontana_wrongsign)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.pcentandep_step,
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.p
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_NoScraping = cms.Path(process.eventFilter_HM * process.NoScraping)
process.Flag_pileupVertexFilterCut = cms.Path(process.eventFilter_HM * process.olvFilter_pp5TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus = cms.Path(process.eventFilter_HM * process.pileUpFilter_pp5TeV_Gplus)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_primaryVertexFilter , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
