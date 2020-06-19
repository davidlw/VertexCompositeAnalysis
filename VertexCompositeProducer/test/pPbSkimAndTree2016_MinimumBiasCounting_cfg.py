import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run2_2016_pA)

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
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PADoubleMuon/AOD/PromptReco-v1/000/286/442/00000/8677067C-53BC-E611-83A0-FA163EDA0C8D.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('80X_dataRun2_v19')

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [ 'HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*' ] # Minimum bias

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping * process.olvFilter_pPb8TeV_dz1p0)

# Luminosity producer
process.lumiInfo = cms.EDProducer('LumiProducerFromBrilcalc',
                                  lumiFile = cms.string("./lumiData.csv"),
                                  throwIfNotFound = cms.bool(False),
                                  doBunchByBunch = cms.bool(False))
process.lumiInfoMB = process.lumiInfo.clone(lumiFile = cms.string("./lumiDataMB.csv"), isTrigger = cms.bool(True))
process.lumiInfoHM = process.lumiInfo.clone(lumiFile = cms.string("./lumiDataHM.csv"), isTrigger = cms.bool(True))
process.lumi_seq = cms.Sequence(process.lumiInfo + process.lumiInfoMB + process.lumiInfoHM)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM * process.lumi_seq )

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.eventAna = particleAna.clone(
  # reconstructed information
  beamSpot = cms.InputTag("offlineBeamSpot"),
  primaryVertices = cms.InputTag("offlinePrimaryVertices"),
  recoParticles = cms.InputTag(""),

  # trigger information
  triggerResults = cms.untracked.InputTag("TriggerResults::HLT"),
  triggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
  triggerInfo = cms.untracked.VPSet([
      #  Double muon triggers
      cms.PSet(path = cms.string('HLT_PAL1DoubleMuOpen_v'), minN = cms.int32(2)), # Dimuons
      # Single muon triggers
      cms.PSet(path = cms.string('HLT_PAL3Mu12_v'), minN = cms.int32(1)), # Electroweak boson
      # Other triggers
      cms.PSet(path = cms.string('HLT_PAFullTracks_Multiplicity120_v')), # High multiplicity
      cms.PSet(path = cms.string('HLT_PAFullTracks_Multiplicity150_v')), # High multiplicity
      cms.PSet(path = cms.string('HLT_PAFullTracks_Multiplicity185_part'), lumiInfo = cms.InputTag("lumiInfoHM", "brilcalc")), # High multiplicity
      cms.PSet(path = cms.string('HLT_PAFullTracks_Multiplicity250_v')), # High multiplicity
      cms.PSet(path = cms.string('HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part'), lumiInfo = cms.InputTag("lumiInfoMB", "brilcalc")), # Minimum bias
  ]),

  #Filter info
  eventFilterResults = cms.untracked.InputTag("TriggerResults"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter',
      'Flag_primaryVertexFilterPA',
      'Flag_NoScraping',
      'Flag_pileupVertexFilterCut',
      'Flag_pileupVertexFilterCutGplus'
  ),
  selectEvents = cms.string("eventFilter_HM_step"),

  # centrality and event plane information
  centralityBin = cms.untracked.InputTag("",""),
  centrality    = cms.untracked.InputTag("pACentrality"),
  eventPlane    = cms.untracked.InputTag(""),

  # luminosity information
  lumiInfo    = cms.untracked.InputTag("lumiInfo", "brilcalc"),
  lumiScalers = cms.untracked.InputTag("scalersRawToDigi"),

  # options
  saveTree  = cms.untracked.bool(True),
  addTrack  = cms.untracked.bool(False),
  addSource = cms.untracked.bool(False),
  addTrgObj = cms.untracked.bool(False),
  dauIDs    = cms.untracked.vint32([])
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('evtana.root'))
process.p = cms.EndPath(process.eventAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.p
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter = cms.Path(process.eventFilter_HM * process.hfCoincFilter)
process.Flag_primaryVertexFilterPA = cms.Path(process.eventFilter_HM * process.primaryVertexFilterPA)
process.Flag_NoScraping = cms.Path(process.eventFilter_HM * process.NoScraping)
process.Flag_pileupVertexFilterCut = cms.Path(process.eventFilter_HM * process.olvFilter_pPb8TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus = cms.Path(process.eventFilter_HM * process.pileUpFilter_pPb8TeV_Gplus)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter , process.Flag_primaryVertexFilterPA , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
