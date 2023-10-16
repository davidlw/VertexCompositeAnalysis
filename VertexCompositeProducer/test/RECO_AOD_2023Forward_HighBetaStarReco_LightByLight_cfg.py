# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: reco --conditions 132X_dataRun3_Prompt_v4 -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent AOD --data --process RECO --scenario pp --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run3 --no_exec --era Run3_2023 --repacked --nThread 8
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
from Configuration.Eras.Modifier_highBetaStar_2018_cff import highBetaStar_2018
from Configuration.ProcessModifiers.egamma_lowPt_exclusive_cff import egamma_lowPt_exclusive
Run3_2023_UPC = cms.ModifierChain(Run3_2023, highBetaStar_2018, egamma_lowPt_exclusive)

process = cms.Process('RECO',Run3_2023_UPC)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_DataMapper_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/hidata/HIRun2023A/HIForward0/RAW/v1/000/374/833/00000/eab8fb72-d975-4085-a0d3-857e5671a4da.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('reco_RAW2DIGI_L1Reco_RECO.root'),
    outputCommands = process.AODEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v4', '')

# Overload the global tag
tags = {"electron_eb_ecalOnly_1To20_0p2To2_mean": "electron_eb_ecalOnly_1To20_0p2To2_mean_LbyL2018V2",
        "pfscecal_eeUncertainty_offline_v2": "pfscecal_eeUncertainty_offline_v2_LbyL2018V2",
        "pfscecal_ebUncertainty_offline_v2": "pfscecal_ebUncertainty_offline_v2_LbyL2018V2",
        "electron_eb_ecalTrk_1To20_0p0002To0p5_sigma": "electron_eb_ecalTrk_1To20_0p0002To0p5_sigma_LbyL2018V2",
        "photon_ee_ecalOnly_1To20_0p2To2_mean": "photon_ee_ecalOnly_1To20_0p2To2_mean_LbyL2018V2",
        "electron_ee_ecalOnly_1To20_0p2To2_mean": "electron_ee_ecalOnly_1To20_0p2To2_mean_LbyL2018V2",
        "electron_ee_ecalTrk_1To20_0p0002To0p5_sigma": "electron_ee_ecalTrk_1To20_0p0002To0p5_sigma_LbyL2018V2",
        "photon_ee_ecalOnly_1To20_0p0002To0p5_sigma": "photon_ee_ecalOnly_1To20_0p0002To0p5_sigma_LbyL2018V2",
        "electron_ee_ecalTrk_1To20_0p2To2_mean": "electron_ee_ecalTrk_1To20_0p2To2_mean_LbyL2018V2",
        "electron_eb_ecalTrk_1To20_0p2To2_mean": "electron_eb_ecalTrk_1To20_0p2To2_mean_LbyL2018V2",
        "electron_eb_ecalOnly_1To20_0p2To2_mean": "electron_eb_ecalOnly_1To20_0p2To2_mean_LbyL2018V2",
        "electron_ee_ecalOnly_1To20_0p0002To0p5_sigma": "electron_ee_ecalOnly_1To20_0p0002To0p5_sigma_LbyL2018V2",
        "photon_eb_ecalOnly_1To20_0p2To2_mean": "photon_eb_ecalOnly_1To20_0p2To2_mean_LbyL2018V2",
        "photon_eb_ecalOnly_1To20_0p0002To0p5_sigma": "photon_eb_ecalOnly_1To20_0p0002To0p5_sigma_LbyL2018V2",
        "pfscecal_ebCorrection_offline_v2": "pfscecal_ebCorrection_offline_v2_LbyL2018V2",
        "electron_eb_ecalOnly_1To20_0p0002To0p5_sigma": "electron_eb_ecalOnly_1To20_0p0002To0p5_sigma_LbyL2018V2",
        "pfscecal_eeCorrection_offline_v2": "pfscecal_eeCorrection_offline_v2_LbyL2018V2"}
for label, tag in tags.items():
    process.GlobalTag.toGet.extend([
        cms.PSet(record = cms.string("GBRDWrapperRcd"),
            tag = cms.string(tag),
            connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
            label = cms.untracked.string(label)
            )
        ]
    )

# Change output collections
process.AODoutput.outputCommands += ['keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*']

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODoutput_step = cms.EndPath(process.AODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = 8
process.options.numberOfStreams = 0

from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
MassReplaceInputTag(process, new="rawDataMapperByLabel", old="rawDataCollector")

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run3

#call to customisation function customisePostEra_Run3 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run3(process)

# End of customisation functions


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
