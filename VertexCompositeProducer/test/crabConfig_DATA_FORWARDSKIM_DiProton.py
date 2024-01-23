# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from WMCore.Configuration import Configuration
#from CRABClient.UserUtilities import getUsername

config = Configuration()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375659_Muon.json'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'

## Submit the muon PDs
config.General.requestName = 'DiProton_HIForward0_HIRun2023A-PromptRec_20231209v1'
config.Data.inputDataset = '/HIForward0/HIRun2023A-PromptReco-v2/AOD'
config.Data.unitsPerJob = 20
config.JobType.maxMemoryMB = 4000
config.JobType.maxJobRuntimeMin = 2100
config.JobType.psetName = 'PbPbSkimAndTree2023_DiProton_ParticleAnalyzer_cfg.py'
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/group/phys_heavyions/davidlw/' 
