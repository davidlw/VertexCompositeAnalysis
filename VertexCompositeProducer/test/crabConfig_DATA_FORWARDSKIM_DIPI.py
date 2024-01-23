from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import getUsername

config = Configuration()

date='2023_10_13'
inputList = 'forward_highbetastar_skim.txt'
jobTag = 'ParticleAnalyzer_DiPi_HIRun2023A_dEdx_'+date+'_v1'
username = getUsername()

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = 'crab_projects/HIFOWARD_HIGHBETASTAR/DIPI/'+date
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PbPbSkimAndTree2023_DiPi_ParticleAnalyzer_cfg.py'
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 350
config.JobType.scriptExe = 'submitScript.sh'
config.JobType.inputFiles = ['emap_2023_newZDC_v3.txt', 'CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db']

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines()
config.Data.totalUnits = len(config.Data.userInputFiles)
#config.Data.inputDataset = '/Alternatively/DefineDataset/InsteadOf/InputFileList'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 14
config.Data.outLFNDirBase = '/store/group/phys_heavyions/' + username + '/CERN/PbPb2023/ParticleAnalyzer/DiPi/' + date + '/' + config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'
