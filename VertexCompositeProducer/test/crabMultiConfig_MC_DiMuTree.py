from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PbPbSkimAndTree2015_DiMuContBoth_mc_cfg.py'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = False

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T1_FR_*','T2_US_*','T2_FR_*','T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {
            "Pythia8_Zmu10mu10_Hydjet_MB": { "PD": "/Pythia8_Zmu10mu10_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13_ext1-v1/AODSIM", "Units": 1, "Memory": 1100, "RunTime": 180 },
            "Pythia8_Ups1SMM_ptUps_00_03_Hydjet_MB": { "PD": "/Pythia8_Ups1SMM_ptUps_00_03_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM", "Units": 1, "Memory": 1100, "RunTime": 180 },
            "Pythia8_Ups1SMM_ptUps_03_06_Hydjet_MB": { "PD": "/Pythia8_Ups1SMM_ptUps_03_06_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM", "Units": 1, "Memory": 1100, "RunTime": 180 },
            "Pythia8_Ups1SMM_ptUps_06_09_Hydjet_MB": { "PD": "/Pythia8_Ups1SMM_ptUps_06_09_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM", "Units": 1, "Memory": 1100, "RunTime": 180 },
            "Pythia8_Ups1SMM_ptUps_09_12_Hydjet_MB": { "PD": "/Pythia8_Ups1SMM_ptUps_09_12_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM", "Units": 1, "Memory": 1100, "RunTime": 180 },
            "Pythia8_Ups1SMM_ptUps_12_15_Hydjet_MB": { "PD": "/Pythia8_Ups1SMM_ptUps_12_15_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM", "Units": 1, "Memory": 1100, "RunTime": 180 },
            "Pythia8_Ups1SMM_ptUps_15_30_Hydjet_MB": { "PD": "/Pythia8_Ups1SMM_ptUps_15_30_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM", "Units": 1, "Memory": 1100, "RunTime": 180 },
            "Pythia8_Ups1SMM_ptUps_30_inf_Hydjet_MB": { "PD": "/Pythia8_Ups1SMM_ptUps_30_inf_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM", "Units": 1, "Memory": 1100, "RunTime": 180 },
          }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_HIFall15_DiMuMC_20190514'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/PbPb2015/TREE/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
