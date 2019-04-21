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
config.JobType.psetName = 'PbPbSkimAndTree2018_DiMuContBoth_mc_cfg.py'
config.JobType.inputFiles = ['HeavyIonRPRcd_PbPb2018_offline.db']

config.section_('Data')
#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=True)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {
            "Ups1S_TuneCP5_HydjetDrumMB_Pythia8": { "PD": "/Ups1SMM_5p02TeV_TuneCP5_Embd/anstahll-Ups1SMM_5p02TeV_TuneCP5_Embd_RECO_20190401-5db5dfa073297cb96791f14c622e83e2/USER", "Units": 40, "Memory": 1800, "RunTime": 720 },
            "Zmu10mu10_TuneCP5_HydjetDrumMB_Pythia8": { "PD": "/ZMM_5p02TeV_TuneCP5_Embd/anstahll-ZMM_5p02TeV_TuneCP5_Embd_RECO_22032019-5db5dfa073297cb96791f14c622e83e2/USER", "Units": 40, "Memory": 1800, "RunTime": 720 },
            #"Ups1S_TuneCP5_HydjetDrumMB_Pythia8": { "PD": "/Upsilon1S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM", "Units": 1, "Memory": 1800, "RunTime": 720 },
            #"Zmu10mu10_TuneCP5_HydjetDrumMB_Pythia8": { "PD": "/Zmu10mu10_pThat-0_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM", "Units": 1, "Memory": 1800, "RunTime": 720 },
            #"Zmu10mu10_TuneCP5_pythia8": { "PD": "/Zmu10mu10_TuneCP5_PbPb_5p02TeV_pythia8/HINPbPbAutumn18DR-NoPUmva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM", "Units": 1, "Memory": 1800, "RunTime": 720 },
            #"DYJetsToLL_TuneCP5_aMCatNLOFXFX-pythia8": { "PD": "/DYJetsToLL_TuneCP5_PbPb_aMCatNLOFXFX-pythia8/HINPbPbAutumn18DR-NoPUmva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM", "Units": 1, "Memory": 1800, "RunTime": 720 },
            }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_HINPbPbAutumn18_DiMuMC_20190421'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/PbPb2018/TREE/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
