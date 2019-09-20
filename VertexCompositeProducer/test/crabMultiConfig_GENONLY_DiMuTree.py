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
config.JobType.psetName = 'pPbSkimAndTree2016_DiMuContBoth_GENONLY_cfg.py'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_KR_*']
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
            "Ups1SToMuMu_pPb-Bst_GENonly": { "PD": "/Upsilon1S_pPb-Bst_8p16-Pythia8/pPb816Spring16GS-pPbBst_GENonly_80X_mcRun2_pA_v4-v4/GEN", "Units": 1, "Memory": 1600, "RunTime": 820 },
            "Ups2SToMuMu_pPb-Bst_GENonly": { "PD": "/Upsilon2S_pPb-Bst_8p16-Pythia8/pPb816Spring16GS-pPbBst_GENonly_80X_mcRun2_pA_v4-v8/GEN", "Units": 1, "Memory": 1600, "RunTime": 820 },
            "Ups3SToMuMu_pPb-Bst_GENonly": { "PD": "/Upsilon3S_pPb-Bst_8p16-Pythia8/pPb816Spring16GS-pPbBst_GENonly_80X_mcRun2_pA_v4-v4/GEN", "Units": 1, "Memory": 1600, "RunTime": 820 },
            "JPsiToMuMu_pPb-Bst_GENonly": { "PD": "/Psi1S_pPb-Bst_8p16-Pythia8/pPb816Spring16GS-pPbBst_GENonly_80X_mcRun2_pA_v4-v4/GEN", "Units": 1, "Memory": 1600, "RunTime": 820 },
            "Psi2SToMuMu_pPb-Bst_GENonly": { "PD": "/Psi2S_pPb-Bst_8p16-Pythia8/pPb816Spring16GS-pPbBst_GENonly_80X_mcRun2_pA_v4-v9/GEN", "Units": 1, "Memory": 1600, "RunTime": 820 },
            "BToPsiToMuMu_pPb-Bst_GENonly": { "PD": "/NonPrPsi1_2S_pPb-Bst_8p16-Pythia8/pPb816Spring16GS-pPbBst_GENonly_80X_mcRun2_pA_v4-v7/GEN", "Units": 1, "Memory": 1600, "RunTime": 820 },
          }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_pPb816Summer16_DiMuGENONLY_20190705'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/pPb2016/TREE/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
