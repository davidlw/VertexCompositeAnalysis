from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'VertexCompositePromptD0Ana_mtd'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 4000
config.JobType.maxJobRuntimeMin = 2750
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
#config.Data.inputDBS = 'phys03'
config.Site.storageSite = 'T2_US_MIT'

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

config.General.requestName = 'LambdaC_mc_mtd_v7'
config.JobType.psetName = '../test/LamC3P_PhaseIIMTD_mc_cfg.py'
config.Data.outputDatasetTag = 'lambdac_mc_mtd_Skim_pca2mm_v7'
config.Data.inputDataset = '/LambdaC_PiKP_prompt_pt1_y4_5p5TeV_TuneCP5_Pythia8/PhaseIIMTDTDRAutumn18DR-NoPU_103X_upgrade2023_realistic_v2-v1/FEVT'
submit(config)
