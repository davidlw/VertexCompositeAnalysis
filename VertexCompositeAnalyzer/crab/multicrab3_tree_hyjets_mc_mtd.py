from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'VertexCompositeHyJetsAna'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
#    config.JobType.maxMemoryMB = 8000
#    config.JobType.maxJobRuntimeMin = 2750
#config.Data.unitsPerJob = 5
#    config.Data.totalUnits = 20
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.publication = False
#config.Data.useParent = True
config.Data.inputDBS = 'phys03'
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

config.General.requestName = 'hyjets_ntp_mc_mtd_v1'
config.JobType.psetName = '../test/HyJets_mc_mtd_ntuple.py'
config.Data.inputDataset = '/Hydjet_5p02TeV_TuneCP5_MTD/yousen-hyjets_mc_mtd_Skim_v1-58f154044687c63a0a655baa770ef20e/USER'
config.Data.outputDatasetTag = 'hyjets_mc_mtd_v1'
submit(config)