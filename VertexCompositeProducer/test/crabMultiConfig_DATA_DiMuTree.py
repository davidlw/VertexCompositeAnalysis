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
config.JobType.psetName = 'ppSkimAndTree2017_DiMuContBoth_cfg.py'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_299996-303819_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys_LowPU.txt'
config.Data.runRange = '299996-303819'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN']
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

dataMap = {}

for i in range(1,6):
    dataMap["HighMultiplicity"+str(i)+"_v2"] = { "PD": "/HighMultiplicityEOF"+str(i)+"/Run2017E-17Nov2017-v1/AOD", "Units": 7, "Memory": 2500, "RunTime": 1220 }
    dataMap["HighMultiplicity"+str(i)+"_v1"] = { "PD": "/HighMultiplicityEOF"+str(i)+"/Run2017C-17Nov2017-v1/AOD", "Units": 7, "Memory": 2500, "RunTime": 1220 }

for i in range(1,10):
    dataMap["MinimumBias"+str(i)+"_v2"] = { "PD": "/L1MinimumBias"+str(i)+"/Run2017E-PromptReco-v1/AOD", "Units": 7, "Memory": 2500, "RunTime": 1220 }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_Run2017LowPU_DiMuMassMin2_20190421'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/pp2017/TREE/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
