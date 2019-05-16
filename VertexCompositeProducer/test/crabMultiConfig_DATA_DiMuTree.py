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
config.JobType.psetName = 'ppSkimAndTree2018_DiMuContBoth_cfg.py'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_318939-319488_13TeV_PromptReco_SpecialCollisions18_JSON_MuonPhys_LOWPU.txt'
config.Data.runRange = '318939-319488'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN']
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

dataMap = {}

dataMap["HighMultiplicity_Run2018B"] = { "PD": "/HighMultiplicityEOF/Run2018B-PromptReco-v2/AOD", "Units": 7, "Memory": 2500, "RunTime": 1220 }
dataMap["HighMultiplicity_Run2018C"] = { "PD": "/HighMultiplicityEOF/Run2018C-PromptReco-v2/AOD", "Units": 7, "Memory": 2500, "RunTime": 1220 }
for i in range(0,3):
    dataMap["HighMultiplicity"+str(i)+"_Run2018C"] = { "PD": "/HighMultiplicityEOF"+str(i)+"/Run2018C-PromptReco-v2/AOD", "Units": 7, "Memory": 2500, "RunTime": 1220 }

for i in range(0,9):
    dataMap["MinimumBias"+str(i)+"_Run2018B"] = { "PD": "/MinimumBias"+str(i)+"/Run2018B-PromptReco-v2/AOD", "Units": 7, "Memory": 2500, "RunTime": 1220 }
    dataMap["MinimumBias"+str(i)+"_Run2018C"] = { "PD": "/MinimumBias"+str(i)+"/Run2018C-PromptReco-v2/AOD", "Units": 7, "Memory": 2500, "RunTime": 1220 }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_Run2018LowPU_DiMuMassMin2_20190514'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/pp2018/TREE/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
