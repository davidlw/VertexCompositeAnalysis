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
config.JobType.psetName = 'pPbSkim2016_DiMuContBoth_cfg.py'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'Cert_285479-286496_HI8TeV_PromptReco_pPb_Pbp_Collisions16_JSON_NoL1T_MuonPhys.txt'
config.Data.runRange = '285479-286496'
config.Data.publication = True
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T2_CH_CERN','T2_BE_IIHE']
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
            "PASingleMuon": { "PD": "/PASingleMuon/PARun2016C-PromptReco-v1/AOD", "Units": 7, "Memory": 2400, "RunTime": 720 },
            "PADoubleMuon": { "PD": "/PADoubleMuon/PARun2016C-PromptReco-v1/AOD", "Units": 7, "Memory": 2400, "RunTime": 720 }
          }
for i in range(1,20):
    dataMap[("PAMinimumBias"+str(i))] = { "PD": ("/PAMinimumBias"+str(i)+"/PARun2016C-PromptReco-v1/AOD"), "Units": 7, "Memory": 2400, "RunTime": 720 }

for i in range(0,7):
    dataMap[("PAHighMultiplicity"+str(i))] = { "PD": ("/PAHighMultiplicity"+str(i)+"/PARun2016C-PromptReco-v1/AOD"), "Units": 7, "Memory": 2400, "RunTime": 720 }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeSkim_'+key+'_PARun2016C_DiMuMassMin2_20190412'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/pPb2016/SKIM/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    submit(config)
