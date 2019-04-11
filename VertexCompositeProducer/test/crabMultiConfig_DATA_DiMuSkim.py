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
config.JobType.psetName = 'PbPbSkim2018_DiMuContBoth_cfg.py'
config.JobType.inputFiles = ['HeavyIonRPRcd_PbPb2018_offline.db']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt'
config.Data.runRange = '326381-327564'
config.Data.publication = True
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_FR_GRIF_LLR'

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
            "HIDoubleMuon_v1": { "PD": "/HIDoubleMuon/HIRun2018A-PromptReco-v1/AOD", "Units": 5, "Memory": 2400, "RunTime": 720 },
            "HIDoubleMuon_v2": { "PD": "/HIDoubleMuon/HIRun2018A-PromptReco-v2/AOD", "Units": 7, "Memory": 2400, "RunTime": 720 },
            "HISingleMuon_v1": { "PD": "/HISingleMuon/HIRun2018A-PromptReco-v1/AOD", "Units": 5, "Memory": 2400, "RunTime": 720 },
            "HISingleMuon_v2": { "PD": "/HISingleMuon/HIRun2018A-PromptReco-v2/AOD", "Units": 7, "Memory": 2400, "RunTime": 720 },
            "HIDoubleMuonPsiPeri_v2": { "PD": "/HIDoubleMuonPsiPeri/HIRun2018A-PromptReco-v2/AOD", "Units": 7, "Memory": 2400, "RunTime": 720 }
            }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeSkim_'+key+'_HIRun2018_DiMuMassMin7_20190407'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/user/%s/RiceHIN/PbPb2018/SKIM/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    submit(config)
