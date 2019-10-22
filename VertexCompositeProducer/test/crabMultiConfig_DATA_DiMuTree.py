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
config.JobType.inputFiles = ['HeavyIonRPRcd_PbPb2018_offline.db']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt'
config.Data.runRange = '326381-327564'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T3_US_Rice'

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
#            "HIDoubleMuon": { "PD": "/HIDoubleMuon/HIRun2018A-04Apr2019-v1/AOD", "Units": 5, "Memory": 2500, "RunTime": 1800, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_cfg.py" },
#            "HISingleMuon": { "PD": "/HISingleMuon/HIRun2018A-04Apr2019-v1/AOD", "Units": 5, "Memory": 2500, "RunTime": 1800, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_cfg.py" },
#            "HIDoubleMuonPsiPeri": { "PD": "/HIDoubleMuonPsiPeri/HIRun2018A-04Apr2019-v1/AOD", "Units": 5, "Memory": 2500, "RunTime": 1800, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_cfg.py" },
#            "HIForward": { "PD": "/HIForward/HIRun2018A-04Apr2019-v1/AOD", "Units": 10, "Memory": 1800, "RunTime": 1400, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_ALLDIMU_cfg.py" },
            }

for i in range(4,7):
    dataMap[("HIMinimumBias"+str(i))] = { "PD": ("/HIMinimumBias"+str(i)+"/HIRun2018A-04Apr2019-v1/AOD"), "Units": 20, "Memory": 2000, "RunTime": 2000, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_JPsiMB_cfg.py" }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_HIRun2018_04Apr2019_DiMuMassMin2p5_Cent30100_20191010'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.JobType.psetName = val["PSet"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/PbPb2018/TREE/%s' % (getUsernameFromSiteDB(), config.General.requestName)

    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
