from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.workArea = 'VertexCompositeAna'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/pPbSkimAndTree2016_D0Both_CutBased_cfg.py'
#config.JobType.psetName = '../test/pPbSkimAndTree2016_D0Both_BDT_cfg.py'
config.JobType.priority = 80

config.section_('Data')
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = 1000
#config.Data.lumiMask = ''
#config.Data.runRange = '285479-286496'
config.Data.publication = False
#config.Data.useParent = True
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Site.whitelist = ['T2_US_Vanderbilt','T2_US_MIT']
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T1_IT_*','T2_IT_*','T2_CH_CERN']
#config.Site.whitelist         = ['T2_CH_CERN','T1_UK_*','T2_UK_*','T3_UK_*','T1_DE_*','T2_DE_*','T1_IT_*','T2_IT_*','T1_FR_*','T2_FR_*','T1_US_*','T2_US_*','T3_US_*']
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_Rice'

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
for i in range(1,2):
    dataMap[("PAHighMultiplicity"+str(i)+"_pPb")] = { "PD": ("/PAHighMultiplicity"+str(i)+"/PARun2016C-PromptReco-v1/AOD"), "Units": 4, "Memory": 4000, "RunTime": 2000, "LumiMask": 'Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt' }
#for i in range(1,2):
#    dataMap[("PAHighMultiplicity"+str(i)+"_Pbp")] = { "PD": ("/PAHighMultiplicity"+str(i)+"/PARun2016C-PromptReco-v1/AOD"), "Units": 1, "Memory": 4000, "RunTime": 2000, "LumiMask": 'Cert_285952-286496_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt' }

#for i in range(1,3):
#    dataMap[("PAMinimumBias"+str(i)+"_pPb")] = { "PD": ("/PAMinimumBias"+str(i)+"/PARun2016C-PromptReco-v1/AOD"), "Units": 10, "Memory": 4000, "RunTime": 2000, "LumiMask": 'Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt' }
#for i in range(1,3):
#    dataMap[("PAMinimumBias"+str(i)+"_Pbp")] = { "PD": ("/PAMinimumBias"+str(i)+"/PARun2016C-PromptReco-v1/AOD"), "Units": 10, "Memory": 4000, "RunTime": 2000, "LumiMask": 'Cert_285952-286496_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt' }

## Submit the muon PDs
for key, val in dataMap.items():
#    config.General.requestName = 'VertexCompositeTree_'+key+'_PARun2016C_D0AndTrackTreeMVA_20190915'
#    config.General.requestName = 'VertexCompositeTree_'+key+'_PARun2016C_D0AndTrackTreeMVA_NoErrHitDA2D_20190914'
#    config.General.requestName = 'VertexCompositeTree_'+key+'_PARun2016C_D0AndTrackTreeMVA_NoPreCuts_20190914v2'
    config.General.requestName = 'VertexCompositeTree_'+key+'_PARun2016C_D0AndTrackTreeMVA_PreCuts_20190916'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.Data.lumiMask = val["LumiMask"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
