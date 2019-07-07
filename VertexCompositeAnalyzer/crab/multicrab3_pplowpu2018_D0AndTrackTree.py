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
config.JobType.psetName = '../test/D0/d0ana_mvatracktree_bdtcut.py'
config.JobType.priority = 100

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
#config.Data.totalUnits = 1
#config.Data.lumiMask = ''
#config.Data.runRange = '285479-286496'
config.Data.publication = False
config.Data.useParent = True
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
config.Site.whitelist         = ['T2_BR_*','T2_CH_CERN','T1_UK_*','T2_UK_*','T3_UK_*','T1_DE_*','T2_DE_*','T1_IT_*','T2_IT_*','T1_FR_*','T2_FR_*','T1_US_*','T2_US_*','T3_US_*']
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

dataMap = {
            "HighMultiplicityEOF": { "PD": "/HighMultiplicityEOF/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
            "HighMultiplicityEOF0": { "PD": "/HighMultiplicityEOF0/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
            "HighMultiplicityEOF1": { "PD": "/HighMultiplicityEOF1/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "HighMultiplicityEOF2": { "PD": "/HighMultiplicityEOF2/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "ZeroBias_pp2018C": { "PD": "/ZeroBias/davidlw-pp2018_Skim_D0Both_default_v1-91cff2ce5544be4dc902cb93cefcf11e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "ZeroBias_pp2018B": { "PD": "/ZeroBias/davidlw-pp2018B_Skim_D0Both_default_v1-91cff2ce5544be4dc902cb93cefcf11e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
          }

#for i in range(0,10):
#    dataMap[("MinimumBias"+str(i)+"_pp2018C")] = { "PD": ("/MinimumBias"+str(i)+"/davidlw-pp_Skim_D0Both_default_v1-91cff2ce5544be4dc902cb93cefcf11e/USER"), "Units": 1, "Memory": 2800, "RunTime": 4000 }
#for i in range(0,10):
#    dataMap[("MinimumBias"+str(i)+"_pp2018B")] = { "PD": ("/MinimumBias"+str(i)+"/davidlw-pp2018B_Skim_D0Both_default_v1-91cff2ce5544be4dc902cb93cefcf11e/USER"), "Units": 1, "Memory": 2800, "RunTime": 4000 }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_Run2018_D0AndTrackTreeMVA_20190704'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
