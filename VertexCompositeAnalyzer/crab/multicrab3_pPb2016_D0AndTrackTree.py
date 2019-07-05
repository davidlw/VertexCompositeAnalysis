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
#config.Data.totalUnits = 10
#config.Data.lumiMask = ''
#config.Data.runRange = '285479-286496'
config.Data.publication = False
config.Data.useParent = True
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN']
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
#            "PAHighMultiplicity1_pPb_test": { "PD": "/PAHighMultiplicity1/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity2_pPb": { "PD": "/PAHighMultiplicity2/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity3_pPb": { "PD": "/PAHighMultiplicity3/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity4_pPb": { "PD": "/PAHighMultiplicity4/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity5_pPb": { "PD": "/PAHighMultiplicity5/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity6_pPb": { "PD": "/PAHighMultiplicity6/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity7_pPb": { "PD": "/PAHighMultiplicity7/davidlw-pPb_Skim_D0Both_default_v1-970f442569e397a6dbe2e80979c0eea4/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity0_pPb": { "PD": "/PAHighMultiplicity0/davidlw-pPb_Skim_D0Both_default_v1-970f442569e397a6dbe2e80979c0eea4/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity1_Pbp": { "PD": "/PAHighMultiplicity1/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity2_Pbp": { "PD": "/PAHighMultiplicity2/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity3_Pbp": { "PD": "/PAHighMultiplicity3/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity4_Pbp": { "PD": "/PAHighMultiplicity4/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity5_Pbp": { "PD": "/PAHighMultiplicity5/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity6_Pbp": { "PD": "/PAHighMultiplicity6/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity7_Pbp": { "PD": "/PAHighMultiplicity7/davidlw-Pbp_Skim_D0Both_default_v1-970f442569e397a6dbe2e80979c0eea4/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
#            "PAHighMultiplicity0_Pbp": { "PD": "/PAHighMultiplicity0/davidlw-Pbp_Skim_D0Both_default_v1-970f442569e397a6dbe2e80979c0eea4/USER", "Units": 1, "Memory": 2800, "RunTime": 4000 },
          }

for i in range(4,7):
    dataMap[("PAMinimumBias"+str(i)+"_pPb")] = { "PD": ("/PAMinimumBias"+str(i)+"/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER"), "Units": 1, "Memory": 2800, "RunTime": 4000 }
#for i in range(1,9):
#    dataMap[("PAMinimumBias"+str(i)+"_Pbp")] = { "PD": ("/PAMinimumBias"+str(i)+"/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER"), "Units": 1, "Memory": 2800, "RunTime": 4000 }
#for i in range(9,21):
#    dataMap[("PAMinimumBias"+str(i)+"_Pbp")] = { "PD": ("/PAMinimumBias"+str(i)+"/davidlw-Pbp_Skim_D0Both_default_v1-970f442569e397a6dbe2e80979c0eea4/USER"), "Units": 1, "Memory": 2800, "RunTime": 4000 }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_PARun2016C_D0AndTrackTreeMVA_20190704'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
