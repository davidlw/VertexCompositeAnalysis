from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException  # updated import for Python 3

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config
config = config()

#inputList = 'filelist_HIForward0.txt'

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['emap_2023_newZDC_v3.txt']

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
#config.Data.inputDBS = 'global'
#config.Data.splitting = 'LumiBased'
#config.Data.totalUnits = 5000
#config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Golden.json'
#config.Data.runRange = '374288-375823'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))  # updated for Python 3

    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))  # updated for Python 3

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {
#            "HIForward": { "PD": "/HIForward/HIRun2018A-04Apr2019-v1/AOD", "Units": 30, "Memory": 1800, "RunTime": 1400, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_ALLDIMU_cfg.py" },
            "OO2025": { "PD": "/IonPhysics/anstahll-crab_Run394153_OORun2025-PromptReco-v1_Run3_2025_UPC_OXY_2025_07_15-0ea11c841dcbf2e6da5b143561e71ec4/USER", "Units": 2, "Memory": 2500, "RunTime": 1400, "PSet": "PbPbSkimAndTree2023_Phi_ParticleAnalyzer_MiniAOD_cfg.py" },
            }

#for i in range(0,1):
#    dataMap[("HIForward"+str(i))] = { "PD": ("/HIForward"+str(i)+"/HIRun2023A-16Jan2024-v1/AOD"), "Units": 25, "Memory": 4000, "RunTime": 2100, "PSet": "PbPbSkimAndTree2023_DiMuCont_ParticleAnalyzer_cfg.py" } # UCC

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'Phi_'+key+'_OOSkimAndTree2025_20250819v1'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.JobType.psetName = val["PSet"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/davidlw/' 

    print("Submitting CRAB job for: " + val["PD"])  # updated for Python 3

    submit(config)
