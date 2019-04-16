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
config.JobType.psetName = 'dimuana.py'

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_FR_*','T2_FR_*','T2_CH_CERN','T2_US_*']
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
            "HIDoubleMuon_v1": { "PD": "/HIDoubleMuon/anstahll-VertexCompositeSkim_HIDoubleMuon_v1_HIRun2018_DiMuMassMin7_20190412-8814e9757927b70b7a2a41b8684a900d/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "HIDoubleMuon_v2": { "PD": "/HIDoubleMuon/anstahll-VertexCompositeSkim_HIDoubleMuon_v2_HIRun2018_DiMuMassMin7_20190412-4d5a32383ad642c1b3ccbb603818f94f/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "HISingleMuon_v1": { "PD": "/HISingleMuon/anstahll-VertexCompositeSkim_HISingleMuon_v1_HIRun2018_DiMuMassMin7_20190412-7e2b6ec6e8df433328b3d7f08ea622ec/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "HISingleMuon_v2": { "PD": "/HISingleMuon/anstahll-VertexCompositeSkim_HISingleMuon_v2_HIRun2018_DiMuMassMin7_20190412-5057dee7bdc499f8cfd24e062bbe1180/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "HIDoubleMuonPsiPeri_v2": { "PD": "/HIDoubleMuonPsiPeri/anstahll-VertexCompositeSkim_HIDoubleMuonPsiPeri_v2_HIRun2018_DiMuMassMin7_20190412-b9b00062ed0736801c215e5f3b5ec90c/USER", "Units": 10, "Memory": 2400, "RunTime": 720 }
            }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_HIRun2018_DiMuMassMin7_20190414'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/PbPb2018/TREE/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    submit(config)
