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
            "PADoubleMuon": { "PD": "/PADoubleMuon/anstahll-VertexCompositeSkim_PADoubleMuon_PARun2016C_DiMuMassMin2_20190412-4a47d065114b374d5497b4cb35cf53a2/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PASingleMuon": { "PD": "/PASingleMuon/anstahll-VertexCompositeSkim_PASingleMuon_PARun2016C_DiMuMassMin2_20190412-8e37ccf40e37cfab655c6fcf90ec04d7/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PASingleMuon_v2": { "PD": "/PASingleMuon/anstahll-VertexCompositeSkim_PASingleMuon_PARun2016C_DiMuMassMin2_20190412-ca8b98cedaf8a4d64d1d48a0ed50024f/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity0": { "PD": "/PAHighMultiplicity0/anstahll-VertexCompositeSkim_PAHighMultiplicity0_PARun2016C_DiMuMassMin2_20190412-f7399c2a174ee6ef7ee6469827d7114f/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity0_v2": { "PD": "/PAHighMultiplicity0/anstahll-VertexCompositeSkim_PAHighMultiplicity0_PARun2016C_DiMuMassMin2_20190412-d4771bc870adfcd441bf0c06583eba83/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity1": { "PD": "/PAHighMultiplicity1/anstahll-VertexCompositeSkim_PAHighMultiplicity1_PARun2016C_DiMuMassMin2_20190412-2d45ed28388e01651a083214a1c990bb/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity2": { "PD": "/PAHighMultiplicity2/anstahll-VertexCompositeSkim_PAHighMultiplicity2_PARun2016C_DiMuMassMin2_20190412-20b8e9aa8941bf044810e6a467a832f0/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity3": { "PD": "/PAHighMultiplicity3/anstahll-VertexCompositeSkim_PAHighMultiplicity3_PARun2016C_DiMuMassMin2_20190412-a6d40e88ec89c63bb06f3f75318748f9/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity4": { "PD": "/PAHighMultiplicity4/anstahll-VertexCompositeSkim_PAHighMultiplicity4_PARun2016C_DiMuMassMin2_20190412-b5c8028a7e0719689e4c419c7b7825bf/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity5": { "PD": "/PAHighMultiplicity5/anstahll-VertexCompositeSkim_PAHighMultiplicity5_PARun2016C_DiMuMassMin2_20190412-46625dcbf2e35c65661adec1b0e15427/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity6": { "PD": "/PAHighMultiplicity6/anstahll-VertexCompositeSkim_PAHighMultiplicity6_PARun2016C_DiMuMassMin2_20190412-a9db28195d9e069f8bb8d36cbef6d10c/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity6_v2": { "PD": "/PAHighMultiplicity6/anstahll-VertexCompositeSkim_PAHighMultiplicity6_PARun2016C_DiMuMassMin2_20190412-f7399c2a174ee6ef7ee6469827d7114f/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAHighMultiplicity7": { "PD": "/PAHighMultiplicity7/anstahll-VertexCompositeSkim_PAHighMultiplicity7_PARun2016C_DiMuMassMin2_20190412-3154688b915ffca6b9b5da6e522b4370/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias1": { "PD": "/PAMinimumBias1/anstahll-VertexCompositeSkim_PAMinimumBias1_PARun2016C_DiMuMassMin2_20190412-a251e958057ac5ccc6a704c4082074fe/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias2": { "PD": "/PAMinimumBias2/anstahll-VertexCompositeSkim_PAMinimumBias2_PARun2016C_DiMuMassMin2_20190412-c689f38ddce83af2d258d8a701777b95/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias3": { "PD": "/PAMinimumBias3/anstahll-VertexCompositeSkim_PAMinimumBias3_PARun2016C_DiMuMassMin2_20190412-f60e46c534b612dd96502a40c2880a55/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias4": { "PD": "/PAMinimumBias4/anstahll-VertexCompositeSkim_PAMinimumBias4_PARun2016C_DiMuMassMin2_20190412-b13bc4d247420a20265ffad758782fe9/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias5": { "PD": "/PAMinimumBias5/anstahll-VertexCompositeSkim_PAMinimumBias5_PARun2016C_DiMuMassMin2_20190412-90798ba6dbe76e655d7f90d1967f1349/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias6": { "PD": "/PAMinimumBias6/anstahll-VertexCompositeSkim_PAMinimumBias6_PARun2016C_DiMuMassMin2_20190412-3908184e90249a470d4b4a16b02b7b1c/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias7": { "PD": "/PAMinimumBias7/anstahll-VertexCompositeSkim_PAMinimumBias7_PARun2016C_DiMuMassMin2_20190412-a3087abefccbeb021b7391ec9d0db759/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias8": { "PD": "/PAMinimumBias8/anstahll-VertexCompositeSkim_PAMinimumBias8_PARun2016C_DiMuMassMin2_20190412-3154688b915ffca6b9b5da6e522b4370/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias9": { "PD": "/PAMinimumBias9/anstahll-VertexCompositeSkim_PAMinimumBias9_PARun2016C_DiMuMassMin2_20190412-5fed557a66354841b445cd850c5c8d97/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias9_v2": { "PD": "/PAMinimumBias9/anstahll-VertexCompositeSkim_PAMinimumBias9_PARun2016C_DiMuMassMin2_20190412-4a47d065114b374d5497b4cb35cf53a2/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias10": { "PD": "/PAMinimumBias10/anstahll-VertexCompositeSkim_PAMinimumBias10_PARun2016C_DiMuMassMin2_20190412-d4771bc870adfcd441bf0c06583eba83/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias10_v2": { "PD": "PAMinimumBias10/anstahll-VertexCompositeSkim_PAMinimumBias10_PARun2016C_DiMuMassMin2_20190412-a6d40e88ec89c63bb06f3f75318748f9/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias11": { "PD": "/PAMinimumBias11/anstahll-VertexCompositeSkim_PAMinimumBias11_PARun2016C_DiMuMassMin2_20190412-ccdb5604392def6bcb0eaefad2a4aa26/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias12": { "PD": "/PAMinimumBias12/anstahll-VertexCompositeSkim_PAMinimumBias12_PARun2016C_DiMuMassMin2_20190412-ca8b98cedaf8a4d64d1d48a0ed50024f/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias12_v2": { "PD": "/PAMinimumBias12/anstahll-VertexCompositeSkim_PAMinimumBias12_PARun2016C_DiMuMassMin2_20190412-20b8e9aa8941bf044810e6a467a832f0/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias13": { "PD": "/PAMinimumBias13/anstahll-VertexCompositeSkim_PAMinimumBias13_PARun2016C_DiMuMassMin2_20190412-fe44f7e1a56db474c9811c8c5dabcdb0/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias13_v2": { "PD": "/PAMinimumBias13/anstahll-VertexCompositeSkim_PAMinimumBias13_PARun2016C_DiMuMassMin2_20190412-2d45ed28388e01651a083214a1c990bb/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias14": { "PD": "/PAMinimumBias14/anstahll-VertexCompositeSkim_PAMinimumBias14_PARun2016C_DiMuMassMin2_20190412-644b1d1504eda34b3f1682ed4aa53851/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias14_v2": { "PD": "/PAMinimumBias14/anstahll-VertexCompositeSkim_PAMinimumBias14_PARun2016C_DiMuMassMin2_20190412-5fed557a66354841b445cd850c5c8d97/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias15": { "PD": "/PAMinimumBias15/anstahll-VertexCompositeSkim_PAMinimumBias15_PARun2016C_DiMuMassMin2_20190412-668a70768a81b521e8d4e9456b20ca4d/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias15_v2": { "PD": "/PAMinimumBias15/anstahll-VertexCompositeSkim_PAMinimumBias15_PARun2016C_DiMuMassMin2_20190412-a9db28195d9e069f8bb8d36cbef6d10c/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias16": { "PD": "/PAMinimumBias16/anstahll-VertexCompositeSkim_PAMinimumBias16_PARun2016C_DiMuMassMin2_20190412-db0313d6fb727cf42d19099a55cb44c5/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias16_v2": { "PD": "/PAMinimumBias16/anstahll-VertexCompositeSkim_PAMinimumBias16_PARun2016C_DiMuMassMin2_20190412-46625dcbf2e35c65661adec1b0e15427/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias17": { "PD": "/PAMinimumBias17/anstahll-VertexCompositeSkim_PAMinimumBias17_PARun2016C_DiMuMassMin2_20190412-043db8602b0f6247afa7c0f9758a8d30/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias17_v2": { "PD": "/PAMinimumBias17/anstahll-VertexCompositeSkim_PAMinimumBias17_PARun2016C_DiMuMassMin2_20190412-b5c8028a7e0719689e4c419c7b7825bf/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias18": { "PD": "/PAMinimumBias18/anstahll-VertexCompositeSkim_PAMinimumBias18_PARun2016C_DiMuMassMin2_20190412-5114cbc3fbbcaf515a6f90200cc37b53/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias18_v2": { "PD": "/PAMinimumBias18/anstahll-VertexCompositeSkim_PAMinimumBias18_PARun2016C_DiMuMassMin2_20190412-fe44f7e1a56db474c9811c8c5dabcdb0/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias19": { "PD": "/PAMinimumBias19/anstahll-VertexCompositeSkim_PAMinimumBias19_PARun2016C_DiMuMassMin2_20190412-e65cf435681cbc97b9d35a84ccf9fba0/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias20": { "PD": "/PAMinimumBias20/anstahll-VertexCompositeSkim_PAMinimumBias20_PARun2016C_DiMuMassMin2_20190412-4a47d065114b374d5497b4cb35cf53a2/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            "PAMinimumBias20_v2": { "PD": "/PAMinimumBias20/anstahll-VertexCompositeSkim_PAMinimumBias20_PARun2016C_DiMuMassMin2_20190412-90798ba6dbe76e655d7f90d1967f1349/USER", "Units": 10, "Memory": 2400, "RunTime": 720 },
            }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeTree_'+key+'_PARun2016C_DiMuMassMin2_20190416'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/pPb2016/TREE/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    submit(config)
