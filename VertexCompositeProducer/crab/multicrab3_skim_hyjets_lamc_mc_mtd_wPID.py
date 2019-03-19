from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'VertexCompositeHyJetsAna_mtd'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 27000
config.JobType.maxJobRuntimeMin = 2750
config.Data.unitsPerJob = 5
config.Data.totalUnits = 5000
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.publication = True
#config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Site.whitelist = ['T2_CH_*' , 'T2_US_*' , 'T2_FR_*']
#config.Site.storageSite = 'T2_US_MIT'
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/MTD/%s/' % (getUsernameFromSiteDB())
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

config.General.requestName = 'HyJets_mc_lamc3p_mtd_wPID_1p5RMS_pt1to1p5_v6'
config.JobType.psetName = '../test/PbPbSkimAndTreePhaseIIMTD_LamC3P_MC_wPID_cfg.py'
config.JobType.pyCfgParams=["pTMin=0.97","pTMax=1.5","yMin=-3.1","yMax=3.1"]
config.Data.inputDataset = '/MinBias_Hydjet_Drume5_5p5TeV_TuneCP5_Pythia8/PhaseIIMTDTDRAutumn18DR-NoPU_103X_upgrade2023_realistic_v2-v2/FEVT'
config.Data.allowNonValidInputDataset = True
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_wPID_1p5RMS_pt1p5to2_v6'
config.JobType.pyCfgParams=["pTMin=1.5","pTMax=2.0","yMin=-3.1","yMax=3.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_wPID_1p5RMS_pt2to2p5_v6'
config.JobType.pyCfgParams=["pTMin=2.0","pTMax=2.5","yMin=-3.1","yMax=3.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_wPID_1p5RMS_pt2p5to3_v6'
config.JobType.pyCfgParams=["pTMin=2.5","pTMax=3.0","yMin=-3.1","yMax=3.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_wPID_1p5RMS_pt3to3p5_v6'
config.JobType.pyCfgParams=["pTMin=3","pTMax=3.5","yMin=-3.1","yMax=3.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_wPID_1p5RMS_pt3p5to4_v6'
config.JobType.pyCfgParams=["pTMin=3.5","pTMax=4.0","yMin=-3.1","yMax=3.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_wPID_1p5RMS_pt4to6_v6'
config.JobType.pyCfgParams=["pTMin=4","pTMax=6","yMin=-3.1","yMax=3.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_wPID_1p5RMS_pt6toinf_v6'
config.JobType.pyCfgParams=["pTMin=6","pTMax=100000.0","yMin=-3.1","yMax=3.1"]
submit(config)
