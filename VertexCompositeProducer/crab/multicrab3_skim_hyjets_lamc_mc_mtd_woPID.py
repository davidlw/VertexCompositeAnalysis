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
config.Data.totalUnits = 200
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

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y1_v3'
config.JobType.psetName = '../test/PbPbSkimAndTreePhaseIIMTD_LamC3P_MC_cfg.py'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-3","yMax=-2.7"]
config.Data.inputDataset = '/MinBias_Hydjet_Drume5_5p5TeV_TuneCP5_Pythia8/PhaseIIMTDTDRAutumn18DR-NoPU_103X_upgrade2023_realistic_v3-v2/FEVT'
config.Data.allowNonValidInputDataset = True
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y2_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-2.7","yMax=-2.4"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y3_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-2.4","yMax=-2.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y4_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-2.1","yMax=-1.8"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y5_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-1.8","yMax=-1.5"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y6_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-1.5","yMax=-1.2"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y7_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-1.2","yMax=-0.9"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y8_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-0.9","yMax=-0.6"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y9_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-0.6","yMax=-0.3"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y10_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=-0.3","yMax=0.0"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y11_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=0.0","yMax=0.3"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y12_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=0.3","yMax=0.6"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y13_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=0.6","yMax=0.9"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y14_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=0.9","yMax=1.2"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y15_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=1.2","yMax=1.5"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y16_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=1.5","yMax=1.8"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y17_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=1.8","yMax=2.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y18_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=2.1","yMax=2.4"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y19_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=2.4","yMax=2.7"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt6toinf_y20_v3'
config.JobType.pyCfgParams=["pTMin=1.0","pTMax=1.3","yMin=2.7","yMax=3.0"]
submit(config)
