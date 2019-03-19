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

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y1_v5'
config.JobType.psetName = '../test/PbPbSkimAndTreePhaseIIMTD_LamC3P_MC_cfg.py'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-3.1","yMax=-2.7"]
config.Data.inputDataset = '/MinBias_Hydjet_Drume5_5p5TeV_TuneCP5_Pythia8/PhaseIIMTDTDRAutumn18DR-NoPU_103X_upgrade2023_realistic_v2-v2/FEVT'
config.Data.allowNonValidInputDataset = True
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y2_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-2.7","yMax=-2.4"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y3_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-2.4","yMax=-2.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y4_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-2.1","yMax=-1.8"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y5_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-1.8","yMax=-1.5"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y6_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-1.5","yMax=-1.2"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y7_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-1.2","yMax=-0.9"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y8_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-0.9","yMax=-0.6"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y9_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-0.6","yMax=-0.3"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y10_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=-0.3","yMax=0.0"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y11_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=0.0","yMax=0.3"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y12_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=0.3","yMax=0.6"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y13_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=0.6","yMax=0.9"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y14_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=0.9","yMax=1.2"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y15_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=1.2","yMax=1.5"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y16_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=1.5","yMax=1.8"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y17_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=1.8","yMax=2.1"]
submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y18_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=2.1","yMax=2.4"]
submit(config)

#config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y19_v5'
#config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=2.4","yMax=2.7"]
#submit(config)

config.General.requestName = 'HyJets_mc_lamc3p_mtd_woPID_pt1p3to1p6_y20_v5'
config.JobType.pyCfgParams=["pTMin=1.3","pTMax=1.6","yMin=2.7","yMax=3.1"]
submit(config)
