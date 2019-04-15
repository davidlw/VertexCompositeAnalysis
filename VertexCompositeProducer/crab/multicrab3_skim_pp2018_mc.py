if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    from CRABClient.UserUtilities import config, getUsernameFromSiteDB
    config = config()

    config.General.workArea = 'VertexCompositeAna'
    config.General.transferOutputs = True
    config.General.transferLogs = False
    config.JobType.pluginName = 'Analysis'
#    config.JobType.maxMemoryMB = 3000
#    config.JobType.maxJobRuntimeMin = 2750
#    config.JobType.psetName = '../test/pPbFlowCorrSkim_2016_D0_cfg.py'
#    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 1000
    config.Data.splitting = 'Automatic'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.publication = True
    config.Site.storageSite = 'T2_US_MIT'
#    config.Site.storageSite = 'T3_US_Rice'

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

    config.General.requestName = 'pp2018MC_Skim_PromptD0Both_v1'
    config.JobType.psetName = '../test/ppSkim2018_D0Both_mc_cfg.py'
    config.Data.outputDatasetTag = 'Skim_D0Both_v1'
    config.Data.inputDataset = '/PrmtD0_pT-1p2_y-2p4_pp_13TeV_Pythia8/RunIILowPUAutumn18DR-102X_upgrade2018_realistic_v15-v1/AODSIM'
    submit(config)

    config.General.requestName = 'pp2018MC_Skim_NonPromptD0Both_v1'
    config.Data.outputDatasetTag = 'Skim_D0Both_v1'
    config.Data.inputDataset = '/NonPrD0_pT-1p2_y-2p4_pp_13TeV_pythia8/RunIILowPUAutumn18DR-102X_upgrade2018_realistic_v15-v1/AODSIM'
    submit(config)
