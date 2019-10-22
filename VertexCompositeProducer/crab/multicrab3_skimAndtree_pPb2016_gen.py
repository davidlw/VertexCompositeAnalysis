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
    config.JobType.priority = 200
#    config.JobType.maxMemoryMB = 3000
#    config.JobType.maxJobRuntimeMin = 2750
#    config.JobType.psetName = '../test/pPbFlowCorrSkim_2016_D0_cfg.py'
    config.JobType.allowUndistributedCMSSW = True
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 5 
    config.Data.splitting = 'FileBased'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = False
#    config.Data.inputDBS = 'phys03'
#    config.Site.storageSite = 'T2_CH_CERN'
    config.Site.storageSite = 'T3_US_Rice'

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

    config.General.requestName = 'pPb2016MC_pPbSkimAndTree_D0Both_NonPromptGEN_20191021_v3'
    config.JobType.psetName = '../test/pPbSkimAndTree2016_D0Both_gen_cfg.py'
    config.Data.outputDatasetTag = 'pPbGEN_SkimAndTree_D0Both_20191021_v3'
    config.Data.inputDataset = '/NonPromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v1/AODSIM'
    submit(config)
