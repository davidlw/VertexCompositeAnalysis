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
    config.JobType.maxMemoryMB = 3000
    config.JobType.maxJobRuntimeMin = 2750
#    config.JobType.psetName = '../test/pPbFlowCorrSkim_2016_D0_cfg.py'
    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 3000
    config.Data.splitting = 'LumiBased'
    config.Data.inputDataset = '/PAHighMultiplicity1/PARun2016C-PromptReco-v1/AOD'
    config.Data.lumiMask = 'Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt'
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

    config.General.requestName = 'pPb2016_pPb_Skim_D0Both_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_D0Both_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_D0Both_v3_test'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DSToKsK_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_DSToKsK_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_DSToKsK_v3_test'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DSToPhiPi_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_DSToPhiPi_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_DSToPhiPi_v3_test'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DPM_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_DPM_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_DPM_v3_test'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_JPsiBoth_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_JPsiBoth_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_JPsiBoth_v3_test'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_MuMuBoth_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_MuMuBoth_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_MuMuBoth_v3_test'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_PhiBoth_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_PhiBoth_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_PhiBoth_v3_test'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_LambdaC_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_LambdaC_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_LambdaC_v3_test'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_V0_b1_v3_test'
    config.JobType.psetName = '../test/pPbSkim2016_V0_cfg.py'
    config.Data.outputDatasetTag = 'pPb_Skim_V0_v3_test'
    submit(config)
