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
    config.JobType.psetName = '../test/pPbSkim2016_DiMuContBoth_cfg.py'
#    config.Data.unitsPerJob = 10
#    config.Data.totalUnits = 2000
    config.Data.splitting = 'Automatic'
    config.Data.lumiMask = 'Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = True
    config.Site.storageSite = 'T2_US_MIT'
#    config.Site.storageSite = 'T2_CH_CERN'

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

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_mb_b1_v1'
    config.Data.inputDataset = '/PAMinimumBias1/PARun2016C-PromptReco-v1/AOD'
    config.Data.lumiMask = 'Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt'
    config.Data.outputDatasetTag = 'pPb_Skim_DiMuContBoth_default_v1'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_mb_b2_v1'
    config.Data.inputDataset = '/PAMinimumBias2/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_mb_b3_v1'
    config.Data.inputDataset = '/PAMinimumBias3/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_mb_b4_v1'
    config.Data.inputDataset = '/PAMinimumBias4/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_mb_b5_v1'
    config.Data.inputDataset = '/PAMinimumBias5/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_mb_b6_v1'
    config.Data.inputDataset = '/PAMinimumBias6/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_mb_b7_v1'
    config.Data.inputDataset = '/PAMinimumBias7/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_mb_b8_v1'
    config.Data.inputDataset = '/PAMinimumBias8/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_DoubleMu_b1_v1'
    config.Data.inputDataset = '/PADoubleMuon/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Skim_DiMuContBoth_default_SingleMu_b1_v1'
    config.Data.inputDataset = '/PASingleMuon/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_mb_b1_v1'
    config.Data.inputDataset = '/PAMinimumBias1/PARun2016C-PromptReco-v1/AOD'
    config.Data.lumiMask = 'Cert_285952-286496_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt'
    config.Data.outputDatasetTag = 'Pbp_Skim_DiMuContBoth_default_v1'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_mb_b2_v1'
    config.Data.inputDataset = '/PAMinimumBias2/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_mb_b3_v1'
    config.Data.inputDataset = '/PAMinimumBias3/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_mb_b4_v1'
    config.Data.inputDataset = '/PAMinimumBias4/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_mb_b5_v1'
    config.Data.inputDataset = '/PAMinimumBias5/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_mb_b6_v1'
    config.Data.inputDataset = '/PAMinimumBias6/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_mb_b7_v1'
    config.Data.inputDataset = '/PAMinimumBias7/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_mb_b8_v1'
    config.Data.inputDataset = '/PAMinimumBias8/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_DoubleMu_b1_v1'
    config.Data.inputDataset = '/PADoubleMuon/PARun2016C-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Skim_DiMuContBoth_default_SingleMu_b1_v1'
    config.Data.inputDataset = '/PASingleMuon/PARun2016C-PromptReco-v1/AOD'
    submit(config)
