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
#    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = -1
    config.Data.splitting = 'Automatic'
#    config.Data.runRange = '285505-285505'
#    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = False
    config.Data.useParent = True
    config.Data.inputDBS = 'phys03'
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

    config.General.requestName = 'TreeDiMContBoth_mb_b1_v4'
    config.JobType.psetName = '../test/JPsi/dimucontana_trainingtree.py'
    config.Data.inputDataset = '/PAMinimumBias1/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    config.Data.outputDatasetTag = 'pPb_Tree_DiMuContBoth_training_v4'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_DM_v4'
    config.Data.inputDataset = '/PADoubleMuon/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_SM_v4'
    config.Data.inputDataset = '/PASingleMuon/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_mb_b2_v4'
    config.Data.inputDataset = '/PAMinimumBias2/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_mb_b3_v4'
    config.Data.inputDataset = '/PAMinimumBias3/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_mb_b4_v4'
    config.Data.inputDataset = '/PAMinimumBias4/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_mb_b5_v4'
    config.Data.inputDataset = '/PAMinimumBias5/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_mb_b6_v4'
    config.Data.inputDataset = '/PAMinimumBias6/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_mb_b7_v4'
    config.Data.inputDataset = '/PAMinimumBias7/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth_mb_b8_v4'
    config.Data.inputDataset = '/PAMinimumBias8/davidlw-pPb_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_mb_b1_v4'
    config.Data.inputDataset = '/PAMinimumBias1/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    config.Data.outputDatasetTag = 'Pbp_Tree_DiMuContBoth_training_v4'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_DM_v4'
    config.Data.inputDataset = '/PADoubleMuon/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_SM_v4'
    config.Data.inputDataset = '/PASingleMuon/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_mb_b2_v4'
    config.Data.inputDataset = '/PAMinimumBias2/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_mb_b3_v4'
    config.Data.inputDataset = '/PAMinimumBias3/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_mb_b4_v4'
    config.Data.inputDataset = '/PAMinimumBias4/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_mb_b5_v4'
    config.Data.inputDataset = '/PAMinimumBias5/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_mb_b6_v4'
    config.Data.inputDataset = '/PAMinimumBias6/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_mb_b7_v4'
    config.Data.inputDataset = '/PAMinimumBias7/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)

    config.General.requestName = 'TreeDiMContBoth2_mb_b8_v4'
    config.Data.inputDataset = '/PAMinimumBias8/davidlw-Pbp_Skim_DiMuContBoth_default_v1-753e5416dcc725930129c4cde0e65308/USER'
    submit(config)
