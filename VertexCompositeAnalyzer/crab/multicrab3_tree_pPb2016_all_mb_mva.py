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
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 20
    config.Data.splitting = 'Automatic'
    config.JobType.psetName = '../test/D0/d0ana_default_mb_ws.py'
#    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0_v2-89e025b59ba99ac07dd655f3dba5c8df/USER'
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

    config.General.requestName = 'pPb2016_pPbMB_D0_BDT_default_hm185_ws_histogram_b1_v1'
    config.Data.inputDataset = '/PAMinimumBias1/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMB_D0_BDT_default_hm185_ws_histogram_b2_v1'
    config.Data.inputDataset = '/PAMinimumBias2/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMB_D0_BDT_default_hm185_ws_histogram_b3_v1'
    config.Data.inputDataset = '/PAMinimumBias3/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMB_D0_BDT_default_hm185_ws_histogram_b4_v1'
    config.Data.inputDataset = '/PAMinimumBias4/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMB_D0_BDT_default_hm185_ws_histogram_b5_v1'
    config.Data.inputDataset = '/PAMinimumBias5/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMB_D0_BDT_default_hm185_ws_histogram_b6_v1'
    config.Data.inputDataset = '/PAMinimumBias6/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMB_D0_BDT_default_hm185_ws_histogram_b7_v1'
    config.Data.inputDataset = '/PAMinimumBias7/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMB_D0_BDT_default_hm185_ws_histogram_b8_v1'
    config.Data.inputDataset = '/PAMinimumBias8/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMB_D0_BDT_default_hm185_ws_histogram_b1_v1'
    config.Data.inputDataset = '/PAMinimumBias1/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMB_D0_BDT_default_hm185_ws_histogram_b2_v1'
    config.Data.inputDataset = '/PAMinimumBias2/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMB_D0_BDT_default_hm185_ws_histogram_b3_v1'
    config.Data.inputDataset = '/PAMinimumBias3/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMB_D0_BDT_default_hm185_ws_histogram_b4_v1'
    config.Data.inputDataset = '/PAMinimumBias4/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMB_D0_BDT_default_hm185_ws_histogram_b5_v1'
    config.Data.inputDataset = '/PAMinimumBias5/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMB_D0_BDT_default_hm185_ws_histogram_b6_v1'
    config.Data.inputDataset = '/PAMinimumBias6/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMB_D0_BDT_default_hm185_ws_histogram_b7_v1'
    config.Data.inputDataset = '/PAMinimumBias7/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMB_D0_BDT_default_hm185_ws_histogram_b8_v1'
    config.Data.inputDataset = '/PAMinimumBias8/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)
