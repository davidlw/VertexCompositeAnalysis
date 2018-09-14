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
    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 20
    config.Data.splitting = 'FileBased'
#    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0_v2-89e025b59ba99ac07dd655f3dba5c8df/USER'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = False
    config.Data.useParent = True
    config.Data.inputDBS = 'phys03'
    config.Site.storageSite = 'T2_US_MIT'

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

#    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_histogram_b1_v7'
#    config.JobType.psetName = '../test/scenarios/d0ana_default.py'
#    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0_v2-89e025b59ba99ac07dd655f3dba5c8df/USER'
#    config.Data.outputDatasetTag = 'pPb_PromptD0_BDT_default_histogram_v7'
#    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_hm185_ws_histogram_b1_v7'
    config.JobType.psetName = '../test/scenarios/d0ana_default_hm185_ws.py'
#    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0_v2-89e025b59ba99ac07dd655f3dba5c8df/USER'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    config.Data.outputDatasetTag = 'pPb_PromptD0_BDT_default_hm185_ws_histogram_v7'
    submit(config)

#    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_hm185_sb_histogram_b1_v7'
#    config.JobType.psetName = '../test/scenarios/d0ana_default_hm185_sb.py'
#    config.Data.outputDatasetTag = 'pPb_PromptD0_BDT_default_hm185_sb_histogram_v7'
#    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_hm185_ws_histogram_b2_v7'
    config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_hm185_ws_histogram_b3_v7'
    config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_hm185_ws_histogram_b4_v7'
    config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_hm185_ws_histogram_b5_v7'
    config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_hm185_ws_histogram_b6_v7'
    config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_PromptD0_BDT_default_hm185_ws_histogram_b1_v7'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    config.Data.outputDatasetTag = 'Pbp_PromptD0_BDT_default_hm185_ws_histogram_v7'
    submit(config)
