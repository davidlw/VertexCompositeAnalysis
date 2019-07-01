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
    config.JobType.maxMemoryMB = 5000
#    config.JobType.maxJobRuntimeMin = 2750
    config.Data.unitsPerJob = 1
    config.Data.totalUnits = -1
    config.Data.splitting = 'FileBased'
#    config.Data.runRange = '285505-285505'
#    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = False
    config.Data.useParent = True
    config.Data.inputDBS = 'phys03'

    config.section_('Site')
    config.Data.ignoreLocality = True
    config.Site.whitelist = ['T3_US_*','T1_US_*','T2_US_*','T2_CH_CERN','T2_IT_*','T1_IT_*','T2_FR_*','T1_FR_*']
    config.Site.storageSite = 'T2_CH_CERN'
    config.JobType.allowUndistributedCMSSW = True
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

    config.General.requestName = 'pPb2016_pPb_Tree_D0_default_BDTCutNewAndTrack_b1_v4'
    config.JobType.psetName = '../test/D0/d0ana_mvatracktree_bdtcut.py'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    config.Data.outputDatasetTag = 'pPb_Tree_D0_default_BDTCutNewAndTrack_v4'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Tree_D0_default_BDTCutNewAndTrack_b2_v4'
    config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Tree_D0_default_BDTCutNewAndTrack_b3_v4'
    config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Tree_D0_default_BDTCutNewAndTrack_b4_v4'
    config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Tree_D0_default_BDTCutNewAndTrack_b5_v4'
    config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_Tree_D0_default_BDTCutNewAndTrack_b6_v4'
    config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Tree_D0_default_BDTCutNewAndTrack_b1_v4'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    config.Data.outputDatasetTag = 'Pbp_Tree_D0_default_BDTCutNewAndTrack_v4'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Tree_D0_default_BDTCutNewAndTrack_b2_v4'
    config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Tree_D0_default_BDTCutNewAndTrack_b3_v4'
    config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Tree_D0_default_BDTCutNewAndTrack_b4_v4'
    config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Tree_D0_default_BDTCutNewAndTrack_b5_v4'
    config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Tree_D0_default_BDTCutNewAndTrack_b6_v4'
    config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
    submit(config)
