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
#    config.JobType.maxJobRuntimeMin = 2750
#    config.JobType.psetName = '../test/d0ana_mc_trainingtree_signal.py'
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 5
    config.Data.splitting = 'FileBased'
#    config.Data.splitting = 'Automatic'
#    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = False
    config.Data.useParent = True
    config.Data.inputDBS = 'phys03'

    config.section_('Site')
    config.Data.ignoreLocality = True
    config.Site.whitelist = ['T1_US_*','T2_US_*','T2_CH_CERN','T2_BE_IIHE']
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

    config.General.requestName = 'pPb2016_pPbMC_PromptD0_MVATree_BDTCutNew_NoPreCuts_v3'
    config.Data.outputDatasetTag = 'pPb_MVATree_BDTCutNew_NoPreCuts_v3'
    config.Data.inputDataset = '/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/davidlw-pPb_Skim_D0Both_v1-b4e336f1831173b15564965986eff7fd/USER'
#    config.JobType.psetName = '../test/D0/d0ana_mc_trainingtree_signal_bdt.py'
    config.JobType.psetName = '../test/D0/d0ana_mc_trainingtree_signal_cutbased.py'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMC_NonPromptD0_MVATree_BDTCutNew_NoPreCuts_v3'
    config.Data.outputDatasetTag = 'pPb_MVATree_BDTCutNew_NoPreCuts_v3'
    config.Data.inputDataset = '/NonPromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/davidlw-pPb_Skim_D0Both_v1-b4e336f1831173b15564965986eff7fd/USER'
    submit(config)
