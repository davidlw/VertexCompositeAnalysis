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
#    config.JobType.psetName = '../test/d0ana_mc_trainingtree_signal.py'
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 10
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

    config.General.requestName = 'pp2018_ppMC_PromptD0_MVATree_BDTCutNew_v3'
    config.Data.outputDatasetTag = 'pp2018_MVATree_BDTCutNew_v3'
    config.Data.inputDataset = '/PrmtD0_pT-1p2_y-2p4_pp_13TeV_Pythia8/davidlw-pp2018_SkimAndNtuple_D0Both_Prompt_v4-2dcf28a0a3592db40c70c00416240936/USER'
    config.JobType.psetName = '../test/D0/d0ana_mc_trainingtree_signal_bdt.py'
    submit(config)

    config.General.requestName = 'pp2018_ppMC_NonPromptD0_MVATree_BDTCutNew_v3'
    config.Data.outputDatasetTag = 'pp2018_MVATree_BDTCutNew_v3'
    config.Data.inputDataset = '/NonPrD0_pT-1p2_y-2p4_pp_13TeV_pythia8/davidlw-pp2018_SkimAndNtuple_D0Both_NonPrompt_v4-fe849d4090efcd3c73e5729216e7482c/USER'
    submit(config)
