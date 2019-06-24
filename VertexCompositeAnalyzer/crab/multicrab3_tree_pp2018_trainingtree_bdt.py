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
#    config.Data.totalUnits = 100
    config.Data.splitting = 'FileBased'
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

    config.General.requestName = 'pp2018_Tree_D0_default_BDTCutNew_b1_v3'
    config.JobType.psetName = '../test/D0/d0ana_trainingtree_bdtcut.py'
    config.Data.inputDataset = '/HighMultiplicityEOF/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER'
    config.Data.outputDatasetTag = 'pp_Tree_D0_default_BDTCutNew_v3'
    submit(config)

    config.General.requestName = 'pp2018_Tree_D0_default_BDTCutNew_b2_v3'
    config.Data.inputDataset = '/HighMultiplicityEOF0/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER'
    submit(config)

    config.General.requestName = 'pp2018_Tree_D0_default_BDTCutNew_b3_v3'
    config.Data.inputDataset = '/HighMultiplicityEOF1/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER'
    submit(config)

    config.General.requestName = 'pp2018_Tree_D0_default_BDTCutNew_b4_v3'
    config.Data.inputDataset = '/HighMultiplicityEOF2/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER'
    submit(config)
