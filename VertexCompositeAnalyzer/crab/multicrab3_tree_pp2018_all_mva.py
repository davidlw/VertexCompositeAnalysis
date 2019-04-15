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
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 20
    config.Data.splitting = 'FileBased'
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

#    config.General.requestName = 'pPb2016_pPb_PromptD0_BDT_default_histogram_b1_v7'
#    config.JobType.psetName = '../test/scenarios/d0ana_default.py'
#    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0_v2-89e025b59ba99ac07dd655f3dba5c8df/USER'
#    config.Data.outputDatasetTag = 'pPb_PromptD0_BDT_default_histogram_v7'
#    submit(config)

    config.General.requestName = 'pp2018_D0_BDT_default_hm100_ws_histogram_b1_v4'
    config.JobType.psetName = '../test/D0/d0ana_default_pphm100_ws.py'
    config.Data.inputDataset = '/HighMultiplicityEOF/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER'
    submit(config)

    config.General.requestName = 'pp2018_D0_BDT_default_hm100_ws_histogram_b2_v4'
    config.Data.inputDataset = '/HighMultiplicityEOF0/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER'
    submit(config)

    config.General.requestName = 'pp2018_D0_BDT_default_hm100_ws_histogram_b3_v4'
    config.Data.inputDataset = '/HighMultiplicityEOF1/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER'
    submit(config)

    config.General.requestName = 'pp2018_D0_BDT_default_hm100_ws_histogram_b4_v4'
    config.Data.inputDataset = '/HighMultiplicityEOF2/davidlw-pp_Skim_D0Both_default_v1-f9516adb0ab4bcb5c36605171e6520f9/USER'
    submit(config)
