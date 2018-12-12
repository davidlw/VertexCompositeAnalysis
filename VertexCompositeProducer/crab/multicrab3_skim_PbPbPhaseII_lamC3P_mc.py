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
    config.JobType.maxMemoryMB = 6000
#    config.JobType.maxJobRuntimeMin = 2750
    config.JobType.psetName = '../test/PbPbSkimPhaseII_LamC3P_MC_cfg.py'
    config.Data.unitsPerJob = 10
    config.Data.totalUnits = 30
    config.Data.splitting = 'FileBased'
#    config.Data.splitting = 'Automatic'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.publication = True
    config.Data.inputDBS = 'phys03'
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

    config.General.requestName = 'PbPbPhaseIISkimLamC3PMB_MC_test_v3'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/davidlw-LambdaC_PiKP_prompt_RECO_102X_test_v2-54694314416b05d17ca16e7e549fe18c/USER'
    config.Data.outputDatasetTag = 'Skim_LamC3P_test_v3'
    submit(config)
