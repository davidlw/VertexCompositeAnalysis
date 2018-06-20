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
    config.Data.totalUnits = 100
    config.Data.splitting = 'FileBased'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.publication = False
    config.Data.useParent = True
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

    config.General.requestName = 'pPb2016_pPb_Tree_D0Both_b1_v4_test'
    config.JobType.psetName = '../test/d0ana.py'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0Both_v3_test-fc6e7002cc4077b6170121cd6a7f86ef/USER'
    config.Data.outputDatasetTag = 'pPb_Tree_D0Both_v4_test'
    submit(config)
