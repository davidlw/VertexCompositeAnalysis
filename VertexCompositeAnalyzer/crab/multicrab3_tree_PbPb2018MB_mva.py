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
    config.JobType.maxMemoryMB = 8000
#    config.JobType.maxJobRuntimeMin = 2750
    config.Data.unitsPerJob = 5
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

    config.General.requestName = 'PbPb2018_D0_BDT_ptsum2p2_pt1_ws_histogram_b0_v2'
    config.JobType.psetName = '../test/D0/d0ana_default_ws_PbPb.py'
    config.Data.inputDataset = '/HIMinimumBias0/davidlw-Skim_D0_Nov20DCS_v1-a1191996f743a6b381131d1462d45506/USER'
    config.Data.outputDatasetTag = 'D0_BDT_ptsum2p2_pt1_ws_histogram_DCSOnly327327_v2'
    submit(config)

    config.General.requestName = 'PbPb2018_D0_BDT_ptsum2p2_pt1_ws_histogram_b1_v2'
    config.Data.inputDataset = '/HIMinimumBias1/davidlw-Skim_D0_DCSOnly327327_v1-5ef0d9856259c20139aa4aaf1006bf0f/USER'
    submit(config)

    config.General.requestName = 'PbPb2018_D0_BDT_ptsum2p2_pt1_ws_histogram_b2_v2'
    config.Data.inputDataset = '/HIMinimumBias2/davidlw-Skim_D0_DCSOnly327327_v1-5ef0d9856259c20139aa4aaf1006bf0f/USER'
    submit(config)

    config.General.requestName = 'PbPb2018_D0_BDT_ptsum2p2_pt1_ws_histogram_b3_v2'
    config.Data.inputDataset = '/HIMinimumBias3/davidlw-Skim_D0_DCSOnly327327_v1-5ef0d9856259c20139aa4aaf1006bf0f/USER'
    submit(config)

    config.General.requestName = 'PbPb2018_D0_BDT_ptsum2p2_pt1_ws_histogram_b7_v2'
    config.Data.inputDataset = '/HIMinimumBias7/davidlw-Skim_D0_DCSOnly327327_v1-c2bc07f5e56345e5d3229c8dec729192/USER'
    submit(config)

