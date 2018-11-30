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
#    config.JobType.maxMemoryMB = 8000
#    config.JobType.maxJobRuntimeMin = 2750
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 20
    config.Data.splitting = 'FileBased'
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

#    config.General.requestName = 'PbPb2018DiMuContBothDCS_training_b1_v1'
#    config.JobType.psetName = '../test/JPsi/dimucontana_trainingtree_PbPb.py'
#    config.Data.inputDataset = '/HIMinimumBias0/davidlw-Skim_DiMuContBoth_Nov20DCS_v1-c4a8626d21e8a6f68e7db48052d953a7/USER'
#    config.Data.outputDatasetTag = 'Tree_DiMuContBoth_training_Nov20DCS_v1'
#    submit(config)

#    config.General.requestName = 'PbPb2018DiMuContBothDCSNov27_training_b1_v1'
#    config.JobType.psetName = '../test/JPsi/dimucontana_trainingtree_PbPb.py'
#    config.Data.inputDataset = '/HIDoubleMuon/davidlw-Skim_DiMuContBoth_DCSOnly327327_v1-5a6e4168a05c7a459d9e19c05b7139bf/USER'
#    config.Data.outputDatasetTag = 'Tree_DiMuContBoth_training_Nov27DCS_v1'
#    submit(config)

    config.General.requestName = 'PbPb2018MBDiMuContBothDCSNov27_training_test_b1_v1'
    config.JobType.psetName = '../test/JPsi/dimucontana_trainingtree_PbPb.py'
    config.Data.inputDataset = '/HIDoubleMuon/davidlw-Skim_DiMuContBoth_DCSOnly327327_test_v1-5a6e4168a05c7a459d9e19c05b7139bf/USER'
    config.Data.outputDatasetTag = 'Tree_DiMuContBoth_training_Nov27DCS_VUtest_v1'
    submit(config)

#    config.General.requestName = 'PbPb2018MBDiMuContBothDCSNov27_training_b1_v1'
#    config.JobType.psetName = '../test/JPsi/dimucontana_trainingtree_PbPb.py'
#    config.Data.inputDataset = '/HIMinimumBias1/davidlw-Skim_DiMuContBoth_DCSOnly327327_v1-5a6e4168a05c7a459d9e19c05b7139bf/USER'
#    config.Data.outputDatasetTag = 'Tree_DiMuContBoth_training_Nov27DCS_v1'
#    submit(config)

#    config.General.requestName = 'PbPb2018MBDiMuContBothDCSNov27_training_b0_v1'
#    config.Data.inputDataset = '/HIMinimumBias0/davidlw-Skim_DiMuContBoth_DCSOnly327327_v1-5a6e4168a05c7a459d9e19c05b7139bf/USER'
#    submit(config)

#    config.General.requestName = 'PbPb2018MBDiMuContBothDCSNov27_training_b2_v1'
#    config.Data.inputDataset = '/HIMinimumBias2/davidlw-Skim_DiMuContBoth_DCSOnly327327_v1-5a6e4168a05c7a459d9e19c05b7139bf/USER'
#    submit(config)

#    config.General.requestName = 'PbPb2018MBDiMuContBothDCSNov27_training_b3_v1'
#    config.Data.inputDataset = '/HIMinimumBias3/davidlw-Skim_DiMuContBoth_DCSOnly327327_v1-5a6e4168a05c7a459d9e19c05b7139bf/USER'
#    submit(config)

#    config.General.requestName = 'PbPb2018MBDiMuContBothDCSNov27_training_b4_v1'
#    config.Data.inputDataset = '/HIMinimumBias4/davidlw-Skim_DiMuContBoth_DCSOnly327327_v1-5a6e4168a05c7a459d9e19c05b7139bf/USER'
#    submit(config)
