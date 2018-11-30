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
#    config.JobType.maxMemoryMB = 6000
#    config.JobType.maxJobRuntimeMin = 2750
    config.JobType.psetName = '../test/PbPbSkim2018_DiMuContBoth_cfg.py'
    config.Data.unitsPerJob = 5
    config.Data.totalUnits = 100
    config.Data.splitting = 'LumiBased'
#    config.Data.splitting = 'Automatic'
#    config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/DCSOnly/json_DCSONLY_HI.txt'
#    config.Data.lumiMask = 'Cert_326381-326618_HI_PromptReco_Collisions18_JSON_MuonPhys.txt'
    config.Data.lumiMask = 'json_DCSONLY_HI_327327.txt'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.publication = True
    config.Site.storageSite = 'T2_US_Vanderbilt'
#    config.Site.storageSite = 'T2_US_MIT'
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

    config.General.requestName = 'PbPb2018SkimDiMuDM_DCSOnly327327_test_v1'
    config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-PromptReco-v2/AOD'
    config.Data.outputDatasetTag = 'Skim_DiMuContBoth_DCSOnly327327_test_v1'
    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuDM_DCSOnly327327_v1'
#    config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-PromptReco-v2/AOD'
#    config.Data.outputDatasetTag = 'Skim_DiMuContBoth_DCSOnly327327_v1'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b0_DCSOnly327327_v1'
#    config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b1_DCSOnly327327_v1'
#    config.Data.inputDataset = '/HIMinimumBias1/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b2_DCSOnly327327_v1'
#    config.Data.inputDataset = '/HIMinimumBias2/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b3_DCSOnly327327_v1'
#    config.Data.inputDataset = '/HIMinimumBias3/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b4_DCSOnly327327_v1'
#    config.Data.inputDataset = '/HIMinimumBias4/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b0_MuonPhys_326381-326618_v1'
#    config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v1/AOD'
#    config.Data.outputDatasetTag = 'Skim_DiMuContBoth_MuonPhys_v1'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b1_MuonPhys_326381-326618_v1'
#    config.Data.inputDataset = '/HIMinimumBias1/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b2_MuonPhys_326381-326618_v1'
#    config.Data.inputDataset = '/HIMinimumBias2/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b3_MuonPhys_326381-326618_v1'
#    config.Data.inputDataset = '/HIMinimumBias3/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b4_MuonPhys_326381-326618_v1'
#    config.Data.inputDataset = '/HIMinimumBias4/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)

#    config.General.requestName = 'PbPb2018SkimDiMuMB_b5_MuonPhys_326381-326618_v1'
#    config.Data.inputDataset = '/HIMinimumBias5/HIRun2018A-PromptReco-v1/AOD'
#    submit(config)
