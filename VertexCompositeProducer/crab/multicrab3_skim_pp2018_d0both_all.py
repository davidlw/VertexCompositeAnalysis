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
    config.JobType.psetName = '../test/ppSkim2018_D0Both_cfg.py'
#    config.Data.unitsPerJob = 10
#    config.Data.totalUnits = 2000
#    config.Data.splitting = 'LumiBased'
    config.Data.splitting = 'Automatic'
#    config.Data.lumiMask = 'json_DCSONLY_SpecialRun.txt'
    config.Data.lumiMask = 'Cert_318939-319488_13TeV_PromptReco_SpecialCollisions18_JSON_LOWPU.txt'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = True
    config.Site.storageSite = 'T2_US_MIT'
#    config.Site.storageSite = 'T2_CH_CERN'

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

    config.General.requestName = 'pp2018_Skim_D0Both_default_b1_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF/Run2018B-PromptReco-v2/AOD'
    config.Data.outputDatasetTag = 'pp_Skim_D0Both_default_v1'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_default_b2_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF/Run2018C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_default_b3_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF0/Run2018C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_default_b4_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF1/Run2018C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_default_b5_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF2/Run2018C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_default_b6_v1'
    config.Data.inputDataset = '/DoubleMuonLowPU/Run2018C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_default_b7_v1'
    config.Data.inputDataset = '/DoubleMuon/Run2018B-PromptReco-v2/AOD'
    submit(config)
