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
    config.JobType.maxJobRuntimeMin = 2750
    config.JobType.psetName = '../test/ppLowPUSkim2017_D0Both_cfg.py'
    config.Data.unitsPerJob = 10
#    config.Data.totalUnits = 2000
    config.Data.splitting = 'LumiBased'
#    config.Data.splitting = 'Automatic'
    config.Data.lumiMask = 'Cert_295969-303819_13TeV_EOY2017ReReco_Collisions17_JSON_LowPU.txt'
#    config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_299996-303819_13TeV_EOY2017ReReco_Collisions17_JSON_LowPU.txt'
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

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b1_v1'
    config.Data.inputDataset = '/ZeroBias1/Run2017C-PromptReco-v2/AOD'
    config.Data.outputDatasetTag = 'pp2017C_Skim_D0Both_default_v1'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b2_v1'
    config.Data.inputDataset = '/ZeroBias2/Run2017C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b3_v1'
    config.Data.inputDataset = '/ZeroBias3/Run2017C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b4_v1'
    config.Data.inputDataset = '/ZeroBias4/Run2017C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b5_v1'
    config.Data.inputDataset = '/ZeroBias5/Run2017C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b6_v1'
    config.Data.inputDataset = '/ZeroBias6/Run2017C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b7_v1'
    config.Data.inputDataset = '/ZeroBias7/Run2017C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b8_v1'
    config.Data.inputDataset = '/ZeroBias8/Run2017C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b9_v1'
    config.Data.inputDataset = '/ZeroBias9/Run2017C-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'pp2017CZB_Skim_D0Both_default_b10_v1'
    config.Data.inputDataset = '/ZeroBias10/Run2017C-PromptReco-v2/AOD'
    submit(config)
