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
    config.JobType.maxMemoryMB = 2800
    config.JobType.maxJobRuntimeMin = 4000
    config.JobType.psetName = '../test/ppSkim2018_D0Both_cfg.py'
    config.JobType.priority = 100
    config.JobType.allowUndistributedCMSSW = True

    config.Data.unitsPerJob = 5
#    config.Data.totalUnits = 2000
    config.Data.splitting = 'LumiBased'
#    config.Data.splitting = 'Automatic'
    config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_318939-319488_13TeV_PromptReco_SpecialCollisions18_JSON_LOWPU.txt'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.allowNonValidInputDataset = True 
    config.Data.ignoreLocality = True
    config.Data.publication = True

    config.section_('Site')
    config.Site.whitelist = ['T1_UK_*','T2_UK_*','T3_UK_*','T2_IT_*','T1_IT_*','T2_FR_*','T1_FR_*','T3_US_*','T1_US_*','T2_US_*','T2_CH_CERN']
    config.Site.storageSite = 'T2_US_MIT'

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

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b0_missing_v2'
    config.Data.inputDataset = '/MinimumBias0/Run2018B-17Sep2018-v1/AOD'
    config.Data.outputDatasetTag = 'pp2018BRereco_Skim_D0Both_default_v1'
    config.Data.lumiMask = 'notprocessedLumis_b0.json'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b1_v2'
    config.Data.inputDataset = '/MinimumBias1/Run2018B-17Sep2018-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b2_missing_v2'
    config.Data.inputDataset = '/MinimumBias2/Run2018B-17Sep2018-v1/AOD'
    config.Data.lumiMask = 'notprocessedLumis_b2.json'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b3_v2'
    config.Data.inputDataset = '/MinimumBias3/Run2018B-17Sep2018-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b4_v2'
    config.Data.inputDataset = '/MinimumBias4/Run2018B-17Sep2018-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b5_v2'
    config.Data.inputDataset = '/MinimumBias5/Run2018B-17Sep2018-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b6_v2'
    config.Data.inputDataset = '/MinimumBias6/Run2018B-17Sep2018-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b7_v2'
    config.Data.inputDataset = '/MinimumBias7/Run2018B-17Sep2018-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b8_missing_v2'
    config.Data.inputDataset = '/MinimumBias8/Run2018B-17Sep2018-v1/AOD'
    config.Data.lumiMask = 'notprocessedLumis_b8.json'
    submit(config)

    config.General.requestName = 'pp2018BMB_Skim_D0Both_default_b9_missing_v2'
    config.Data.inputDataset = '/MinimumBias9/Run2018B-17Sep2018-v1/AOD'
    config.Data.lumiMask = 'notprocessedLumis_b9.json'
    submit(config)
