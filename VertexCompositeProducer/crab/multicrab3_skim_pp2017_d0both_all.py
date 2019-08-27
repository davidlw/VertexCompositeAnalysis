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
    config.JobType.maxMemoryMB = 2000
    config.JobType.maxJobRuntimeMin = 2000
    config.JobType.psetName = '../test/ppLowPUSkim2017_D0Both_cfg.py'
    config.Data.unitsPerJob = 5
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

    config.section_('Site')
    config.Data.ignoreLocality = True
    config.Site.whitelist         = ['T2_CH_CERN','T1_UK_*','T2_UK_*','T3_UK_*','T1_DE_*','T2_DE_*','T1_IT_*','T2_IT_*','T1_FR_*','T2_FR_*','T1_US_*','T2_US_*','T3_US_*']

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

    config.General.requestName = 'pp2017C_Skim_D0Both_default_b1_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF1/Run2017C-17Nov2017-v1/AOD'
    config.Data.outputDatasetTag = 'pp2017C_Skim_D0Both_default_v1'
    submit(config)

    config.General.requestName = 'pp2017C_Skim_D0Both_default_b2_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF2/Run2017C-17Nov2017-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2017C_Skim_D0Both_default_b3_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF3/Run2017C-17Nov2017-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2017C_Skim_D0Both_default_b4_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF4/Run2017C-17Nov2017-v1/AOD'
    submit(config)

    config.General.requestName = 'pp2017C_Skim_D0Both_default_b5_v1'
    config.Data.inputDataset = '/HighMultiplicityEOF5/Run2017C-17Nov2017-v1/AOD'
    submit(config)
