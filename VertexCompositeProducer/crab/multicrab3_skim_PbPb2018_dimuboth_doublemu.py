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
    config.JobType.inputFiles=['../test/HeavyIonRPRcd_PbPb2018_offline.db']
#    config.Data.unitsPerJob = 20
#    config.Data.totalUnits = 100
#    config.Data.splitting = 'LumiBased'
    config.Data.splitting = 'Automatic'
    config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327560_HI_PromptReco_Collisions18_JSON_MuonPhys.txt'
    #config.Data.lumiMask = 'Cert_326381-327489_HI_PromptReco_Collisions18_JSON_MuonPhys.txt'
#    config.Data.lumiMask = 'jsondiff_327560_327489.txt'
#    config.Data.lumiMask = 'json_DCSONLY_HI_327327.txt'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.publication = True
#    config.Site.storageSite = 'T2_US_Vanderbilt'
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

    config.General.requestName = 'PbPbDiMu_v1_MP327560_v2'
    config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-PromptReco-v1/AOD'
    config.Data.outputDatasetTag = '2018Skim_DiMuCont_MuonPhysics_v1_NoEvtSel'
    submit(config)

    config.General.requestName = 'PbPbDiMu_v2_MP327560_v2'
    config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-PromptReco-v2/AOD'
    config.Data.outputDatasetTag = '2018Skim_DiMuCont_MuonPhysics_v2_NoEvtSel'
    submit(config)
