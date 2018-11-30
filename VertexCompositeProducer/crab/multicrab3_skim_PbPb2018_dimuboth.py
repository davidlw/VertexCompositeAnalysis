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
    config.JobType.psetName = '../test/PbPbSkimAndTree2018_DiMuCont_cfg.py'
    config.Data.unitsPerJob = 20
#    config.Data.totalUnits = 100
    config.Data.splitting = 'LumiBased'
#    config.Data.splitting = 'Automatic'
#    config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/DCSOnly/json_DCSONLY_HI.txt'
    config.Data.lumiMask = 'Cert_326381-326618_HI_PromptReco_Collisions18_JSON_MuonPhys.txt'
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

    config.General.requestName = 'PbPbDiMuMB_rb9_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat9/HIRun2018A-PromptReco-v1/AOD'
    config.Data.outputDatasetTag = 'Skim_DiMuCont_MuonPhysics_v1'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb8_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat8/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb7_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat7/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb6_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat6/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb5_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat5/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb4_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat4/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb3_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat3/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat2/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb11_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat11/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb10_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat10/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb1_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat1/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb0v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat0/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_rb0_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBiasReducedFormat0/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b9v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias9/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b9_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias9/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b8v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias8/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b8_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias8/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b7v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias7/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b7_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias7/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b6v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias6/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b6_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias6/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b5v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias5/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b5_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias5/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b4v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias4/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b4_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias4/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b3v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias3/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b3_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias3/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b2v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias2/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias2/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b19v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias19/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b19_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias19/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b18v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias18/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b18_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias18/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b17v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias17/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b17_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias17/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b16v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias16/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b16_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias16/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b15v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias15/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b15_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias15/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b14v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias14/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b14_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias14/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b13v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias13/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b13_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias13/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b12v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias12/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b12_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias12/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b11v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias11/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b11_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias11/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b10v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias10/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b10_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias10/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b1v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias1/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b1_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias1/HIRun2018A-PromptReco-v1/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b0v2_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v2/AOD'
    submit(config)

    config.General.requestName = 'PbPbDiMuMB_b0_MP326618_v1'
    config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v1/AOD'
    submit(config)
