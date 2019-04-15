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
#    config.JobType.psetName = '../test/pPbFlowCorrSkim_2016_D0_cfg.py'
    config.Data.unitsPerJob = 2
#    config.Data.totalUnits = 1000
    config.Data.splitting = 'FileBased'
#    config.Data.splitting = 'Automatic'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.inputDBS = 'phys03'
    config.Data.publication = True
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

    config.General.requestName = 'pp2018_Skim_D0Both_Prompt_b1_v3'
    config.JobType.psetName = '../test/ppSkim2018_D0Both_mc_cfg.py'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_Prompt_b1_v3'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-prompt_pt1p2_y2p4_pp2018_RECO_b1_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_Prompt_b2_v3'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_Prompt_b2_v3'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-prompt_pt1p2_y2p4_pp2018_RECO_b2_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_Prompt_b3_v3'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_Prompt_b3_v3'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-prompt_pt1p2_y2p4_pp2018_RECO_b3_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b1_v2'
    config.JobType.psetName = '../test/ppSkim2018_D0Both_mc_cfg.py'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b1_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b1_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b2_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b2_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b2_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b3_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b3_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b3_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b4_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b4_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b4_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b5_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b5_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b5_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b6_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b6_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b6_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b7_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b7_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b7_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b8_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b8_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b8_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b9_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b9_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b9_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)

    config.General.requestName = 'pp2018_Skim_D0Both_NonPrompt_b10_v2'
    config.Data.outputDatasetTag = 'pp2018_Skim_D0Both_NonPrompt_b10_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_pp2018_RECO_b10_v2-992f8f980a519dc2ad8fd8d00bdc2fa2/USER'
    submit(config)
