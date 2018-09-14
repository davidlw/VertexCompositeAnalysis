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
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 20
    config.Data.splitting = 'FileBased'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = False
    config.Data.useParent = True
    config.Data.inputDBS = 'phys03'
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

    config.General.requestName = 'pPb2016_pPb_PromptJPsi_BDT_default2_hm185_ws_histogram_b1_v3'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    config.JobType.psetName = '../test/JPsi/jpsiana_pPb_default2_hm185_ws.py'
    config.Data.outputDatasetTag = 'pPb_PromptJPsi_BDT_default2_hm185_ws_histogram_v3'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptJPsi_BDT_default2_hm185_ws_histogram_b2_v3'
    config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-pPb_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptJPsi_BDT_default2_hm185_ws_histogram_b3_v3'
    config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-pPb_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptJPsi_BDT_default2_hm185_ws_histogram_b4_v3'
    config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-pPb_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptJPsi_BDT_default2_hm185_ws_histogram_b5_v3'
    config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-pPb_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_pPb_PromptJPsi_BDT_default2_hm185_ws_histogram_b6_v3'
    config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-pPb_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_PromptJPsi_BDT_default2_hm185_ws_histogram_b1_v3'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-Pbp_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    config.JobType.psetName = '../test/JPsi/jpsiana_Pbp_default2_hm185_ws.py'
    config.Data.outputDatasetTag = 'Pbp_PromptJPsi_BDT_default2_hm185_ws_histogram_v3'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_PromptJPsi_BDT_default2_hm185_ws_histogram_b2_v3'
    config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-Pbp_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_PromptJPsi_BDT_default2_hm185_ws_histogram_b3_v3'
    config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-Pbp_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_PromptJPsi_BDT_default2_hm185_ws_histogram_b4_v3'
    config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-Pbp_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_PromptJPsi_BDT_default2_hm185_ws_histogram_b5_v3'
    config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-Pbp_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_PromptJPsi_BDT_default2_hm185_ws_histogram_b6_v3'
    config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-Pbp_Skim_JPsiBoth_default_v1-e0e26e1bcab095b6ddfe93551db7333d/USER'
    submit(config)
