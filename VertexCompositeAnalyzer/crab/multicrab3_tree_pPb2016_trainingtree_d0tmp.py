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
#    config.JobType.maxMemoryMB = 5000
#    config.JobType.maxJobRuntimeMin = 2750
    config.Data.unitsPerJob = 2
    config.Data.totalUnits = 650
    config.Data.splitting = 'FileBased'
#    config.Data.runRange = '285505-285505'
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

#    config.General.requestName = 'pPb2016_pPb_Tree_D0_default_NoBDTCutScale07_run285505_b1_v1'
#    config.JobType.psetName = '../test/D0/d0ana_trainingtree.py'
#    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER'
#    config.Data.outputDatasetTag = 'pPb_Tree_D0_default_NoBDTCutScale07_run285505_v1'
#    submit(config)

    config.General.requestName = 'pPb2016_pPb_Tree_D0_default_pt6to8_b1_v5'
    config.JobType.psetName = '../test/D0/d0ana_trainingtree_ptcut.py'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-pPb_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER' #1000 files
    config.Data.outputDatasetTag = 'pPb_Tree_D0_default_pt6to8_v5'
    submit(config)

    config.General.requestName = 'pPb2016_Pbp_Tree_D0_default_pt1p5to4_b1_v5'
    config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-Pbp_Skim_D0Both_default_v1-1c60067745b1f8eb361c5a1e6ce2795e/USER' #1185 files
    config.Data.outputDatasetTag = 'Pbp_Tree_D0_default_pt1p5to4_v5'
    submit(config)
