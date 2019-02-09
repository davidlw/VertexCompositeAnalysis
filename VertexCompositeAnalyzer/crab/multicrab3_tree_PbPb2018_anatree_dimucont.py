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
#    config.JobType.maxMemoryMB = 8000
#    config.JobType.maxJobRuntimeMin = 2750
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 20
    config.Data.splitting = 'FileBased'
#    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = False
    config.Data.useParent = True
    config.Data.inputDBS = 'phys03'
    config.Site.storageSite = 'T2_CH_CERN'

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

    config.General.requestName = 'PbPb2018DiMuCont_anatree_b1_v1'
    config.JobType.psetName = '../test/JPsi/dimucontana_anatree_PbPb.py'
    config.Data.inputDataset = '/HIDoubleMuon/davidlw-2018Skim_DiMuCont_MuonPhysics_v2_NoEvtSel-825b193719cc9e390b7eb5f175a9d062/USER'
    config.Data.outputDatasetTag = 'AnaTree_DiMuContBoth_mass7toinf_v2'
    submit(config)

    config.General.requestName = 'PbPb2018DiMuCont_anatree_b2_v1'
    config.Data.inputDataset = '/HIDoubleMuon/davidlw-2018Skim_DiMuCont_MuonPhysics_v1_NoEvtSel-64fcce520b592387adbd248afe57eaa4/USER'
    config.Data.outputDatasetTag = 'AnaTree_DiMuContBoth_mass7toinf_v1'
    submit(config)

    config.General.requestName = 'PbPb2018DiMuCont_anatree_b1_peripheral_v1'
    config.JobType.psetName = '../test/JPsi/dimucontana_anatree_PbPb_peripheral.py'
    config.Data.inputDataset = '/HIDoubleMuonPsiPeri/davidlw-2018Skim_DiMuCont_MuonPhysics_v2_NoEvtSel-64fcce520b592387adbd248afe57eaa4/USER'
    config.Data.outputDatasetTag = 'AnaTree_DiMuContBoth_peripheral_v1'
    submit(config)
