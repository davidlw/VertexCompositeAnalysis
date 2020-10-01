import FWCore.ParameterSet.Config as cms

analysisSkimContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *',
      'keep recoVertexs_offlinePrimaryVertices__*',
      'keep recoVertexs_*_vertices_*',
      'keep recoBeamSpot_offlineBeamSpot__*',
      'keep edmTriggerResults_TriggerResults__*',
      'keep recoDeDxDataedmValueMap_*__*SKIM',
      'keep triggerTriggerEvent_hltTriggerSummaryAOD__*',
      'keep Global*BlkBXVector_gtStage2Digis__*',
      'keep recoCentrality_*__*',
      'keep int_centralityBin_*_*',
      'keep recoEvtPlanes_hiEvtPlaneFlat__*',
      'keep OnlineLuminosityRecord_onlineMetaDataDigis__*',
      'keep LumiScalerss_scalersRawToDigi__*',
      'keep *_externalLHEProducer_*_*',
      'keep GenEventInfoProduct_generator__*',
      'keep recoGenParticles_genParticles__*',
      'keep patGenericParticles_*_*_*',
      'keep intedmValueMap_nTracks__*',
      'keep LumiInfo_*_*_*',
      )
    )
