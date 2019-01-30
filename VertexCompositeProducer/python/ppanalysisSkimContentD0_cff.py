import FWCore.ParameterSet.Config as cms

analysisSkimContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *',
      # event
#      'keep *_offlinePrimaryVertices_*_*',
      'keep *_offlinePrimaryVerticesRecovery_*_*',
#      'keep *_allVertices*_*_*',
#      'keep *_pixelVertices_*_*',
#      'keep *_pixel3Vertices_*_*',                                     
#      'keep *_sortedGoodVertices_*_*',
#      'keep *_bestVertex_*_*',
#      'keep *_offlineBeamSpot_*_*',
      'keep *_TriggerResults_*_*',
#      'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
#      'keep *_l1extraParticles_*_*',
      # jet
#      'keep *_towerMaker_*_*',
#      'keep *patJet*_selectedPatJets*_*_*',
#      'keep *PFJet*_ak5PFJets_*_*',
      # photons
#      'keep *_correctedHybridSuperClusters*_*_*',
#      'keep *_correctedMulti5x5SuperClustersWithPreshower*_*_*',
#      'keep *_photons*_*_*',
#      'keep *_conversions*_*_*',
      # muon
#      'keep recoMuons_muons_*_*',      
#      'keep *_hiCentrality_*_*',
      'keep *_hiEvtPlane_*_*',
      'keep *_hiEvtPlaneFlat_*_*',
      'keep *_centralityBin_*_*',
      # V0
      'keep *_generalD0Candidates*_D0_*',
      'keep *_generalD0Candidates*_MVAValues*_*',
      # LamC3P
      'keep *_generalLamC3PCandidates*_LamC3P_*',
      'keep *_generalLamC3PCandidates*_MVAValues*_*',
      # tracks
#      'keep recoTracks_generalTracks*_*_*',  
#      'keep recoTracks_pixelTracks_*_*',
#      'keep recoTracks_refitterForDeDx_*_*',
#      'keep recoTracks_selectTracks*_*_*',
#      'keep recoTracks_allTracks*_*_*',
#      'keep recoTracks_allVertices*_*_*',
#      'keep recoRecoChargedCandidates_allTracks*_*_*',
#      'keep *_dedx*_*_*',
#      'keep *DeDx*_*_*_*',
#      'keep *_pACentrality_*_*',
      # clusters
#      'keep *_siPixelClusters_*_*',       
      # PFCandidate
#      'keep *_particleFlow_*_*',
      # mc (if present)
#      'keep GenEventInfoProduct_*_*_*',
#      'keep *_*GenJet*_*_*',
#      'keep reco*_genParticles_*_*',
#      'keep *_generator_*_*',
#      'keep *SimVertex*_*_*_*',
      # missing ET
#      'keep *_pfMet*_*_*',
      # sim track matching
#      'keep *_trackingParticleRecoTrackAsssociation_*_*',
#      'keep TrackingParticles_mergedtruth_MergedTrackTruth_*'
      )
    )
