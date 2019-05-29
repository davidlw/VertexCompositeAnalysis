import FWCore.ParameterSet.Config as cms

from ZDC2018Gain_cfg import gainPreset

QWzdcreco = cms.EDProducer('QWZDC2018RecHit',
        srcDetId = cms.untracked.InputTag("zdcdigi", "DetId"),
        srcfC = cms.untracked.InputTag("zdcdigi", "regularfC"),
        ZDCCalib = gainPreset
        )

