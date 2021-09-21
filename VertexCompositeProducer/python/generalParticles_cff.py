import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cfi import *

generalLambdaCandidatesNew = generalParticles.clone(
    pdgId = cms.int32(3122),
    mass = cms.double(1.115683),
    charge = cms.int32(0),
    doSwap = cms.bool(False),
    width = cms.double(0.03),
    vtxSortByTrkSize = cms.bool(False),

    preSelection = cms.string(""
       "charge==0"
       ),
    pocaSelection = cms.string(""
       "pt > 0."
       "&& userFloat('dca') > 0 && userFloat('dca') < 0.5"
       ),
    postSelection = cms.string(""
       "userFloat('normChi2') < 7.0"
       ),
    finalSelection = cms.string(""
       "userFloat('lVtxSig') > 4.5"
       "&& cos(userFloat('angle3D')) > 0.999"
       "&& abs(rapidity) < 2.4"
       ),
#
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1),
           selection = cms.string(
              "pt>0.0 && abs(eta)<2.4"
              "&& quality('loose')"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
              )
           ),
        cms.PSet(pdgId = cms.int32(2212), charge = cms.int32(+1),
           selection = cms.string(
              "pt>0.0 && abs(eta)<2.4"
              "&& quality('loose')"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
              )
           )
    ])
  )

generalAntiLambdaCandidatesNew = generalLambdaCandidatesNew.clone(

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(+1),
           selection = cms.string(
              "pt>0.0 && abs(eta)<2.4"
              "&& quality('loose')"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
              )
           ),
        cms.PSet(pdgId = cms.int32(2212), charge = cms.int32(-1),
           selection = cms.string(
              "pt>0.0 && abs(eta)<2.4"
              "&& quality('loose')"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
              )
           )
    ])
)

generalKshortCandidatesNew = generalParticles.clone(
    pdgId = cms.int32(310),
    mass = cms.double(0.497611),
    charge = cms.int32(0),
    doSwap = cms.bool(False),
    width = cms.double(0.05),
    vtxSortByTrkSize = cms.bool(False),

    preSelection = cms.string(""
       "charge==0"
       ),
    pocaSelection = cms.string(""
       "pt > 0."
       "&& userFloat('dca') > 0 && userFloat('dca') < 0.5"
       ),
    postSelection = cms.string(""
       "userFloat('normChi2') < 7.0"
       ),
    finalSelection = cms.string(""
       "userFloat('lVtxSig') > 4.5"
       "&& cos(userFloat('angle3D')) > 0.999"
       "&& abs(rapidity) < 2.4"
       ),
#
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1),
           selection = cms.string(
              "pt>0.0 && abs(eta)<2.4"
              "&& quality('loose')"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
#              '&& abs(userFloat("dzSig")) < 50.0 && abs(userFloat("dxySig")) < 50.0'
              )
           ),
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(+1),
           selection = cms.string(
              "pt>0.0 && abs(eta)<2.4"
              "&& quality('loose')"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
#              '&& abs(userFloat("dzSig")) < 50.0 && abs(userFloat("dxySig")) < 50.0'
              )
           )
    ])
  )

generalXiCandidatesNew = generalParticles.clone(

    mass = cms.double(1.32171),
    charge = cms.int32(1),
    width = cms.double(0.05),
    pdgId = cms.int32(3312),
    doSwap = cms.bool(False),
    vtxSortByTrkSize = cms.bool(False),

    preMassSelection = cms.string("abs(charge)==1"),

    postSelection = cms.string(""
       "userFloat('normChi2') < 7.0"
       ),

    finalSelection = cms.string(""
       "userFloat('lVtxSig') > 2." # 3 
       "&& cos(userFloat('angle3D')) > 0.9995"
       "&& abs(rapidity) < 2.4"
     ),

#    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2'),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(3122), source = cms.InputTag('generalLambdaCandidatesNew'), finalSelection = cms.string("userFloat('lVtxSig') > 3.0")), # 5
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1),
          selection = cms.string("pt>0. && abs(eta)<2.4"
              "&& quality('loose')"
              ),
          finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
#              '&& abs(userFloat("dzSig")) < 10.0 && abs(userFloat("dxySig")) < 10.0'
          )
        ),
    ])
)

generalAntiXiCandidatesNew = generalXiCandidatesNew.clone(

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(3122), source = cms.InputTag('generalAntiLambdaCandidatesNew'), finalSelection = cms.string("userFloat('lVtxSig') > 3.0")), # 5
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(+1),
          selection = cms.string("pt>0. && abs(eta)<2.4"
              "&& quality('loose')"
              ),
          finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
#              '&& abs(userFloat("dzSig")) < 10.0 && abs(userFloat("dxySig")) < 10.0'
          )
        ),
    ])
)

generalOmegaCandidatesNew = generalXiCandidatesNew.clone(
    pdgId = cms.int32(3334),
    mass = cms.double(1.67245),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(3122), source = cms.InputTag('generalLambdaCandidatesNew'), finalSelection = cms.string("userFloat('lVtxSig') > 3.0")), # 5
        cms.PSet(pdgId = cms.int32(321), charge = cms.int32(-1),
          selection = cms.string("pt>0. && abs(eta)<2.4"
              "&& quality('loose')"
              ),
          finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
#              '&& abs(userFloat("dzSig")) < 10.0 && abs(userFloat("dxySig")) < 10.0'
          )
        ),
    ])
)

generalAntiOmegaCandidatesNew = generalOmegaCandidatesNew.clone(

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(3122), source = cms.InputTag('generalAntiLambdaCandidatesNew'), finalSelection = cms.string("userFloat('lVtxSig') > 3.0")), # 5
        cms.PSet(pdgId = cms.int32(321), charge = cms.int32(+1),
          selection = cms.string("pt>0. && abs(eta)<2.4"
              "&& quality('loose')"
              ),
          finalSelection = cms.string(''
              'abs(userFloat("dzSig")) > 2.0 && abs(userFloat("dxySig")) > 2.0'
#              '&& abs(userFloat("dzSig")) < 10.0 && abs(userFloat("dxySig")) < 10.0'
          )
        ),
    ])
)
