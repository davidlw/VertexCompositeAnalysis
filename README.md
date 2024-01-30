# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_13_2_6_patch2

cd CMSSW_13_2_6_patch2/src

cmsenv

git clone -b ParticleFitter_13_2_X https://github.com/davidlw/VertexCompositeAnalysis

cd VertexCompositeAnalysis

scram b -j8

cd VertexCompositeProducer/test

# Phi reconstruction
cmsRun PbPbSkimAndTree2023_Phi_ParticleAnalyzer_cfg.py

# V0 reconstruction
cmsRun PbPbSkimAndTree2023_V0_ParticleAnalyzer_cfg.py
