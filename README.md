# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_13_2_5_patch3

cd CMSSW_13_2_5_patch3/src

cmsenv

git clone -b ParticleFitter_13_2_X https://github.com/davidlw/VertexCompositeAnalysis

cd VertexCompositeAnalysis

scram b -j8

cd VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2023_Phi_ParticleAnalyzer_cfg.py
