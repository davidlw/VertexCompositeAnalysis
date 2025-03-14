# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_15_0_0

cd CMSSW_15_0_0/src

cmsenv

git cms-merge-topic stahlleiton:ParticleAnalyzer_CMSSW_15_0_0

git clone git@github.com:stahlleiton/VertexCompositeAnalysis.git -b ParticleFitter_15_0_X

scram b -j8

cd VertexCompositeAnalysis/VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2025_SingleObject_MC_ParticleAnalyzer_cfg.py
