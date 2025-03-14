# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_15_0_0

cd CMSSW_15_0_0/src

cmsenv

git clone -b ParticleFitter_15_0_X https://github.com/stahlleiton/VertexCompositeAnalysis
git cms-merge-topic stahlleiton:ParticleAnalyzer_CMSSW_15_0_0
scram b -j8

cd VertexCompositeAnalysis/VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2024_DiKa_ParticleAnalyzer_cfg.py
