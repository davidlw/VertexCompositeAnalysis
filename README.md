# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_10_6_4_patch1

cd CMSSW_10_6_4_patch1/src

cmsenv

git clone -b ParticleFitter_10_6_X https://github.com/davidlw/VertexCompositeAnalysis

cd VertexCompositeAnalysis

scram b -j8

cd VertexCompositeProducer/test

## V0 reconstruction inside jets in pp using MINIAOD
cmsRun ppRun2UL_V0Both_MiniAOD_cfg.py 
