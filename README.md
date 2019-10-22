# VertexCompositeAnalysis

Example of setting up and running dimuon tree

cmsrel CMSSW_10_3_3_patch1
cd CMSSW_10_3_3_patch1/src
cmsenv
git clone -b 10_3_X https://github.com/davidlw/VertexCompositeAnalysis
cd VertexCompositeAnalysis
scram b -j8
cd VertexCompositeProducer/test
cmsRun PbPbSkimAndTree2018_DiMuContBoth_mc_cfg.py 
