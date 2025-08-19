# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_15_0_10_patch1

cd CMSSW_15_0_10_patch1/src

cmsenv

git cms-addpkg DataFormats/PatCandidates
git remote add cmssw git@github.com:stahlleiton/cmssw.git
git fetch cmssw ParticleAnalyzer_CMSSW_15_0_0
git cherry-pick 91562810cee10bac976afe057a878addc088493e

git clone git@github.com:davidlw/VertexCompositeAnalysis.git -b ParticleFitter_15_0_X

scram b -j8

cd VertexCompositeAnalysis/VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2023_Phi_ParticleAnalyzer_MiniAOD_cfg.py
