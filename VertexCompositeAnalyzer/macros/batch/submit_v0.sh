#!/bin/sh

echo "setup cmssw"
cd ~davidlw/CMSSW/CMSSW_10_6_4_patch1/src/
eval `scramv1 runtime -sh`
cd ~davidlw/CMSSW/CMSSW_10_6_4_patch1/src/VertexCompositeAnalysis/VertexCompositeAnalyzer/macros/batch
echo PWD: $PWD

root -b ../jetv0.C+\(\"list_v0_$1\"\)
