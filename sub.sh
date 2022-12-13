#!/bin/bash

echo "Setup CMSSW (Import ROOT version)"
cd /afs/cern.ch/user/d/ddesouza/CMSSW_12_5_0/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/d/ddesouza/news/pPbskims
echo "Submit skim jobs at "
echo PWD: $PWD

root -l -b -q "jettrackcorrelation_analyzer.C(\"$1\",$2)"
