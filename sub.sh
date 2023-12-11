#!/bin/bash

echo "Setup CMSSW (Import ROOT version)"
cd /afs/cern.ch/work/d/ddesouza/UIC/pPbMultAna/CMSSW_13_0_5/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/d/ddesouza/UIC/pPbMultAna/CMSSW_13_0_5/src/dijetMult/pgoing
echo "Submit skim jobs at "
echo PWD: $PWD

root -l -b -q "jettrackcorrelation_analyzer.C(\"$1\",\"$2\",$3,$4,$5)"
