#!/bin/bash

echo "Submit the files at "
echo PWD: $PWD

root -l -b -q "jettrackcorrelation_analyzer.C(\"$1\",$2)"
