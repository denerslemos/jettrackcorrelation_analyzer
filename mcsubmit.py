#!/usr/bin/env python3

'''Add important stuff'''
import os.path
import optparse

outputfolder = "/eos/user/d/ddesouza/Dijethistos"

os.system("python3 Clean_pPb.py")

''' Inputs for the skim code '''
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--input', dest='input', help='add folder on: '+str(outputfolder), default='fwdbkw', type='string') # MC to be implemented tomorrow
parser.add_option('-s', '--side', dest='side', help='choose the side: pgoing or Pbgoing', default='pgoing', type='string')
(opt, args) = parser.parse_args()
inPut = opt.input
sideFiles = opt.side

#os.system("mkdir -p "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles))

if sideFiles == "pgoing" and inPut == "EPOSPYTHIA":
    os.system("condor_submit MC_embedded_Pbgoing_submission_nopthatcut.sub")

if sideFiles == "Pbgoing" and inPut == "EPOSPYTHIA":
    os.system("condor_submit MC_embedded_pgoing_submission_nopthatcut.sub")

if sideFiles == "pgoing" and inPut == "PYTHIA":
    os.system("condor_submit MC_unembedded_Pbgoing_submission_nopthatcut.sub")

if sideFiles == "Pbgoing" and inPut == "PYTHIA":
    os.system("condor_submit MC_unembedded_pgoing_submission_nopthatcut.sub")
