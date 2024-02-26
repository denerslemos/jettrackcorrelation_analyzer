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

os.system("mkdir -p "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles))

if sideFiles == "pgoing":
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/pgoing/HM250_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM250 -f workday -c 1 -n 100 -s outHM250")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD1_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD1 -f workday -c 1 -n 100 -s outHM185PD1")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD2_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD2 -f workday -c 1 -n 100 -s outHM185PD2")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD3_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD3 -f workday -c 1 -n 100 -s outHM185PD3")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD4_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD4 -f workday -c 1 -n 100 -s outHM185PD4")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD5_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD5 -f workday -c 1 -n 100 -s outHM185PD5")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD6_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD6 -f workday -c 1 -n 100 -s outHM185PD6")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD1_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD1 -f workday -c 1 -n 100 -s outMBPD1")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD2_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD2 -f workday -c 1 -n 100 -s outMBPD2")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD3_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD3 -f workday -c 1 -n 100 -s outMBPD3")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD4_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD4 -f workday -c 1 -n 100 -s outMBPD4")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD5_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD5 -f workday -c 1 -n 100 -s outMBPD5")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD6_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD6 -f workday -c 1 -n 100 -s outMBPD6")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD7_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD7 -f workday -c 1 -n 100 -s outMBPD7")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD8_pgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD8 -f workday -c 1 -n 100 -s outMBPD8")
	os.system("condor_submit MC_embedded_Pbgoing_submission_nopthatcut.sub")

if sideFiles == "Pbgoing":
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/Pbgoing/HM250_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM250 -f workday -c 1 -n 100 -s outHM250")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD1_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD1 -f workday -c 1 -n 100 -s outHM185PD1")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD2_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD2 -f workday -c 1 -n 100 -s outHM185PD2")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD3_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD3 -f workday -c 1 -n 100 -s outHM185PD3")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD4_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD4 -f workday -c 1 -n 100 -s outHM185PD4")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD5_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD5 -f workday -c 1 -n 100 -s outHM185PD5")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD6_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/HM185PD6 -f workday -c 1 -n 100 -s outHM185PD6")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD1_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD1 -f workday -c 1 -n 100 -s outMBPD1")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD2_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD2 -f workday -c 1 -n 100 -s outMBPD2")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD3_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD3 -f workday -c 1 -n 100 -s outMBPD3")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD4_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD4 -f workday -c 1 -n 100 -s outMBPD4")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD5_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD5 -f workday -c 1 -n 100 -s outMBPD5")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD6_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD6 -f workday -c 1 -n 100 -s outMBPD6")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD7_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD7 -f workday -c 1 -n 100 -s outMBPD7")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD8_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD8 -f workday -c 1 -n 100 -s outMBPD8")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD9_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD9 -f workday -c 1 -n 100 -s outMBPD9")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD10_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD10 -f workday -c 1 -n 100 -s outMBPD10")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD11_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD11 -f workday -c 1 -n 100 -s outMBPD11")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD12_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD12 -f workday -c 1 -n 100 -s outMBPD12")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD13_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD13 -f workday -c 1 -n 100 -s outMBPD13")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD14_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD14 -f workday -c 1 -n 100 -s outMBPD14")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD15_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD15 -f workday -c 1 -n 100 -s outMBPD15")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD16_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD16 -f workday -c 1 -n 100 -s outMBPD16")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD17_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD17 -f workday -c 1 -n 100 -s outMBPD17")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD18_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD18 -f workday -c 1 -n 100 -s outMBPD18")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD19_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD19 -f workday -c 1 -n 100 -s outMBPD19")
	os.system("python3 HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD20_Pbgoing -o "+str(outputfolder)+"/"+str(inPut)+"/"+str(sideFiles)+"/MBPD20 -f workday -c 1 -n 100 -s outMBPD20")
	os.system("condor_submit MC_embedded_pgoing_submission_nopthatcut.sub")