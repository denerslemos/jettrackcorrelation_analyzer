#!/usr/bin/env python

'''Add important stuff'''
import os.path
import optparse


''' Inputs for the skim code '''
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--input', dest='input', help='input: datashapes3p0 or datashapes1p6 or dataflow1p3 or datafwdbkw or MCUnemb or MCUnembnocut or MCEmb or MCEmbnocut', default='datashapes1p6', type='string') # MC to be implemented tomorrow
parser.add_option('-s', '--side', dest='side', help='choose the side: pgoing or Pbgoing', default='pgoing', type='string')
(opt, args) = parser.parse_args()
inPut = opt.input
sideFiles = opt.side

if inPut == "MCUnemb" and sideFiles == "pgoing":
	os.system("condor_submit MC_unembedded_pgoing_submission.sub")

if inPut == "MCUnemb" and sideFiles == "Pbgoing":
	os.system("condor_submit MC_unembedded_Pbgoing_submission.sub")
	
if inPut == "MCUnembnocut" and sideFiles == "pgoing":
	os.system("condor_submit MC_unembedded_pgoing_submission_nopthatcut.sub")

if inPut == "MCUnembnocut" and sideFiles == "Pbgoing":
	os.system("condor_submit MC_unembedded_Pbgoing_submission_nopthatcut.sub")

if inPut == "MCEmb" and sideFiles == "pgoing":
	os.system("condor_submit MC_embedded_pgoing_submission.sub")

if inPut == "MCEmb" and sideFiles == "Pbgoing":
	os.system("condor_submit MC_embedded_Pbgoing_submission.sub")
	
if inPut == "MCEmbnocut" and sideFiles == "pgoing":
	os.system("condor_submit MC_embedded_pgoing_submission_nopthatcut.sub")

if inPut == "MCEmbnocut" and sideFiles == "Pbgoing":
	os.system("condor_submit MC_embedded_Pbgoing_submission_nopthatcut.sub")

if inPut == "datashapes1p6" and sideFiles == "pgoing":
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/pgoing/HM250_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/HM250 -f testmatch -c 16 -n 100 -s outHM250")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD1_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/HM185PD1 -f testmatch -c 16 -n 15 -s outHM185PD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD2_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/HM185PD2 -f testmatch -c 16 -n 15 -s outHM185PD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD3_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/HM185PD3 -f testmatch -c 16 -n 15 -s outHM185PD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD4_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/HM185PD4 -f testmatch -c 16 -n 15 -s outHM185PD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD5_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/HM185PD5 -f testmatch -c 16 -n 15 -s outHM185PD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD6_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/HM185PD6 -f testmatch -c 16 -n 15 -s outHM185PD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD1_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/MBPD1 -f testmatch -c 16 -n 100 -s outMBPD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD2_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/MBPD2 -f testmatch -c 16 -n 100 -s outMBPD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD3_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/MBPD3 -f testmatch -c 16 -n 100 -s outMBPD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD4_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/MBPD4 -f testmatch -c 16 -n 100 -s outMBPD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD5_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/MBPD5 -f testmatch -c 16 -n 100 -s outMBPD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD6_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/MBPD6 -f testmatch -c 16 -n 100 -s outMBPD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD7_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/MBPD7 -f testmatch -c 16 -n 100 -s outMBPD7")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD8_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/pgoing/MBPD8 -f testmatch -c 16 -n 100 -s outMBPD8")

if inPut == "datashapes1p6" and sideFiles == "Pbgoing":
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/Pbgoing/HM250_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/HM250 -f testmatch -c 16 -n 10 -s outHM250")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD1_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/HM185PD1 -f testmatch -c 16 -n 100 -s outHM185PD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD2_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/HM185PD2 -f testmatch -c 16 -n 100 -s outHM185PD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD3_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/HM185PD3 -f testmatch -c 16 -n 100 -s outHM185PD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD4_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/HM185PD4 -f testmatch -c 16 -n 100 -s outHM185PD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD5_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/HM185PD5 -f testmatch -c 16 -n 100 -s outHM185PD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD6_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/HM185PD6 -f testmatch -c 16 -n 100 -s outHM185PD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD1_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD1 -f testmatch -c 16 -n 30 -s outMBPD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD2_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD2 -f testmatch -c 16 -n 30 -s outMBPD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD3_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD3 -f testmatch -c 16 -n 30 -s outMBPD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD4_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD4 -f testmatch -c 16 -n 30 -s outMBPD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD5_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD5 -f testmatch -c 16 -n 30 -s outMBPD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD6_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD6 -f testmatch -c 16 -n 30 -s outMBPD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD7_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD7 -f testmatch -c 16 -n 30 -s outMBPD7")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD8_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD8 -f testmatch -c 16 -n 30 -s outMBPD8")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD9_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD9 -f testmatch -c 16 -n 30 -s outMBPD9")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD10_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD10 -f testmatch -c 16 -n 30 -s outMBPD10")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD11_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD11 -f testmatch -c 16 -n 30 -s outMBPD11")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD12_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD12 -f testmatch -c 16 -n 30 -s outMBPD12")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD13_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD13 -f testmatch -c 16 -n 30 -s outMBPD13")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD14_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD14 -f testmatch -c 16 -n 30 -s outMBPD14")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD15_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD15 -f testmatch -c 16 -n 30 -s outMBPD15")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD16_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD16 -f testmatch -c 16 -n 30 -s outMBPD16")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD17_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD17 -f testmatch -c 16 -n 30 -s outMBPD17")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD18_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD18 -f testmatch -c 16 -n 30 -s outMBPD18")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD19_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD19 -f testmatch -c 16 -n 30 -s outMBPD19")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD20_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p6shapes/Pbgoing/MBPD20 -f testmatch -c 16 -n 30 -s outMBPD20")

if inPut == "datashapes3p0" and sideFiles == "pgoing":
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/pgoing/HM250_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/HM250 -f tomorrow -c 4 -n 10 -s outHM250")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD1_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/HM185PD1 -f testmatch -c 16 -n 15 -s outHM185PD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD2_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/HM185PD2 -f testmatch -c 16 -n 15 -s outHM185PD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD3_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/HM185PD3 -f testmatch -c 16 -n 15 -s outHM185PD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD4_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/HM185PD4 -f testmatch -c 16 -n 15 -s outHM185PD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD5_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/HM185PD5 -f testmatch -c 16 -n 15 -s outHM185PD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD6_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/HM185PD6 -f testmatch -c 16 -n 15 -s outHM185PD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD1_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/MBPD1 -f testmatch -c 16 -n 100 -s outMBPD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD2_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/MBPD2 -f testmatch -c 16 -n 100 -s outMBPD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD3_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/MBPD3 -f testmatch -c 16 -n 100 -s outMBPD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD4_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/MBPD4 -f testmatch -c 16 -n 100 -s outMBPD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD5_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/MBPD5 -f testmatch -c 16 -n 100 -s outMBPD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD6_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/MBPD6 -f testmatch -c 16 -n 100 -s outMBPD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD7_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/MBPD7 -f testmatch -c 16 -n 100 -s outMBPD7")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD8_pgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/pgoing/MBPD8 -f testmatch -c 16 -n 100 -s outMBPD8")

if inPut == "datashapes3p0" and sideFiles == "Pbgoing":
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/Pbgoing/HM250_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/HM250 -f testmatch -c 16 -n 10 -s outHM250")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD1_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/HM185PD1 -f testmatch -c 16 -n 100 -s outHM185PD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD2_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/HM185PD2 -f testmatch -c 16 -n 100 -s outHM185PD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD3_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/HM185PD3 -f testmatch -c 16 -n 100 -s outHM185PD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD4_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/HM185PD4 -f testmatch -c 16 -n 100 -s outHM185PD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD5_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/HM185PD5 -f testmatch -c 16 -n 100 -s outHM185PD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD6_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/HM185PD6 -f testmatch -c 16 -n 100 -s outHM185PD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD1_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD1 -f testmatch -c 16 -n 30 -s outMBPD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD2_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD2 -f testmatch -c 16 -n 30 -s outMBPD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD3_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD3 -f testmatch -c 16 -n 30 -s outMBPD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD4_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD4 -f testmatch -c 16 -n 30 -s outMBPD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD5_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD5 -f testmatch -c 16 -n 30 -s outMBPD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD6_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD6 -f testmatch -c 16 -n 30 -s outMBPD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD7_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD7 -f testmatch -c 16 -n 30 -s outMBPD7")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD8_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD8 -f testmatch -c 16 -n 30 -s outMBPD8")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD9_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD9 -f testmatch -c 16 -n 30 -s outMBPD9")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD10_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD10 -f testmatch -c 16 -n 30 -s outMBPD10")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD11_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD11 -f testmatch -c 16 -n 30 -s outMBPD11")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD12_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD12 -f testmatch -c 16 -n 30 -s outMBPD12")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD13_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD13 -f testmatch -c 16 -n 30 -s outMBPD13")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD14_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD14 -f testmatch -c 16 -n 30 -s outMBPD14")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD15_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD15 -f testmatch -c 16 -n 30 -s outMBPD15")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD16_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD16 -f testmatch -c 16 -n 30 -s outMBPD16")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD17_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD17 -f testmatch -c 16 -n 30 -s outMBPD17")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD18_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD18 -f testmatch -c 16 -n 30 -s outMBPD18")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD19_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD19 -f testmatch -c 16 -n 30 -s outMBPD19")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD20_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta3p0shapes/Pbgoing/MBPD20 -f testmatch -c 16 -n 30 -s outMBPD20")

if inPut == "dataflow1p3" and sideFiles == "pgoing":
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/pgoing/HM250_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/HM250 -f testmatch -c 16 -n 10 -s outHM250")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD1_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/HM185PD1 -f testmatch -c 16 -n 100 -s outHM185PD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD2_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/HM185PD2 -f testmatch -c 16 -n 100 -s outHM185PD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD3_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/HM185PD3 -f testmatch -c 16 -n 100 -s outHM185PD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD4_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/HM185PD4 -f testmatch -c 16 -n 100 -s outHM185PD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD5_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/HM185PD5 -f testmatch -c 16 -n 100 -s outHM185PD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD6_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/HM185PD6 -f testmatch -c 16 -n 100 -s outHM185PD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD1_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/MBPD1 -f testmatch -c 16 -n 100 -s outMBPD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD2_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/MBPD2 -f testmatch -c 16 -n 100 -s outMBPD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD3_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/MBPD3 -f testmatch -c 16 -n 100 -s outMBPD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD4_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/MBPD4 -f testmatch -c 16 -n 100 -s outMBPD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD5_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/MBPD5 -f testmatch -c 16 -n 100 -s outMBPD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD6_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/MBPD6 -f testmatch -c 16 -n 100 -s outMBPD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD7_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/MBPD7 -f testmatch -c 16 -n 100 -s outMBPD7")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD8_pgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/pgoing/MBPD8 -f testmatch -c 16 -n 100 -s outMBPD8")

if inPut == "dataflow1p3" and sideFiles == "Pbgoing":
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/Pbgoing/HM250_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/HM250 -f testmatch -c 16 -n 10 -s outHM250")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD1_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/HM185PD1 -f testmatch -c 16 -n 100 -s outHM185PD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD2_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/HM185PD2 -f testmatch -c 16 -n 100 -s outHM185PD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD3_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/HM185PD3 -f testmatch -c 16 -n 100 -s outHM185PD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD4_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/HM185PD4 -f testmatch -c 16 -n 100 -s outHM185PD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD5_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/HM185PD5 -f testmatch -c 16 -n 100 -s outHM185PD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD6_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/HM185PD6 -f testmatch -c 16 -n 100 -s outHM185PD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD1_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD1 -f testmatch -c 16 -n 100 -s outMBPD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD2_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD2 -f testmatch -c 16 -n 100 -s outMBPD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD3_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD3 -f testmatch -c 16 -n 100 -s outMBPD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD4_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD4 -f testmatch -c 16 -n 100 -s outMBPD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD5_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD5 -f testmatch -c 16 -n 100 -s outMBPD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD6_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD6 -f testmatch -c 16 -n 100 -s outMBPD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD7_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD7 -f testmatch -c 16 -n 100 -s outMBPD7")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD8_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD8 -f testmatch -c 16 -n 100 -s outMBPD8")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD9_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD9 -f testmatch -c 16 -n 100 -s outMBPD9")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD10_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD10 -f testmatch -c 16 -n 100 -s outMBPD10")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD11_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD11 -f testmatch -c 16 -n 100 -s outMBPD11")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD12_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD12 -f testmatch -c 16 -n 100 -s outMBPD12")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD13_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD13 -f testmatch -c 16 -n 100 -s outMBPD13")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD14_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD14 -f testmatch -c 16 -n 100 -s outMBPD14")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD15_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD15 -f testmatch -c 16 -n 100 -s outMBPD15")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD16_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD16 -f testmatch -c 16 -n 100 -s outMBPD16")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD17_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD17 -f testmatch -c 16 -n 100 -s outMBPD17")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD18_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD18 -f testmatch -c 16 -n 100 -s outMBPD18")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD19_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD19 -f testmatch -c 16 -n 100 -s outMBPD19")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD20_Pbgoing -o /eos/user/d/ddesouza/JetResults/eta1p3flow/Pbgoing/MBPD20 -f testmatch -c 16 -n 100 -s outMBPD20")

if inPut == "datafwdbkw" and sideFiles == "pgoing":
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/pgoing/HM250_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/HM250 -f workday -c 2 -n 100 -s outHM250")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD1_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/HM185PD1 -f workday -c 2 -n 100 -s outHM185PD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD2_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/HM185PD2 -f workday -c 2 -n 100 -s outHM185PD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD3_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/HM185PD3 -f workday -c 2 -n 100 -s outHM185PD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD4_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/HM185PD4 -f workday -c 2 -n 100 -s outHM185PD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD5_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/HM185PD5 -f workday -c 2 -n 100 -s outHM185PD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/pgoing/HM185_PD6_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/HM185PD6 -f workday -c 2 -n 100 -s outHM185PD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD1_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/MBPD1 -f workday -c 2 -n 100 -s outMBPD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD2_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/MBPD2 -f workday -c 2 -n 100 -s outMBPD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD3_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/MBPD3 -f workday -c 2 -n 100 -s outMBPD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD4_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/MBPD4 -f workday -c 2 -n 100 -s outMBPD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD5_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/MBPD5 -f workday -c 2 -n 100 -s outMBPD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD6_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/MBPD6 -f workday -c 2 -n 100 -s outMBPD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD7_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/MBPD7 -f workday -c 2 -n 100 -s outMBPD7")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/pgoing/MB_PD8_pgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/pgoing/MBPD8 -f workday -c 2 -n 100 -s outMBPD8")

if inPut == "datafwdbkw" and sideFiles == "Pbgoing":
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM250/Pbgoing/HM250_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/HM250 -f workday -c 2 -n 100 -s outHM250")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD1_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/HM185PD1 -f workday -c 2 -n 100 -s outHM185PD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD2_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/HM185PD2 -f workday -c 2 -n 100 -s outHM185PD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD3_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/HM185PD3 -f workday -c 2 -n 100 -s outHM185PD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD4_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/HM185PD4 -f workday -c 2 -n 100 -s outHM185PD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD5_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/HM185PD5 -f workday -c 2 -n 100 -s outHM185PD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_HM185/Pbgoing/HM185_PD6_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/HM185PD6 -f workday -c 2 -n 100 -s outHM185PD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD1_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD1 -f workday -c 2 -n 100 -s outMBPD1")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD2_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD2 -f workday -c 2 -n 100 -s outMBPD2")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD3_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD3 -f workday -c 2 -n 100 -s outMBPD3")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD4_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD4 -f workday -c 2 -n 100 -s outMBPD4")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD5_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD5 -f workday -c 2 -n 100 -s outMBPD5")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD6_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD6 -f workday -c 2 -n 100 -s outMBPD6")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD7_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD7 -f workday -c 2 -n 100 -s outMBPD7")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD8_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD8 -f workday -c 2 -n 100 -s outMBPD8")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD9_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD9 -f workday -c 2 -n 100 -s outMBPD9")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD10_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD10 -f workday -c 2 -n 100 -s outMBPD10")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD11_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD11 -f workday -c 2 -n 100 -s outMBPD11")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD12_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD12 -f workday -c 2 -n 100 -s outMBPD12")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD13_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD13 -f workday -c 2 -n 100 -s outMBPD13")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD14_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD14 -f workday -c 2 -n 100 -s outMBPD14")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD15_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD15 -f workday -c 2 -n 100 -s outMBPD15")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD16_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD16 -f workday -c 2 -n 100 -s outMBPD16")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD17_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD17 -f workday -c 2 -n 100 -s outMBPD17")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD18_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD18 -f workday -c 2 -n 100 -s outMBPD18")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD19_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD19 -f workday -c 2 -n 100 -s outMBPD19")
	os.system("python HTCondor_submit_data.py -i files_input/pPb_8160/DATA_MB/Pbgoing/MB_PD20_Pbgoing -o /eos/user/d/ddesouza/JetResults/fwdbkw/Pbgoing/MBPD20 -f workday -c 2 -n 100 -s outMBPD20")
