#!/usr/bin/env python3

import os.path

os.system('rm -r cond/*')
os.system('rm -r out*.sub')
os.system('rm -r *.cc')
os.system('rm files_input/pPb_8160/DATA_PAEGJet/pgoing/*part*.txt')
os.system('rm files_input/pPb_8160/DATA_PAEGJet/Pbgoing/*part*.txt')
os.system('rm files_input/pPb_8160/DATA_MB/pgoing/*part*.txt')
os.system('rm files_input/pPb_8160/DATA_MB/Pbgoing/*part*.txt')
os.system('rm files_input/pPb_8160/DATA_HM185/pgoing/*part*.txt')
os.system('rm files_input/pPb_8160/DATA_HM185/Pbgoing/*part*.txt')
os.system('rm files_input/pPb_8160/DATA_HM250/pgoing/*part*.txt')
os.system('rm files_input/pPb_8160/DATA_HM250/Pbgoing/*part*.txt')
os.system('mkdir -p cond')
