#! /usr/bin/python
#
# This script

# Usage:    bu2s_utils/bu2s_2h5.py ~/Desktop/fst1.h5 . FSTtimeSeries.txt.bz2
# ~/flaxmans/bu2s/bu2s_utils/bu2s_2h555.py ~/flaxmans/bu2s/runs/sm.txt /media/schimar/schimar2/bu2s_h5/sm/dXY_sm.h5 . dXY



import h5py
from sys import argv
import os
import pandas as pd
import bz2
#from re import findall


h5fileName = argv[2]
direc = argv[3]
param = argv[4]


skiprows = None
header = None
sep= ' '
eng = 'c'    # or 'python'

#
if param == 'Fst':
    parameter = 'FSTtimeSeries.txt.bz2'
elif param == 'Fix':
    parameter = 'FixationLog.txt.bz2'
elif param == 'NeutAF':
    parameter = 'NeutralAlleleFrequencies.txt.bz2'
elif param == 'SelAF':
    parameter = 'SelectedAlleleFrequencies.txt.bz2'
elif param == 'dXY':
    parameter = 'DXYtimeSeries.csv'
    skiprows = [0]
    sep = ','
elif param == 'NeutAFspect':
    parameter = 'NeutralAlleleFrequencySpectrum.csv'
    skiprows = [0]
    sep = ','
elif param == 'SelAFspect':
    parameter = 'SelectedAlleleFrequencySpectrum.csv'
    skiprows = [0]
    sep = ','
elif param == 'LDneutAvg':
    parameter = 'LDneutSitesAvg.txt.bz2'
elif param == 'LDselAvg':
    parameter = 'LDselSitesAvg.txt.bz2'
elif param == 'effMig':
    parameter = 'EffectiveMigrationRates.txt.bz2'
elif param == "map":
    parameter = 'Map.m.bz2'
elif param == 'afts':
    parameter = 'AlleleFreqTimeSeries.txt.bz2'

    # NOTE:
    # CHECK COL NUMBERS FOR THE OTHER FILES!
    # redo the dtypes for columns, (e.g. the number of generations is a float...)   (?)


#dirs = [d for d in os.listdir(direc) if os.path.isdir(os.path.join(direc, d))]   # this one takes forever... let's just index then:

lsdir = os.listdir(direc)#[2:]#1002]

# loop through runs to be added
runLst = list()
with open(argv[1], 'rb') as file:
    for i, line in enumerate(file):
        if line[0] == '#':
            continue
        else:
            runLst.append(line.split('\n')[0])
        #print i, runLst[i]
    file.close()




# go through folders and write hdf5 file
with h5py.File(h5fileName, 'w') as outFile:
    runs = outFile.create_group("runs")
    for i, runDir in enumerate(runLst):
        files = os.listdir(runDir)   # ls within respective run directory
        if runDir in lsdir:
            if runDir in runs:
                continue
            else:
                runs.create_group(runDir)
                if os.stat(runDir + '/' + parameter).st_size != 14:   ## (an empty FixationLog.txt.bz2 has size of 14 bytes)
                    tmpPath = str(runDir + '/' + parameter)
                    #
                    tmpFile =  pd.read_table(tmpPath, sep= sep, header= header, engine= eng, skiprows= skiprows, na_values= '-NAN')
                    #
                    runs[runDir].create_dataset(param, data= tmpFile.transpose()) # NOTE: transpose here
                else:
                    continue


    # here or after 'with..': create another group for the params.txt file (and write it in there)
    outFile.close()


