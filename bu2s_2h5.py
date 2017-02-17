#! /usr/bin/python
#
# This script

# Usage:


import h5py
from sys import argv
from os import listdir
import pandas as pd
import bz2
#from re import findall


h5fileName = argv[1]
ls = listdir(argv[2])


with h5py.File(h5fileName, 'w') as outFile:
    runs = outFile.create_group("runs")
    for i, runDir in enumerate(ls):
        files = listdir(runDir)
        for f in files:
            if ls[i] in runs:
                continue
            else:
                runs.create_group(ls[i])
                #tempFst = pd.read_table("FSTtimeSeries.txt.bz2")
                tmpPathFst = str(ls[i] + '/FSTtimeSeries.txt.bz2')
                runs[ls[i]].create_dataset("Fst", data= pd.read_table(tmpPathFst, sep= ' ', header= None)) #, names= ["totalGenerationsElapsed", "locusID", "Fst", "allele_frequencies", "S_MAX1", "S_MAX0", "chromosomeMembership", "MAP", "locType"]))


    # here or after 'with..': create another group for the params.txt file (and write it in there)
    outFile.close()





        #runs = outFile.create_group("runs")

