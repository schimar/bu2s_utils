#! /usr/bin/python
#
# This script reads through the JAFSdata folder for a given bu2s run, and returns the data for all respective files, so they can be written to one single file per run.
# NOTE: you have to be in the respective (JAfSdata) folder, as the script calls 'ls' in the folder.


# Usage: ./numPerPatch2df.py > outfile.txt

#from sys import argv
from os import listdir
from re import findall

lsFold = listdir('.')               #argv[1])

patt = 'JAFSnumPerPatch[0-9]+\.R'

for file in lsFold:
    if findall(patt, file):
        match = findall(patt, file)[0]
        nGen = match.split('.')[0].split('Patch')[1]
        with open(match, 'rb') as matchFile:
            for line in matchFile:
                patch0, patch1 = line.split('\n')[0].split('(')[1].split(')')[0].split(',')
                print nGen + '\t' + patch0 + '\t' + patch1
            matchFile.close()








