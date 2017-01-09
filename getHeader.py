#! /usr/bin/python
#
# This script

# reads a table containing header information for the bu2s output (see

# Usage: ./clean_runInfo.py <input-file_name.txt> <new_file_name.txt>


from sys import argv
#import re
#import shutil
#import tempfile

file = argv[1]

with open("headers.txt", 'r') as hdrFile:
    hdrs = dict()
    for line in hdrFile:
        print line

    hdrFile.close()






