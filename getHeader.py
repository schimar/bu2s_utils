#! /usr/bin/python
#
# This script simply reads a file ("headers.txt") containing header information for a set of files from bu2s and prints said header info for the supplied file name (2nd argument).

# Usage: ./getHeader.py <file_name>


from sys import argv

file = argv[1]

with open("headers.txt", 'r') as hdrFile:
    hdrs = dict()
    for line in hdrFile:
        line = line.split('\n')[0].split(' ')
        hdrs[line[0]] = line[1:len(line)]
        #print hdrs[line[0]]
    hdrFile.close()


print ' '.join(hdrs[file])


