#!/usr/bin/env python

# find positions where sequences do not all agree

import sys

positions = []
state = 0

if len(sys.argv) < 2:
    print "Usage: {:s} fasta_file".format(sys.argv[0])
    exit(1)

file = open(sys.argv[1])
file.readline()
firstline = file.readline()
length = len(firstline)

lineno = 3
for line in file:
    if state == 0:
        if line[0] != '>':
            exit("Malformed file, line {:d}".format(lineno))
        state = 1
    if state == 1:
        if len(line) != length:
            exit("Malformed file, line {:d}".format(lineno))
        else:
            for pos in range(length - 1):
                if firstline[pos] != line[pos] and pos+1 not in positions:
                    positions.append(pos+1)
        state = 0
    lineno += 1

for pos in positions:
    print pos
