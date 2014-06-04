#!/usr/bin/env python

# prints fastq reads between given lengths

import sys
import os

if len(sys.argv) < 3:
    print "Usage: {:s} min_length max_length [filename]".format(sys.argv[0])
    exit(1)

min_length = int(sys.argv[1])
max_length = int(sys.argv[2])
isfile = False
if len(sys.argv) > 3:
    isfile = True
    input_file = open(sys.argv[3])
else:
    input_file = sys.stdin

state = 0
lineno = 1
try:
    for line in input_file:
        if state == 0:
            if line[0] != '@':
                exit("Malformed file, line {:d}, @ expected".format(lineno))
            topheader = line
            state = 1
        elif state == 1:
            sequence = line
            state = 2
        elif state == 2:
            if line[0] != '+':
                exit("Malformed file, line {:d}, + expected".format(lineno))
            midheader = line
            state = 3
        elif state == 3:
            if len(line) != len(sequence):
                exit("Malformed file, line {:d} wrong length".format(lineno))
            quality = line
            seq_length = len(sequence) - 1
            if seq_length >= min_length and seq_length <= max_length:
                sys.stdout.write(topheader)
                sys.stdout.write(sequence)
                sys.stdout.write(midheader)
                sys.stdout.write(quality)
            state = 0
        lineno += 1
except IOError as e:
    exit()
finally:
    if isfile:
        input_file.close()
