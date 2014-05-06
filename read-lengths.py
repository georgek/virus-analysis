#!/usr/bin/env python

import sys

lengths = []
state = 0
lineno = 1

if len(sys.argv) < 2:
    print "Usage: {:s} fastq_file".format(sys.argv[0])
    exit(1)

file = open(sys.argv[1])
for line in file:
    if state == 0:
        if line[0] != '@':
            exit("Malformed file, line {:d}".format(lineno))
        state = 1
    elif state == 1:
        lengths.append(len(line))
        state = 2
    elif state == 2:
        if line[0] != '+':
            exit("Malformed file, line {:d}".format(lineno))
        state = 3
    elif state == 3:
        if len(line) != lengths[-1]:
            exit("Malformed file, line {:d}".format(lineno))
        state = 0
    lineno += 1

print "Min: {:.1f}, max: {:.1f}, mean: {:.1f}".format(
    min(lengths), max(lengths), sum(lengths)/len(lengths))

