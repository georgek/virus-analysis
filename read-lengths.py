#!/usr/bin/env python

import sys

lengths = []
state = 0
lineno = 1

if len(sys.argv) < 1:
    print "Usage: {:s} [fastq_file]".format(sys.argv[0])
    exit(1)

isfile = False
if len(sys.argv) > 1:
    file = open(sys.argv[1])
    isfile = True
else:
    file = sys.stdin

for line in file:
    if state == 0:
        if line[0] != '@':
            exit("Malformed file, line {:d}".format(lineno))
        state = 1
    elif state == 1:
        lengths.append(len(line) - 1)
        state = 2
    elif state == 2:
        if line[0] != '+':
            exit("Malformed file, line {:d}".format(lineno))
        state = 3
    elif state == 3:
        if (len(line) - 1) != lengths[-1]:
            exit("Malformed file, line {:d}".format(lineno))
        state = 0
    lineno += 1

# print "Min: {:.1f}, max: {:.1f}, mean: {:.1f}".format(
#     min(lengths), max(lengths), sum(lengths)/len(lengths))

if isfile:
    file.close()

for length in lengths:
    print length

