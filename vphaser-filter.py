#!/usr/bin/env python

# filters variants below threshold from vphaser file

import sys
import os

if len(sys.argv) < 2:
    print "Usage: {:s} threshold [filename]".format(os.path.basename(sys.argv[0]))
    exit(1)

try:
    threshold = float(sys.argv[1])
except ValueError as e:
    print "Threshold must be a number."
    exit()

isfile = False
if len(sys.argv) > 2:
    isfile = True
    input_file = open(sys.argv[2])
else:
    input_file = sys.stdin

nsnps = 0
nlps = 0

try:
    for line in input_file:
        split = line.split()
        if line[0] == '#':
            if split[1] == "Summary:":
                sys.stdout.write("{:s} {:s} {:s} {:d};\t {:s} {:d}\n"
                                 .format(split[0], split[1], split[2], nsnps,
                                         split[4], nlps))
            else:
                sys.stdout.write(line)
        elif float(split[5]) >= threshold:
            if split[4] == "snp":
                nsnps += 1
            elif split[4] == "lp":
                nlps += 1
            sys.stdout.write(line)
finally:
    input_file.close()

# weird hack to stop pipe errors
try:
    sys.stdout.close()
except:
    pass
try:
    sys.stderr.close()
except:
    pass
