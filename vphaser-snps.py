#!/usr/bin/env python

import sys

snpscount = {}
snpsfreq = {}

for arg in sys.argv[1:]:
    lines = 0
    file = open(arg, 'r')
    filesnps = {}
    for line in file:
        if line[0] == '#': continue
        lines = lines + 1
        if line.split()[4] != "snp": continue
        pos = int(line.split()[0])
        if pos in filesnps: continue
        filesnps[pos] = True
        if pos in snpscount:
            snpscount[pos] = snpscount[pos] + 1
        else:
            snpscount[pos] = 1
        freq = float(line.split()[5])
        if pos in snpsfreq:
            snpsfreq[pos].append(freq)
        else:
            snpsfreq[pos] = [freq]

    file.close()

for pos in sorted(snpscount, key=snpscount.get, reverse=True):
    meanfreq = sum(snpsfreq[pos])/len(snpsfreq[pos])
    maxfreq = max(snpsfreq[pos])
    minfreq = min(snpsfreq[pos])
    print "{:5d}: {:2d} (min: {:3.3f}, max: {:3.3f}, mean: {:3.3f})".format(
        pos, snpscount[pos], minfreq, maxfreq, meanfreq)

