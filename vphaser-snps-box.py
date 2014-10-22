#!/usr/bin/env python

import sys
import os.path
import operator
from sets import Set

usage = "Usage: {:s} sample_sheet genbank_ref_file segn vphaser_dir"

if len(sys.argv) < 5:
    print(usage.format(os.path.basename(sys.argv[0])))
    exit(1)
else:
    ss_file = sys.argv[1]
    gb_file = sys.argv[2]
    segn =    sys.argv[3]
    vp_dir =  sys.argv[4]

segn = int(segn)

try:
    sample_sheet = open(ss_file)
except IOError as e:
    exit(e)

samples = []
hosts = Set()
for line in sample_sheet:
    split = line.split(",")
    samples.append((split[0],split[1]))
    hosts.add(split[0])
sample_sheet.close()

try:
    genbank_file = open(gb_file)
except IOError as e:
    exit(e)

segments = []
for line in genbank_file:
    split = line.split(",")
    segments.append(split[0])
genbank_file.close()

snpcounts = {}

for i in range(len(samples)):
    sample_filename = vp_dir + "/" + samples[i][0] + "-d" + samples[i][1] + "/" + segments[segn] + ".fdr.var.txt"
    if not os.path.isfile(sample_filename): continue
    vp_file = open(sample_filename)
    for line in vp_file:
        if line[0] == '#': continue
        if line.split()[4] != "snp": continue
        pos = int(line.split()[0])
        freq = float(line.split()[5])
        if pos in snpcounts:
            snpcounts[pos][i] += freq
        else:
            snpcounts[pos] = [0] * len(samples)
            snpcounts[pos][i] += freq
    vp_file.close()

onehost = 0
onehostndays = 0
allhosts = 0
for pos, snpcount in snpcounts.iteritems():
    poshosts = {}
    for i in range(len(snpcount)):
        if snpcount[i] > 0:
            if samples[i][0] in poshosts:
                poshosts[samples[i][0]] += 1
            else:
                poshosts[samples[i][0]] = 1
    if len(poshosts) == 1 and poshosts.values()[0] > 1:
        onehostndays += 1
    if len(poshosts) == 1:
        onehost += 1
    if len(poshosts) == len(hosts):
        allhosts += 1

# print "one host:", onehost
# print "all hosts:", allhosts
# print "overall:", len(snpcounts)

# print "| ", onehost, " | ", onehostndays, " | ", allhosts, " | ", len(snpcounts), " |"

sigcounts = []
for pos, snpcount in snpcounts.iteritems():
    if snpcount.count(0) < len(samples) - 2:
        sigcounts.append((pos,snpcount))

sigcounts = sorted(sigcounts,key=lambda(c):map(lambda(x):x>0, c[1]))
for i in range(len(samples)):
    print("{:s}-d{:s}".format(samples[i][0],samples[i][1])),
print
for pos,count in sigcounts:
    print(pos),
    for c in count:
        print(c),
    print
