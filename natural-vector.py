#!/usr/bin/env python

# calculates "natural vector" of a sequence

import sys
import os.path
import argparse

basetable = {'A':0b0001,
             'C':0b0010,
             'G':0b0100,
             'T':0b1000,
             'U':0b1000}
basetable['R'] = basetable['A'] | basetable['G']
basetable['Y'] = basetable['C'] | basetable['T']
basetable['K'] = basetable['G'] | basetable['T']
basetable['M'] = basetable['A'] | basetable['C']
basetable['S'] = basetable['C'] | basetable['G']
basetable['W'] = basetable['A'] | basetable['T']
basetable['B'] = ~basetable['A'] & 0xF
basetable['D'] = ~basetable['C'] & 0xF
basetable['H'] = ~basetable['G'] & 0xF
basetable['V'] = ~basetable['T'] & 0xF
basetable['N'] = 0b1111

def w(b1,b2):
    return float(bin(b1 & b2).count("1"))/bin(b2).count("1")

isfile = False
if len(sys.argv) > 1:
    isfile = True
    input_file = open(sys.argv[1])
else:
    input_file = sys.stdin

position = 1
sequence = []
for line in input_file:
    if line[0] == '>' or line[0] == '\n':
        continue
    for char in line[:-1]:
        sequence.append(basetable[char])

if isfile:
    input_file.close()

bases = map(lambda(c): basetable[c], ['A','C','G','T'])
counts = {}
means = {}
D = {}
for base in bases:
    counts[base] = 0
    means[base] = 0.0
    D[base] = 0.0

for i in range(len(sequence)):
    for base in bases:
        if base & sequence[i]:
            means[base] *= counts[base]/(counts[base]+w(base,sequence[i]))
            means[base] += (w(base,sequence[i])*(i+1))/(counts[base]+w(base,sequence[i]))
            counts[base] += w(base,sequence[i])

for i in range(len(sequence)):
    if sequence[i] in bases:
        D[sequence[i]] += pow(i + 1 - means[sequence[i]], 2)/(counts[sequence[i]]*len(sequence))

sys.stdout.write("(")
for base in bases:
    sys.stdout.write("{:>9.3f},{:>9.3f},{:>9.3f}".format(counts[base],means[base],D[base]))
    if base != bases[-1]:
        sys.stdout.write(",")
sys.stdout.write(")\n")
