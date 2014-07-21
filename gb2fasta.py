#!/usr/bin/env python

# makes fasta files from genbank references

import sys
import os.path
from Bio import SeqIO

def pathname(path):
    name = os.path.dirname(path)
    if name:
        return name
    else:
        return '.'

def fasta(name, seq):
    "Formats the name and sequence as a fasta record."
    lines = [seq[i:i+80] for i in range(0, len(seq), 80)]
    return ">{:s}\n{:s}\n".format(name, "\n".join(lines))


usage = "Usage: {:s} genbank_ref_file output_file"

if len(sys.argv) < 3:
    print(usage.format(os.path.basename(sys.argv[0])))
    exit(1)
else:
    gb_file = sys.argv[1]
    output_file = sys.argv[2]

if os.path.isfile(output_file):
    exit("File {:s} exists!".format(output_file))

try:
    gb = open(gb_file)
    output = open(output_file, 'w')
except IOError as e:
    exit(e)

for line in gb:
    split = line[:-1].split(',')
    genbank = open(pathname(gb_file) + '/' + split[1])
    print genbank
    gbrec = SeqIO.read(genbank, "genbank")
    output.write(fasta(split[0], str(gbrec.seq)))
    genbank.close()

output.close()
gb.close()
