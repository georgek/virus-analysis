#!/usr/bin/env python

# generates data for variant graphs

import sys
import os.path
import math
import argparse
from collections import namedtuple
import sqlite3

sql_consensus_table = """
CREATE TEMP TABLE consensus AS
SELECT position, 
CASE MAX(SUM(Af)+SUM(Ar), SUM(Cf)+SUM(Cr), 
         SUM(Tf)+SUM(Tr), SUM(Gf)+SUM(Gr))
   WHEN SUM(Af)+SUM(Ar) THEN 'A'
   WHEN SUM(Cf)+SUM(Cr) THEN 'C'
   WHEN SUM(Gf)+SUM(Gr) THEN 'G'
   WHEN SUM(Tf)+SUM(Tr) THEN 'T'
END nuc, chromosome
FROM nucleotides
GROUP BY chromosome, position;
"""

Sample = namedtuple("Sample", ["animal", "day", "seqnames", "sequences"])

bases = ['A','C','G','T']
def argmax(list):
    maxpos = 0
    maxval = list[0]
    for i in range(1,len(list)):
        if list[i] > maxval:
            maxpos = i
            maxval = list[i]
    return maxpos

prominence_thresh = 4
def filter_bias(pos):
    fwd_max = max(pos["Af"],pos["Cf"],pos["Gf"],pos["Tf"])
    rev_max = max(pos["Ar"],pos["Cr"],pos["Gr"],pos["Tr"])
    fwd_sum = pos["Af"]+pos["Cf"]+pos["Gf"]+pos["Tf"]
    rev_sum = pos["Ar"]+pos["Cr"]+pos["Gr"]+pos["Tr"]
    fwd_prominence = float(fwd_max) / fwd_sum if fwd_sum > prominence_thresh else 0
    rev_prominence = float(rev_max) / rev_sum if rev_sum > prominence_thresh else 0
    if fwd_prominence > rev_prominence:
        new_pos = [pos["Af"],pos["Cf"],pos["Gf"],pos["Tf"]]
        ratio = float(rev_sum) / fwd_sum if fwd_sum > 0 else 0
    else:
        new_pos = [pos["Ar"],pos["Cr"],pos["Gr"],pos["Tr"]]
        ratio = float(fwd_sum) / rev_sum if rev_sum > 0 else 0
    return map(lambda n: n * (1 + ratio), new_pos)

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Make distance matrices between each sample from database.")
parser.add_argument("database", type=str, help="The database file.")
parser.add_argument("output_dir", type=str, help="Output directory.")
parser.add_argument("-d", "--delimiters", type=str, dest="delimiter", 
                    default=",", help="Use as delimiter (default ,).")

args = parser.parse_args()
# ----- end command line parsing -----

if os.path.isfile(args.database):
    db = sqlite3.connect(args.database)
else:
    parser.print_usage()
    print("database file must exist!")
    exit(1)

db.row_factory = sqlite3.Row
c = db.cursor()

samples = []
try:
    c.execute("SELECT * FROM chromosomes ORDER BY id ASC;")
    chromosomes = c.fetchall()
    c.execute("SELECT * FROM animals ORDER BY id ASC;")
    animals = c.fetchall()
    for animal in animals:
        c.execute("SELECT DISTINCT day FROM nucleotides WHERE animal = ? ORDER BY day ASC;",
                  (animal["id"],))
        days = c.fetchall()
        for day in days:
            sample = Sample(animal["name"], int(day["day"]), [], [])
            for chromosome in chromosomes:
                sample.seqnames.append(chromosome["name"])
                c.execute("""SELECT position,Af,Ar,Cf,Cr,Gf,Gr,Tf,Tr,D FROM nucleotides
                             WHERE animal = ? AND day = ? AND chromosome = ?
                             ORDER BY position ASC;""",
                          (animal["id"], day["day"], chromosome["id"]))
                sample.sequences.append(c.fetchall())
            samples.append(sample)
except sqlite3.DatabaseError as e:
    print "sqlite3 error: {:s} ({:s})".format(e, args.database)
    exit(1)

delim = args.delimiter.decode("string_escape")
c.execute(sql_consensus_table)
for i in range(len(chromosomes)):
    c.execute("SELECT position,nuc FROM consensus WHERE chromosome = ?;", (i+1,))
    consensus = c.fetchall()
    for sample in samples:
        filename = args.output_dir + "/{:s}-d{:d}-{:s}".format(
            sample.animal,sample.day,chromosomes[i]["name"])
        if os.path.isfile(filename):
            print "Error, file exists ({:s}).".format(filename)
            exit(1)
        out = open(filename, 'w')
        out.write('position{0:s}nucleotide{0:s}amount\n'.format(args.delimiter))
        for cons,pos in zip(consensus,sample.sequences[i]):
            assert(cons["position"] == pos["position"])
            filt = filter_bias(pos)
            out.write(str(cons["position"]))
            out.write(args.delimiter)
            if cons["nuc"] == 'A':
                altb = [bases[1],bases[2],bases[3]]
                altv = [filt[1],filt[2],filt[3]]
            elif cons["nuc"] == 'C':
                altb = [bases[0],bases[2],bases[3]]
                altv = [filt[0],filt[2],filt[3]]
            elif cons["nuc"] == 'G':
                altb = [bases[0],bases[1],bases[3]]
                altv = [filt[0],filt[1],filt[3]]
            elif cons["nuc"] == 'T':
                altb = [bases[0],bases[1],bases[2]]
                altv = [filt[0],filt[1],filt[2]]
            alt = argmax(altv)
            out.write(altb[alt])
            out.write(args.delimiter)
            if sum(filt) > 0:
                out.write(str(altv[alt]/sum(filt)))
            else:
                out.write('0')
            out.write('\n')
        out.write('\n')
        out.close()

db.close()
