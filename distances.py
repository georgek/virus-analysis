#!/usr/bin/env python

# get consensus sequences from database

import sys
import os.path
import math
import argparse
from collections import namedtuple
import sqlite3

Sample = namedtuple("Sample", ["animal", "day", "seqnames", "sequences"])

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

def euclidean_distance(pos1, pos2):
    return math.sqrt(sum(map(lambda x,y: pow((x-y), 2), pos1, pos2)))

position_distance_threshold = 100
def position_distance(pos1, pos2):
    sum1 = sum(pos1)
    sum2 = sum(pos2)
    if sum1 > position_distance_threshold and sum2 > position_distance_threshold:
        norm1 = map(lambda n: float(n)/sum1 if sum1 > 0 else 0, pos1)
        norm2 = map(lambda n: float(n)/sum2 if sum2 > 0 else 0, pos2)
        return euclidean_distance(norm1, norm2) / math.sqrt(2)
    else:
        return 0

def sequence_distance(seq1, seq2):
    filt1 = map(filter_bias, seq1)
    filt2 = map(filter_bias, seq2)
    dist = sum(map(position_distance, filt1, filt2)) / min(len(filt1), len(filt2))
    return (-0.74 * math.log(1 - (4.0/3 * dist))) + 0

def cat(seqs):
    return reduce(lambda s1,s2: s1 + s2, seqs)

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
    c.execute("select * from chromosomes order by id asc;")
    chromosomes = c.fetchall()
    c.execute("select * from animals order by id asc;")
    animals = c.fetchall()
    for animal in animals:
        c.execute("select distinct day from nucleotides where animal = ? order by day asc;",
                  str(animal["id"]))
        days = c.fetchall()
        for day in days:
            sample = Sample(animal["name"], int(day["day"]), [], [])
            for chromosome in chromosomes:
                sample.seqnames.append(chromosome["name"])
                c.execute("""select Af,Ar,Cf,Cr,Gf,Gr,Tf,Tr,D from nucleotides
                             where animal = ? and day = ? and chromosome = ?
                             order by position asc;""",
                          (animal["id"], day["day"], chromosome["id"]))
                sample.sequences.append(c.fetchall())
            samples.append(sample)
except sqlite3.DatabaseError as e:
    print "sqlite3 error: {:s} ({:s})".format(e, args.database)
    exit(1)

delim = args.delimiter.decode("string_escape")
for i in range(len(chromosomes)):
    filename = args.output_dir + "/{:s}.dist".format(chromosomes[i]["name"])
    if os.path.isfile(filename):
        print "Error, file exists ({:s}).".format(filename)
        exit(1)
    out = open(filename, 'w')
    out.write(delim.join(map(lambda s: "{:s}-d{:d}".format(s.animal, s.day), samples)))
    out.write('\n')
    for sample in samples:
        out.write(delim.join(map(lambda s: str(sequence_distance(sample.sequences[i],
                                                                 s.sequences[i])),
                                 samples)))
        out.write('\n')
    out.close()

filename = args.output_dir + "/{:s}.dist".format("all")
if os.path.isfile(filename):
    print "Error, file exists ({:s}).".format(filename)
    exit(1)
out = open(filename, 'w')
out.write(delim.join(map(lambda s: "{:s}-d{:d}".format(s.animal, s.day), samples)))
out.write('\n')
for sample in samples:
    out.write(delim.join(map(lambda s: str(sequence_distance(cat(sample.sequences),
                                                             cat(s.sequences))),
                             samples)))
    out.write('\n')
out.close()

db.close()
