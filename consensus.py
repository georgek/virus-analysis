#!/usr/bin/env python

# get consensus sequences from database

import sys
import os.path
import argparse
from collections import namedtuple
import sqlite3

sql_protein_table = """CREATE TEMP TABLE protein AS
SELECT (position-:start)/3+1 AS pos,
Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,Ile,
Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val,STOP
FROM amino_acids WHERE chromosome = :chrid AND position >= :start AND
position <= :end AND (position-:start)%3 = 0 AND
animal = :animal AND day = :day
ORDER BY pos ASC;"""

sql_protein_table_rev = """CREATE TEMP TABLE protein AS
SELECT (:start-position)/3+1 AS pos,
Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,Ile,
Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val,STOP
FROM amino_acids_rev WHERE chromosome = :chrid AND position >= :start AND
position <= :end AND (:start-position)%3 = 0 AND
animal = :animal AND day = :day
ORDER BY pos ASC;"""

def fasta(name, seq):
    "Formats the name and sequence as a fasta record."
    lines = [seq[i:i+80] for i in range(0, len(seq), 80)]
    return ">{:s}\n{:s}\n".format(name, "\n".join(lines))

def count_to_char(count):
    "Returns character corresponding to max count of form [a,c,g,t]."
    if sum(count) > 0:
        return ['A','C','G','T'][count.index(max(count))]
    else:
        return 'N'

def count_to_char_acid(count):
    "Returns character corresponding to max count of form [Ala,Arg,...,Val,STOP]."
    if sum(count) > 0:
        return ['A','R','N','D','C','Q','E','G',
                'H','I','L','K','M','F','P','S',
                'T','W','Y','V','*'][count.index(max(count))]
    else:
        return 'X'

def nuc_consensus(animalid, day, chrid):
    "Returns a nucleotide consensus."
    c.execute("select A,C,G,T from nucleotides_nd where animal = ? and day = ? and chromosome = ? order by position asc",
              (animalid, day, chrid))
    counts = c.fetchall()
    seq = ""
    for count in counts:
        seq += count_to_char([count["a"],count["c"],count["g"],count["t"]])
    return seq

def prot_consensus(animalid, day, geneid):
    "Returns a protein consensus."
    c.execute("select chromosome, start, end, strand from cds_regions join cds on (cds.id = cds) where gene = ?;", (geneid,))
    cds_regions = c.fetchall()
    carry = 0
    seq = ""
    for cds_region in cds_regions:
        if cds_region["strand"] > 0:
            real_start = cds_region["start"] + carry
            real_end = cds_region["end"]
        else:
            real_start = cds_region["start"]
            real_end = cds_region["end"] - carry
        c.execute(sql_protein_table if cds_region["strand"] > 0 else sql_protein_table_rev,
                  {"start": real_start, "chrid": cds_region["chromosome"], "end": real_end,
                   "animal": animalid, "day": day})
        c.execute("select * from protein;")
        counts = c.fetchall()
        for count in counts:
            seq += count_to_char_acid([count["Ala"],count["Arg"],count["Asn"],
                                       count["Asp"],count["Cys"],count["Gln"],
                                       count["Glu"],count["Gly"],count["His"],
                                       count["Ile"],count["Leu"],count["Lys"],
                                       count["Met"],count["Phe"],count["Pro"],
                                       count["Ser"],count["Thr"],count["Trp"],
                                       count["Tyr"],count["Val"],count["STOP"]])
        carry = (3-(real_end - real_start + 1)%3)%3
        if carry > 0:
            seq = seq[0:-1]+'X'
        c.execute("DROP TABLE protein;")
    return seq

def write_sample_files(output_dir, samples):
    for sample in samples:
        filename = output_dir + "/{:s}-d{:d}.fasta".format(sample.animal, sample.day).replace(" ", "_")
        if os.path.isfile(filename):
            print "Error, file exists ({:s}).".format(filename)
            exit(1)
        out = open(filename, 'w')
        for i in range(len(sample.sequences)):
            out.write(fasta(sample.seqnames[i], sample.sequences[i]))
        out.close()

def write_sequence_files(output_dir, samples):
    for i in range(len(samples[0].sequences)):
        filename = output_dir + "/{:s}.fasta".format(samples[0].seqnames[i]).replace(" ", "_")
        if os.path.isfile(filename):
            print "Error, file exists ({:s}).".format(filename)
            exit(1)
        out = open(filename, 'w')
        for j in range(len(samples)):
            out.write(fasta("{:s}-d{:d}".format(samples[j].animal,samples[j].day),
                            samples[j].sequences[i]))
        out.close()

Sample = namedtuple("Sample", ["animal", "day", "seqnames", "sequences"])

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Get consensus sequences from a database.")
parser.add_argument("database", type=str, help="The database file.")
parser.add_argument("output_dir", type=str, help="Output directory.")

grouping = parser.add_mutually_exclusive_group()
grouping.add_argument("-s", "--seq",
                      action="store_const", dest="grouping", const="seq",
                      help="One file per sequence (default).")
grouping.add_argument("-S", "--sample",
                      action="store_const", dest="grouping", const="sample",
                      help="One file per sample.")
parser.set_defaults(grouping="seq")

seq_type = parser.add_mutually_exclusive_group()
seq_type.add_argument("-n", "--nucleotide",
                      action="store_const", dest="seq_type", const="nuc",
                      help="Make nucleotide sequences (default).")
seq_type.add_argument("-p", "--protein",
                      action="store_const", dest="seq_type", const="prot",
                      help="Make protein sequences.")
parser.set_defaults(seq_type="nuc")

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
    c.execute("select * from chromosomes;")
    chromosomes = c.fetchall()
    c.execute("select * from genes;")
    genes = c.fetchall()
    c.execute("select * from animals;")
    animals = c.fetchall()
    for animal in animals:
        c.execute("select distinct day from nucleotides where animal = ?;", (animal["id"],))
        days = c.fetchall()
        for day in days:
            sample = Sample(animal["name"], day["day"], [], [])
            if args.seq_type == "nuc":
                for chromosome in chromosomes:
                    sample.seqnames.append(chromosome["name"])
                    sample.sequences.append(nuc_consensus(animal["id"], day["day"], chromosome["id"]))
            else:
                for gene in genes:
                    sample.seqnames.append(gene["product"])
                    sample.sequences.append(prot_consensus(animal["id"], day["day"], gene["id"]))
            samples.append(sample)

except sqlite3.DatabaseError as e:
    print "sqlite3 error: {:s} ({:s})".format(e, args.database)
    exit(1)

if args.grouping == "sample":
    write_sample_files(args.output_dir, samples)
else:
    write_sequence_files(args.output_dir, samples)

db.close()
