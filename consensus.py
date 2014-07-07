#!/usr/bin/env python

# get consensus sequences from database

import sys
import os.path
import argparse
import sqlite3

def fasta(name, seq):
    "Formats the name and sequence as a fasta record."
    lines = [seq[i:i+80] for i in range(0, len(seq), 80)]
    return ">{:s}\n{:s}".format(name, "\n".join(lines))

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

if os.path.isfile(args.database):
    db = sqlite3.connect(args.database)
else:
    parser.print_usage()
    print("database file must exist!")
    exit(1)

db.row_factory = sqlite3.Row
c = db.cursor()

try:
    c.execute("select * from animals;")
    animals = c.fetchall()
    samples = []
    for animal in animals:
        for row in c.execute("select distinct day from nucleotides where animal = ?;", str(animal["id"])):
            print animal['id'], row['day']
        
except sqlite3.DatabaseError as e:
    print "sqlite3 error: {:s} ({:s})".format(e, args.database)
    exit(1)



db.close()
