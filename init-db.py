#!/usr/bin/env python

# initialises a database from genbank references, the names of the genbank
# files are used as the names of the segments

import sys
import os.path
from collections import namedtuple
import sqlite3
from Bio import SeqIO

Segment = namedtuple("Segment", ['name','genbank'])

if len(sys.argv) < 3:
    print("Usage: {:s} database_filename genbank_references ...".format(sys.argv[0]))
else:
    db_file = sys.argv[1]
    gb_files = sys.argv[2:]

if os.path.isfile(db_file):
    exit("File {:s} exists!".format(db_file))

segments = []
for gb_file in gb_files:
    try:
        input_file = open(gb_file)
        segments.append(Segment(gb_file,SeqIO.read(input_file, 'genbank')))
    except IOError as e:
        exit(e)
    finally:
        input_file.close()

db = sqlite3.connect(db_file)

c = db.cursor()
c.execute("create table chromosomes(name unique, length);")
c.execute("create table gene(name unique, product);")
c.execute("create table cds(chromosome, gene, start integer, end integer);")
db.commit()

for segment in segments:
    # print("{:s} = {:d}bp".format(segment.name, len(segment.genbank.seq)))
    c.execute("insert into chromosomes values(?, ?);", 
              (segment.name, len(segment.genbank.seq)))
    chrid = c.lastrowid
    for feature in segment.genbank.features:
        if feature.type == 'CDS':
            gene = feature.qualifiers['gene'][0]
            product = feature.qualifiers['product'][0]
            c.execute("insert or ignore into gene values(?, ?);", (gene, product))
            c.execute("select rowid from gene where name = ?;", (gene,))
            geneid = c.fetchone()[0]
            location = feature.location
            for part in location.parts:
                # print("{:d},{:d},{:s}".format(part.start+1, part.end, gene))
                c.execute("insert into cds values(?,?,?,?);",
                          (chrid, geneid, part.start+1, part.end))

db.commit()
db.close()
