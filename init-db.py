#!/usr/bin/env python

# initialises a database from genbank references, the names of the genbank
# files are used as the names of the segments

import sys
import os.path
from collections import namedtuple
import sqlite3
from Bio import SeqIO

usage = "Usage: {:s} database_filename sample_sheet genbank_references..."

help = """A sqlite3 database will be created called database_filename. sample_sheet
should be a csv containing the animal names and days and read files.
genbank_references should be references in genbank format, the filename of
each is the name of the chromosome/segment."""

Segment = namedtuple("Segment", ['name','genbank'])

if len(sys.argv) < 4:
    print(usage.format(os.path.basename(sys.argv[0])))
    exit(1)
else:
    db_file = sys.argv[1]
    ss_file = sys.argv[2]
    gb_files = sys.argv[3:]

if os.path.isfile(db_file):
    exit("File {:s} exists!".format(db_file))

try:
    sample_sheet = open(ss_file)
except IOError as e:
    exit(e)

segments = []
for gb_file in gb_files:
    try:
        genbank_file = open(gb_file)
        segments.append(Segment(gb_file, SeqIO.read(genbank_file, 'genbank')))
    except IOError as e:
        exit(e)
    finally:
        genbank_file.close()

db = sqlite3.connect(db_file)

c = db.cursor()
c.execute("create table chromosomes(name unique, length integer);")
c.execute("create table genes(name unique, product);")
c.execute("create table cds(chromosome integer, gene integer, start integer, end integer);")
c.execute("create table animals(name, ndays integer);")
c.execute("create table read_data(animal integer, day integer, nreads integer, paired integer);")
db.commit()

animal_ndays = {}
ss_line_no = 1
for line in sample_sheet:
    columns = line.split(',')
    if len(columns) not in [3,4]:
        exit("Sample sheet malformed (line %d)." % ss_line_no)
    ss_line_no += 1
    animal = columns[0]
    if animal in animal_ndays:
        animal_ndays[animal] += 1
    else:
        animal_ndays[animal] = 1
sample_sheet.close()

animalids = {}
for animal in animal_ndays:
    c.execute("insert into animals values(?,?);", 
              (animal, animal_ndays[animal]))
    animalids[animal] = c.lastrowid

try:
    sample_sheet = open(ss_file)
except IOError as e:
    exit(e)
for line in sample_sheet:
    columns = line.split(',')
    animal = columns[0].strip()
    day = columns[1].strip()
    nreads = 0
    paired = 0
    read_file = columns[2].strip()
    try:
        reads = open(read_file)
        print("Counting reads in %s..." % read_file)
        for read in reads:
            nreads += 1
        reads.close()
    except IOError as e:
        exit(e)
    if len(columns) > 3:
        nreads2 = 0
        paired = 1
        read_file = columns[3].strip()
        try:
            reads = open(read_file)
            print("Counting reads in %s..." % read_file)
            for read in reads:
                nreads2 += 1
            reads.close()
        except IOError as e:
            exit(e)
        if nreads != nreads2:
            exit("Paired end reads don't match (%s day %s)." % (animal, day))
    c.execute("insert into read_data values(?,?,?,?)",
              (animalids[animal], day, nreads/4, paired))
sample_sheet.close()

genes = {}
for segment in segments:
    # print("{:s} = {:d}bp".format(segment.name, len(segment.genbank.seq)))
    c.execute("insert into chromosomes values(?, ?);", 
              (segment.name, len(segment.genbank.seq)))
    chrid = c.lastrowid
    for feature in segment.genbank.features:
        if feature.type == 'CDS':
            gene = feature.qualifiers['gene'][0]
            product = feature.qualifiers['product'][0]
            if gene in genes:
                geneid = genes[gene]
            else:
                c.execute("insert into genes values(?, ?);", (gene, product))
                geneid = c.lastrowid
                genes[gene] = geneid
            location = feature.location
            for part in location.parts:
                # print("{:d},{:d},{:s}".format(part.start+1, part.end, gene))
                c.execute("insert into cds values(?,?,?,?);",
                          (chrid, geneid, part.start+1, part.end))

db.commit()
db.close()
