#!/usr/bin/env python

# makes fastq files from bas.h5 files

import sys
import os
from pbcore.io import BasH5Reader, FastqWriter

if len(sys.argv) != 3:
    print "Usage: {:s} bas.h5_file output_prefix".format(sys.argv[0])
    exit(1)

input_filename = sys.argv[1]
output_prefix = sys.argv[2]

bas = BasH5Reader(input_filename)

filenames = {}
writers = {}
filenames['raw'] = output_prefix + ".fastq"
filenames['subread'] = output_prefix + ".subreads.fastq"
filenames['ccs'] = output_prefix + ".ccs.fastq"
for filetype in filenames:
    if os.path.isfile(filenames[filetype]):
        exit("Error: file {:s} exists!".format(filenames[filetype]))
    else:
        writers[filetype] = FastqWriter(filenames[filetype])

for zmw in bas:
    if len(zmw.read()) > 0:
        writers['raw'].writeRecord(zmw.read().readName,
                                   zmw.read().basecalls(),
                                   zmw.read().QualityValue())

    for subread in zmw.subreads:
        if len(subread) > 0:
            writers['subread'].writeRecord(subread.readName,
                                           subread.basecalls(),
                                           subread.QualityValue())

    if zmw.ccsRead is not None:
        writers['ccs'].writeRecord(zmw.ccsRead.readName,
                                   zmw.ccsRead.basecalls(),
                                   zmw.ccsRead.QualityValue())

for filetype in writers:
    writers[filetype].close()

bas.close()
