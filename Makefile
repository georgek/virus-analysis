CC = c99
CFLAGS ?= -pedantic -Wall -O2 -g -I$(HOME)/include -L$(HOME)/lib

all: bam2db bam2db2 bam2db3 trim myccs bam-dr bampos2readpos bam2fasta distances

bam2db: bam2db.c
	$(CC) $(CFLAGS) bam2db.c -o bam2db -lsqlite3

bam2db2: bam2db2.c errors.h
	$(CC) $(CFLAGS) bam2db2.c -o bam2db2 -lsqlite3 -lbam -lz

bam2db3: bam2db3.c errors.h
	$(CC) $(CFLAGS) bam2db3.c utils.c -o bam2db3 -lsqlite3 -lbam -lz

trim: trim.c
	$(CC) $(CFLAGS) trim.c -o trim

myccs: myccs.c
	$(CC) $(CFLAGS) -pthread myccs.c -o myccs -lbam -lz

bam-dr:	bam-dr.c
	$(CC) $(CFLAGS) -pthread bam-dr.c -o bam-dr -lbam -lz

bamtest: bamtest.c
	$(CC) $(CFLAGS) -pthread bamtest.c -o bamtest -lbam -lz

bampos2readpos: bampos2readpos.c
	$(CC) $(CFLAGS) -pthread bampos2readpos.c -o bampos2readpos -lbam -lz

bam2fasta: bam2fasta.c
	$(CC) $(CFLAGS) -pthread bam2fasta.c -o bam2fasta -lbam -lz

distances: distances.c
	$(CC) $(CFLAGS) distances.c utils.c -o distances -lsqlite3

clean:
	rm -f bam2db bam2db2 bam2db3 trim myccs bam-dr bamtest bampos2readpos
