CC = gcc --std=c99
CPPFLAGS = -I$(HOME)/include
CFLAGS = -pedantic -Wall -O2 -g
LDFLAGS = -L$(HOME)/lib
LDLIBS = -pthread -lsqlite3 -lbam -lz

PROGS = bam2db3 trim myccs bam-dr bampos2readpos bam2fasta distances

all: $(PROGS)

bam2db3: LDLIBS = -pthread -lsqlite3 -lhts -lz
bam2db3: bam2db3.o utils.o

trim: trim.o

myccs: myccs.o

bam-dr:	bam-dr.o

bamtest: bamtest.o

bampos2readpos: bampos2readpos.o

bam2fasta: bam2fasta.o

distances: distances.o utils.o

clean:
	rm -f $(PROGS)
