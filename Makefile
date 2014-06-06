CFLAGS ?= -Wall -O2 -g -I$(HOME)/include -L$(HOME)/lib

all: bam2db bam2db2 trim myccs bam-dr

bam2db: bam2db.c
	gcc $(CFLAGS) bam2db.c -o bam2db -lsqlite3

bam2db2: bam2db2.c
	gcc $(CFLAGS) bam2db2.c -o bam2db2 -lsqlite3 -lbam -lz

trim: trim.c
	gcc $(CFLAGS) trim.c -o trim

myccs: myccs.c
	gcc $(CFLAGS) -pthread myccs.c -o myccs -lbam -lz

bam-dr:	bam-dr.c
	gcc $(CFLAGS) -pthread bam-dr.c -o bam-dr -lbam -lz

bamtest: bamtest.c
	gcc $(CFLAGS) -pthread bamtest.c -o bamtest -lbam -lz

clean:
	rm -f bam2db trim myccs bam-dr
