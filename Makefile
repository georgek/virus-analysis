CFLAGS ?= -Wall -O2 -g -I$(HOME)/include -L$(HOME)/lib

all: bam2db trim myccs

bam2db: bam2db.c
	gcc $(CFLAGS) bam2db.c -o bam2db -lsqlite3

trim: trim.c
	gcc $(CFLAGS) trim.c -o trim

myccs: myccs.c
	gcc $(CFLAGS) -pthread myccs.c -o myccs -lbam -lz

clean:
	rm -f bam2db trim myccs
