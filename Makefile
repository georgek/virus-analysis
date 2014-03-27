CFLAGS ?= -Wall -O2 -g

all: bam2db trim

bam2db: bam2db.c
	gcc $(CFLAGS) bam2db.c -o bam2db -lsqlite3

trim: trim.c
	gcc $(CFLAGS) trim.c -o trim

clean:
	rm -f bam2db
