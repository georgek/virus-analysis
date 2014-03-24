CFLAGS ?= -Wall -O2 -g

bam2db:bam2db.c
	gcc $(CFLAGS) bam2db.c -o bam2db -lsqlite3

clean:
	rm -f bam2db
