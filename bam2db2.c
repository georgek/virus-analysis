#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sqlite3.h>

#include <bam/sam.h>
#include <bam/bam.h>

static sqlite3 *db = NULL;
static sqlite3_stmt *stmt;
static samfile_t *bamin;
static bam1_t *b;
static char *dbname, *bamfilename, *animal;
static int day;

int animal_id(void* ret_val, int ncols, char **vals, char **names)
{
     if (1 == ncols){
          *(int *)ret_val = strtoul(vals[0], NULL, 10);
          return 0;
     }
     else {
          return 1;
     }
}

int main(int argc, char *argv[])
{
     int dbopencode = 0;
     char *errormessage = NULL;
     uint8_t *dr = NULL;
     uint32_t dr_count = 0;

     b = bam_init1();

     if (argc < 5) {
          fprintf(stderr, "Usage: %s database bamfile animal day\n", argv[0]);
          exit(1);
     }
     else {
          printf("DB name: %s\n", argv[1]);
          dbname = argv[1];
          bamfilename = argv[2];
          animal = argv[3];
          day = strtoul(argv[3], NULL, 10);
     }
     /* try to open db */
     dbopencode = sqlite3_open(argv[1], &db);
     if (SQLITE_OK == dbopencode) {
          printf("OK\n");
     }
     else {
          printf("SQLite3 error: %d\n", dbopencode);
          exit(1);
     }
     /* increases speed of insertion (means if the program crashes the db is
      * left invalid) */
     sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL, &errormessage);

     /* try to open bam file */
     bamin = samopen(bamfilename, "rb", NULL);
     if (!bamin) {
          fprintf(stderr, "Error opening bamfile %s\n", bamfilename);
          exit(1);
     }

     

     while (samread(bamin, b) >= 0) {
          dr = bam_aux_get(b, "DR");
          if (dr) {
               dr_count += bam_aux2i(dr);
          }
     }

     samclose(bamin);
     bam_destroy1(b);

     sqlite3_close(db);

     printf("Duplicates removed: %u\n", dr_count);

     return 0;
}
