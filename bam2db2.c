#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include <sqlite3.h>

#include <bam/sam.h>
#include <bam/bam.h>

static sqlite3 *db = NULL;
static sqlite3_stmt *stmt;
static samfile_t *bamin;
static bam1_t *bam_read;
static char *dbname, *bamfilename, *animal;
static int day;

char sql_select_chrid[] = "select rowid from chromosomes where name = ?;";

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
     /* uint8_t *dr = NULL; */
     /* uint32_t dr_count = 0; */
     int i;
     int n_chr;
     sqlite3_int64 *chr_rowids = NULL;
     int res_code;

     int32_t tid;

     unsigned *chr_mapped = NULL;
     unsigned *chr_unmapped = NULL;

     bam_read = bam_init1();

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

     /* try to open bam file */
     bamin = samopen(bamfilename, "rb", NULL);
     if (!bamin) {
          fprintf(stderr, "Error opening bamfile %s\n", bamfilename);
          exit(1);
     }

     /* get rowids of chromosomes */
     n_chr = bamin->header->n_targets;
     chr_rowids = malloc(sizeof(sqlite3_int64) * n_chr);
     printf("targets: %d\n", n_chr);
     sqlite3_prepare_v2(db, sql_select_chrid, strlen(sql_select_chrid)+1,
                        &stmt, NULL);
     for (i = 0; i < n_chr; ++i) {
          sqlite3_bind_text(stmt, 1, bamin->header->target_name[i], -1,
                            SQLITE_STATIC);
          res_code = sqlite3_step(stmt);
          if (res_code != SQLITE_ROW) {
               fprintf(stderr, "SQLite3 error: %d!\n", res_code);
               exit(1);
          }
          chr_rowids[i] = sqlite3_column_int64(stmt, 0);
          sqlite3_reset(stmt);
          printf("%s, %d, %ld\n", bamin->header->target_name[i],
                 i, (long)chr_rowids[i]);
     }
     sqlite3_finalize(stmt);
     
     /* increases speed of insertion (means if the program crashes the db is
      * left invalid) */
     sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errormessage);

     /* while (samread(bamin, bam_read) >= 0) { */
     /*      dr = bam_aux_get(bam_read, "DR"); */
     /*      if (dr) { */
     /*           dr_count += bam_aux2i(dr); */
     /*      } */
     /* } */

     chr_mapped = calloc(n_chr + 1, sizeof(unsigned));
     chr_unmapped = calloc(n_chr + 1, sizeof(unsigned));
     while (samread(bamin, bam_read) >= 0) {
          tid = bam_read->core.tid;
          if (bam_read->core.flag & BAM_FUNMAP) {
               /* unmapped */
               chr_unmapped[tid < 0 ? n_chr : tid]++;
          }
          else {
               chr_mapped[tid < 0 ? n_chr : tid]++;
          }
     }
     for (i = 0; i < n_chr; ++i) {
          printf("%-4s %5"PRIu32" %8u %8u\n",
                 bamin->header->target_name[i], bamin->header->target_len[i],
                 chr_mapped[i], chr_unmapped[i]);
     }
     printf("%-4s %5d %8u %8u\n", "*", 0, 
            chr_mapped[n_chr], chr_unmapped[n_chr]);

     samclose(bamin);
     bam_destroy1(bam_read);
     free(chr_rowids);
     free(chr_mapped);
     free(chr_unmapped);

     sqlite3_close(db);

     /* printf("Duplicates removed: %u\n", dr_count); */

     return 0;
}
