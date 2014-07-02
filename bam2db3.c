#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <limits.h>

#include <sqlite3.h>

#include <bam/sam.h>
#include <bam/bam.h>

#include "errors.h"

static const size_t win_len = 1000;

char sql_select_chrid[] = "SELECT id FROM chromosomes WHERE name = ?;";
char sql_select_animalid[] = "SELECT id FROM animals WHERE name = ?;";
char sql_nuc_insert[] = "INSERT INTO nucleotides VALUES("
     "?1, ?2, ?3, "
     "?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)";
char sql_codon_insert[] = "INSERT INTO codons VALUES("
     " ?1, ?2, ?3, ?4,"
     " ?5, ?6, ?7, ?8, ?9,?10,?11,?12,?13,?14,?15,?16,?17,?18,?19,?20,"
     "?21,?22,?23,?24,?25,?26,?27,?28,?29,?30,?31,?32,?33,?34,?35,?36,"
     "?37,?38,?39,?40,?41,?42,?43,?44,?45,?46,?47,?48,?49,?50,?51,?52,"
     "?53,?54,?55,?56,?57,?58,?59,?60,?61,?62,?63,?64,?65,?66,?67,?68,"
     "?69);";

void memerror(void *ptr)
{
     if (!ptr) {
          printf("Couldn't allocate memory.\n");
          exit(MEM_ERROR);
     }
}

/* nucleotide counts for a reference position */
typedef struct pos_nucs {
     sqlite3_int64 forward[4];
     sqlite3_int64 reverse[4];
     sqlite3_int64 dels;
} PosNucs;

/* codon counts for a cds position */
typedef struct pos_cods {
     sqlite3_int64 codons[4][4][4];
     sqlite3_int64 dels;
} PosCods;

static size_t buffer_init_size = 64;
typedef struct {
     size_t beg, size;
     PosNucs *nucs;
     PosCods *cods;
} Buffer;

int init_buffer(Buffer *buf)
{
     buf->beg = 0;
     buf->size = buffer_init_size;
     buf->nucs = calloc(buffer_init_size, sizeof(PosNucs));
     if (!buf->nucs) {
          return -1;
     }
     buf->cods = calloc(buffer_init_size, sizeof(PosCods));
     if (!buf->cods) {
          return -1;
     }
     return 0;
}

void free_buffer(Buffer *buf)
{
     buf->beg = 0;
     buf->size = 0;
     free(buf->nucs);
     free(buf->cods);
}

/* makes buffer at least as big as new_size */
int resize_buffer(Buffer *buf, size_t new_size)
{
     size_t old_size = buf->size, diff;
     while (buf->size < new_size) {
          buf->size <<= 1;
     }
     buf->nucs = realloc(buf->nucs, buf->size * sizeof(PosNucs));
     if (!buf->nucs) {
          return -1;
     }
     buf->cods = realloc(buf->cods, buf->size * sizeof(PosNucs));
     if (!buf->cods) {
          return -1;
     }
     diff = buf->size - old_size;
     if (diff > 0) {
          memset(buf->nucs + old_size, 0, diff);
          /* can use memcpy because max_size is doubled */
          memcpy(buf->nucs + buf->beg + diff,
                 buf->nucs + buf->beg, old_size - buf->beg);

          memset(buf->cods + old_size, 0, diff);
          memcpy(buf->cods + buf->beg + diff,
                 buf->cods + buf->beg, old_size - buf->beg);
     }

     return 0;
}

void buffer_insert_read(Buffer *buf, bam1_t *read)
{
     
}

/* print sqlite error */
static int nerrors = 0, maxerrors = 10;
void sql_error(char **errormessage)
{
     if (*errormessage) {
          fprintf(stderr, "SQLite3 error: %s\n", *errormessage);
          sqlite3_free(*errormessage);
          *errormessage = NULL;
          nerrors++;
          if (nerrors >= maxerrors) {
               fprintf(stderr, "(Stopping after %d errors.)\n", maxerrors);
               exit(SQL_ERROR);
          }
     }
}

void sqlite3_step_onerow(sqlite3_stmt *stmt)
{
     int res_code;
     res_code = sqlite3_step(stmt);
     if (res_code != SQLITE_ROW) {
          fprintf(stderr, "SQLite3 error: %d (%s)\n",
                  res_code, sqlite3_sql(stmt));
          exit(SQL_ERROR);
     }
}

static int beg_to_db(sqlite3 *db, sqlite3_stmt *stmt, int32_t pos, Buffer *buf)
{
     int res_code;
     size_t beg = buf->beg;

     sqlite3_bind_int64(stmt, 4, pos + 1);
     sqlite3_bind_int64(stmt, 5, buf->nucs[beg].forward[0]);
     sqlite3_bind_int64(stmt, 6, buf->nucs[beg].forward[1]);
     sqlite3_bind_int64(stmt, 7, buf->nucs[beg].forward[2]);
     sqlite3_bind_int64(stmt, 8, buf->nucs[beg].forward[3]);
     sqlite3_bind_int64(stmt, 9, buf->nucs[beg].reverse[0]);
     sqlite3_bind_int64(stmt, 10, buf->nucs[beg].reverse[1]);
     sqlite3_bind_int64(stmt, 11, buf->nucs[beg].reverse[2]);
     sqlite3_bind_int64(stmt, 12, buf->nucs[beg].reverse[3]);
     sqlite3_bind_int64(stmt, 13, buf->nucs[beg].dels);

     res_code = sqlite3_step(stmt);
     if (res_code != SQLITE_DONE) {
          fprintf(stderr, "SQLite3 error: %d!\n", res_code);
          exit(SQL_ERROR);
     }

     sqlite3_reset(stmt);
     return 0;
}

int main(int argc, char *argv[])
{
     int i;

     char *dbname;
     char *bamfilename;
     char *animal;
     int day;

     sqlite3 *db = NULL;
     sqlite3_stmt *stmt;
     int dbopencode = 0;
     char *errormessage = NULL;

     samfile_t *samin;

     int32_t n_chr;
     sqlite3_int64 *db_chrids = NULL;
     sqlite3_int64 animal_id;

     bam1_t read;
     int32_t cur_tid, cur_pos, read_len;
     Buffer buf;
     memset(&read, 0, sizeof(bam1_t));

     if (argc < 5) {
          printf("Usage: %s database bamfile animal day\n", argv[0]);
          exit(ARG_ERROR);
     }
     else {
          fprintf(stderr, "DB name: %s... ", argv[1]);
          dbname = argv[1];
          bamfilename = argv[2];
          animal = argv[3];
          day = strtoul(argv[4], NULL, 10);
     }
     /* try to open db */
     dbopencode = sqlite3_open(dbname, &db);
     if (SQLITE_OK == dbopencode) {
          fprintf(stderr,"OK\n");
     }
     else {
          printf("SQLite3 error: %d\n", dbopencode);
          exit(SQL_ERROR);
     }

     /* try to open sam/bam file */
     samin = samopen(bamfilename, "rb", NULL);
     if (!samin) {
          fprintf(stderr, "Error opening file %s\n", bamfilename);
          exit(BAM_ERROR);
     }

     /* get chromosome info */
     n_chr = samin->header->n_targets;
     db_chrids = malloc(n_chr * sizeof(sqlite3_int64));
     memerror(db_chrids);
     sqlite3_prepare_v2(db, sql_select_chrid, sizeof(sql_select_chrid),
                       &stmt, NULL);
     for (i = 0; i < n_chr; i++) {
          sqlite3_bind_text(stmt, 1, samin->header->target_name[i], -1,
                            SQLITE_STATIC);
          sqlite3_step_onerow(stmt);
          db_chrids[i] = sqlite3_column_int64(stmt, 0);
          sqlite3_reset(stmt);
     }
     sqlite3_finalize(stmt);

     /* get id of animal */
     sqlite3_prepare_v2(db, sql_select_animalid, sizeof(sql_select_animalid),
                        &stmt, NULL);
     sqlite3_bind_text(stmt, 1, animal, -1, SQLITE_STATIC);
     sqlite3_step_onerow(stmt);
     animal_id = sqlite3_column_int64(stmt, 0);
     sqlite3_finalize(stmt);

     /* increases speed of insertion (means if the program crashes the db is */
     /* left invalid) */
     sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errormessage);

     /* start prepared statement for nucleotide insert */
     sqlite3_prepare_v2(db, sql_nuc_insert, sizeof(sql_nuc_insert),
                        &stmt, NULL);
     sqlite3_bind_int64(stmt, 1, animal_id);
     sqlite3_bind_int(stmt, 2, day);

     sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &errormessage);
     sql_error(&errormessage);

     cur_tid = cur_pos = 0;
     init_buffer(&buf);
     while(samread(samin, &read) > 0) {
          if (cur_tid < read.core.tid) {
               /* flush buffer */
               while (cur_pos < samin->header->target_len[cur_tid]) {
                    beg_to_db(db, stmt, cur_pos, &buf);
                    cur_pos++;
                    buf.beg++;
               }
               cur_tid = read.core.tid;
          }
          while (cur_pos < read.core.pos) {
               beg_to_db(db, stmt, cur_pos, &buf);
               cur_pos++;
               buf.beg++;
          }
          read_len = bam_cigar2qlen(&read.core, bam1_cigar(&read));
          if (read_len > buf.size) {
               resize_buffer(&buf, read_len);
          }
          buffer_insert_read(&buf, &read);
     }

     sqlite3_exec(db, "COMMIT TRANSACTION;", NULL, NULL, &errormessage);
     sql_error(&errormessage);
     sqlite3_finalize(stmt);

     free(db_chrids);
     samclose(samin);
     sqlite3_close(db);

     return 0;
}
