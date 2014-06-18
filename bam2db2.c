#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <limits.h>

#include <sqlite3.h>

#include <bam/sam.h>
#include <bam/bam.h>

#include "errors.h"

static const size_t win_len = 500;

typedef struct pos_nucs {
     sqlite3_int64 forward[4];
     sqlite3_int64 reverse[4];
     sqlite3_int64 dels;
} PosNucs;

char sql_select_chrid[] = "SELECT id FROM chromosomes WHERE name = ?;";
char sql_select_animalid[] = "SELECT id FROM animals WHERE name = ?;";
char sql_nuc_insert[] = "INSERT INTO nucleotides VALUES(?1, ?2, ?3, "
     "?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)";

typedef struct win_data {
     int32_t beg;
     PosNucs *posarr;
} WinData;

/* print sqlite error */
static int nerrors = 0, maxerrors = 10;
void do_error(char **errormessage)
{
     if (*errormessage) {
          fprintf(stderr, "SQLite3 error: %s\n", *errormessage);
          sqlite3_free(*errormessage);
          *errormessage = NULL;
          nerrors++;
          if (nerrors >= maxerrors) {
               fprintf(stderr, "(Stopping after %d errors.)\n", maxerrors);
               exit(SQLITE_ERROR);
          }
     }
}

/* callback for bam_fetch() */
static int fetch_func(const bam1_t *b, void *data)
{
     bam_plbuf_t *buf = data;
     if (b->core.flag & BAM_FPROPER_PAIR) {
          bam_plbuf_push(b, buf);
     }
     return 0;
}

/* callback for bam_plbuf_init() */
static int pileup_func(uint32_t tid, uint32_t pos, int n,
                       const bam_pileup1_t *pl, void *data)
{
     WinData *win = data;
     int win_pos, i;
     sqlite_int64 *nucs;

     win_pos = pos - win->beg;
     if (win_pos >= 0 && win_pos < win_len) {
          for (i = 0; i < n; ++i) {
               if (pl[i].is_del) {
                    win->posarr[win_pos].dels++;
               }
               else {
                    if (pl[i].b->core.flag & BAM_FREVERSE) {
                         nucs = win->posarr[win_pos].reverse;
                    }
                    else {
                         nucs = win->posarr[win_pos].forward;
                    }
                    switch(bam1_seqi(bam1_seq(pl[i].b), pl[i].qpos)) {
                    case 1:
                         nucs[0]++;
                         break;
                    case 2:
                         nucs[1]++;
                         break;
                    case 4:
                         nucs[2]++;
                         break;
                    case 8:
                         nucs[3]++;
                         break;
                    }
               }
          }
     }

     return 0;
}

static int window_to_db(sqlite3 *db, sqlite3_stmt *stmt,
                        sqlite3_int64 chr_id, uint32_t chr_len, WinData *win) 
{
     int i;
     int res_code;
     size_t win_len_used = chr_len - win->beg;
     win_len_used = win_len_used < win_len ? win_len_used : win_len;

     for (i = 0; i < win_len_used; ++i) {
          sqlite3_bind_int64(stmt, 3, chr_id);
          sqlite3_bind_int64(stmt, 4, win->beg + i + 1);
          sqlite3_bind_int64(stmt, 5, win->posarr[i].forward[0]);
          sqlite3_bind_int64(stmt, 6, win->posarr[i].forward[1]);
          sqlite3_bind_int64(stmt, 7, win->posarr[i].forward[2]);
          sqlite3_bind_int64(stmt, 8, win->posarr[i].forward[3]);
          sqlite3_bind_int64(stmt, 9, win->posarr[i].reverse[0]);
          sqlite3_bind_int64(stmt, 10, win->posarr[i].reverse[1]);
          sqlite3_bind_int64(stmt, 11, win->posarr[i].reverse[2]);
          sqlite3_bind_int64(stmt, 12, win->posarr[i].reverse[3]);
          sqlite3_bind_int64(stmt, 13, win->posarr[i].dels);

          res_code = sqlite3_step(stmt);
          if (res_code != SQLITE_DONE) {
               fprintf(stderr, "SQLite3 error: %d!\n", res_code);
               exit(SQL_ERROR);
          }

          sqlite3_reset(stmt);
     }
     return 0;
}
     
int main(int argc, char *argv[])
{
     /* char *progname; */
     /* int verbose; */

     char *dbname;
     char *bamfilename;
     char *animal;
     int day;

     sqlite3 *db = NULL;
     sqlite3_stmt *stmt;
     samfile_t *bamin;
     bam_index_t *bamidx;

     int dbopencode = 0;
     char *errormessage = NULL;
     int n_chr;
     sqlite3_int64 *chr_ids = NULL;
     sqlite3_int64 animal_id;
     int res_code;

     bam_plbuf_t *buf;
     size_t i, j;
     WinData win;

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

     /* try to open bam file */
     bamin = samopen(bamfilename, "rb", NULL);
     if (!bamin) {
          fprintf(stderr, "Error opening bamfile %s\n", bamfilename);
          exit(BAM_ERROR);
     }
     /* try to open index */
     bamidx = bam_index_load(bamfilename);
     if (!bamidx) {
          fprintf(stderr, "Error opening index for %s\n", bamfilename);
          exit(BAM_ERROR);
     }

     /* get ids of chromosomes */
     n_chr = bamin->header->n_targets;
     chr_ids = malloc(sizeof(sqlite3_int64) * n_chr);
     sqlite3_prepare_v2(db, sql_select_chrid, sizeof(sql_select_chrid),
                        &stmt, NULL);
     for (i = 0; i < n_chr; ++i) {
          sqlite3_bind_text(stmt, 1, bamin->header->target_name[i], -1,
                            SQLITE_STATIC);
          res_code = sqlite3_step(stmt);
          if (res_code != SQLITE_ROW) {
               fprintf(stderr, "SQLite3 error: %d!\n", res_code);
               exit(SQL_ERROR);
          }
          chr_ids[i] = sqlite3_column_int64(stmt, 0);
          sqlite3_reset(stmt);
     }
     sqlite3_finalize(stmt);
     /* get id of animal */
     sqlite3_prepare_v2(db, sql_select_animalid, sizeof(sql_select_animalid),
                        &stmt, NULL);
     sqlite3_bind_text(stmt, 1, animal, -1, SQLITE_STATIC);
     res_code = sqlite3_step(stmt);
     if (res_code != SQLITE_ROW) {
          fprintf(stderr, "SQLite3 error: %d!\n", res_code);
          exit(SQL_ERROR);
     }
     animal_id = sqlite3_column_int64(stmt, 0);
     sqlite3_finalize(stmt);
     
     /* increases speed of insertion (means if the program crashes the db is */
     /* left invalid) */
     sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errormessage);

     win.posarr = calloc(win_len, sizeof(PosNucs));

     buf = bam_plbuf_init(&pileup_func, &win);
     /* disable maximum pileup depth */
     bam_plp_set_maxcnt(buf->iter, INT_MAX);

     /* start prepared statement for nucleotide insert */
     sqlite3_prepare_v2(db, sql_nuc_insert, sizeof(sql_nuc_insert),
                        &stmt, NULL);
     sqlite3_bind_int64(stmt, 1, animal_id);
     sqlite3_bind_int(stmt, 2, day);
     
     sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &errormessage);
     do_error(&errormessage);
     for (i = 0; i < n_chr; ++i) {
          fprintf(stderr, "%s", bamin->header->target_name[i]);
          for (j = 0; j < bamin->header->target_len[i]; j += win_len) {
               win.beg = j;
               bam_fetch(bamin->x.bam, bamidx,
                         i, j, j + win_len,
                         buf, &fetch_func);
               bam_plbuf_push(0, buf);
               window_to_db(db, stmt,
                            chr_ids[i], bamin->header->target_len[i], &win);
               bam_plbuf_reset(buf);
               memset(win.posarr, 0, sizeof(PosNucs)*win_len);
               fprintf(stderr, ".");
          }
          fprintf(stderr, "\n");
     }
     sqlite3_exec(db, "COMMIT TRANSACTION;", NULL, NULL, &errormessage);
     do_error(&errormessage);
     sqlite3_finalize(stmt);

     free(win.posarr);
     bam_index_destroy(bamidx);
     samclose(bamin);
     free(chr_ids);
     sqlite3_close(db);

     return 0;
}
