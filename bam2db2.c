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

char sql_count_cds[] = "SELECT count(*) FROM cds WHERE chromosome = ?;";
char sql_select_cds[] = "SELECT id FROM cds WHERE chromosome = ?;";
char sql_count_cds_regions[] = "SELECT count(*) FROM cds_regions where cds = ?;";
char sql_select_cds_regions[] = "SELECT start, end FROM cds_regions where cds = ?;";
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

/* region in a chromosome that is a cds (beg and end incluive) */
typedef struct cds_region {
     sqlite3_int64 ref_beg, ref_end;
} CDSRegion;

typedef struct cds {
     sqlite_int64 id;
     sqlite_int64 nregions;
     CDSRegion *regions;
} CDS;

typedef struct chromosome {
     sqlite3_int64 id;
     sqlite3_int64 ncds;
     CDS *cds;
} Chromosome;

/* a bit of a cds that overlaps with the window
 *
 * ref_beg,ref_end are the positions in the reference of the leftmost
 * nucleotide in the codon (first for forward frames, last for reverse)
 *
 * cds_beg,cds_end are the positions of the window in cds coordinates */
typedef struct cds_window {
     CDS *cds;
     sqlite3_int64 ref_beg, reg_end;
     sqlite3_int64 cds_beg, cds_end;
     PosCods *codarrays;
} CDSWin;

/* window of reference we are looking at */
typedef struct win_data {
     int32_t beg;
     PosNucs *nucarrays;
     CDSWin *regions;
} WinData;

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
          fprintf(stderr, "SQLite3 error: %d (%s)\n", res_code, sqlite3_sql(stmt));
          exit(SQL_ERROR);
     }
}

/* gets chromosome information from database */
Chromosome *get_chromosomes(sqlite3 *db, int32_t n_chr, char **chr_names)
{
     int32_t i, j, k;
     Chromosome *chromosomes;
     sqlite3_stmt *id_stmt, *ncds_stmt, *cds_stmt, *ncdsr_stmt, *cdsr_stmt;

     chromosomes = malloc(sizeof(Chromosome) * n_chr);

     sqlite3_prepare_v2(db, sql_select_chrid, sizeof(sql_select_chrid),
                        &id_stmt, NULL);
     sqlite3_prepare_v2(db, sql_count_cds, sizeof(sql_count_cds),
                        &ncds_stmt, NULL);
     sqlite3_prepare_v2(db, sql_select_cds, sizeof(sql_select_cds),
                        &cds_stmt, NULL);
     sqlite3_prepare_v2(db, sql_count_cds_regions, sizeof(sql_count_cds_regions),
                        &ncdsr_stmt, NULL);
     sqlite3_prepare_v2(db, sql_select_cds_regions, sizeof(sql_select_cds_regions),
                        &cdsr_stmt, NULL);

     for (i = 0; i < n_chr; ++i) {
          /* id */
          sqlite3_bind_text(id_stmt, 1, chr_names[i], -1,
                            SQLITE_STATIC);
          sqlite3_step_onerow(id_stmt);
          chromosomes[i].id = sqlite3_column_int64(id_stmt, 0);
          sqlite3_reset(id_stmt);

          /* cds */
          sqlite3_bind_int64(ncds_stmt, 1, chromosomes[i].id);
          sqlite3_step_onerow(ncds_stmt);
          chromosomes[i].ncds = sqlite3_column_int64(ncds_stmt, 0);
          sqlite3_reset(ncds_stmt);

          chromosomes[i].cds =
               malloc(sizeof(CDS) * chromosomes[i].ncds);
          sqlite3_bind_int64(cds_stmt, 1, chromosomes[i].id);
          for (j = 0; j < chromosomes[i].ncds; ++j) {
               sqlite3_step_onerow(cds_stmt);
               chromosomes[i].cds[j].id = sqlite3_column_int64(cds_stmt, 0);
               /* cds regions */
               sqlite3_bind_int64(ncdsr_stmt, 1, chromosomes[i].cds[j].id);
               sqlite3_step_onerow(ncdsr_stmt);
               chromosomes[i].cds[j].nregions = sqlite3_column_int64(ncdsr_stmt, 0);
               sqlite3_reset(ncdsr_stmt);

               chromosomes[i].cds[j].regions =
                    malloc(sizeof(CDSRegion) * chromosomes[i].cds[j].nregions);
               sqlite3_bind_int64(cdsr_stmt, 1, chromosomes[i].cds[j].id);
               for (k = 0; k < chromosomes[i].cds[j].nregions; ++k) {
                    sqlite3_step_onerow(cdsr_stmt);
                    chromosomes[i].cds[j].regions[k].ref_beg =
                         sqlite3_column_int64(cdsr_stmt, 0);
                    chromosomes[i].cds[j].regions[k].ref_end =
                         sqlite3_column_int64(cdsr_stmt, 1);
               }
               sqlite3_reset(cdsr_stmt);
          }
          sqlite3_reset(cds_stmt);
     }

     sqlite3_finalize(cdsr_stmt);
     sqlite3_finalize(ncdsr_stmt);
     sqlite3_finalize(cds_stmt);
     sqlite3_finalize(ncds_stmt);
     sqlite3_finalize(id_stmt);
     return chromosomes;
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
                    win->nucarrays[win_pos].dels++;
               }
               else {
                    if (pl[i].b->core.flag & BAM_FREVERSE) {
                         nucs = win->nucarrays[win_pos].reverse;
                    }
                    else {
                         nucs = win->nucarrays[win_pos].forward;
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
                        uint32_t chr_len, WinData *win)
{
     int i;
     int res_code;
     size_t win_len_used = chr_len - win->beg;
     win_len_used = win_len_used < win_len ? win_len_used : win_len;

     for (i = 0; i < win_len_used; ++i) {
          sqlite3_bind_int64(stmt, 4, win->beg + i + 1);
          sqlite3_bind_int64(stmt, 5, win->nucarrays[i].forward[0]);
          sqlite3_bind_int64(stmt, 6, win->nucarrays[i].forward[1]);
          sqlite3_bind_int64(stmt, 7, win->nucarrays[i].forward[2]);
          sqlite3_bind_int64(stmt, 8, win->nucarrays[i].forward[3]);
          sqlite3_bind_int64(stmt, 9, win->nucarrays[i].reverse[0]);
          sqlite3_bind_int64(stmt, 10, win->nucarrays[i].reverse[1]);
          sqlite3_bind_int64(stmt, 11, win->nucarrays[i].reverse[2]);
          sqlite3_bind_int64(stmt, 12, win->nucarrays[i].reverse[3]);
          sqlite3_bind_int64(stmt, 13, win->nucarrays[i].dels);

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
     int32_t n_chr;
     Chromosome *chromosomes = NULL;
     sqlite3_int64 animal_id;

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

     /* get chromosome info */
     n_chr = bamin->header->n_targets;
     chromosomes = get_chromosomes(db, n_chr, bamin->header->target_name);

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

     win.nucarrays = calloc(win_len, sizeof(PosNucs));

     buf = bam_plbuf_init(&pileup_func, &win);
     /* disable maximum pileup depth */
     bam_plp_set_maxcnt(buf->iter, INT_MAX);

     /* start prepared statement for nucleotide insert */
     sqlite3_prepare_v2(db, sql_nuc_insert, sizeof(sql_nuc_insert),
                        &stmt, NULL);
     sqlite3_bind_int64(stmt, 1, animal_id);
     sqlite3_bind_int(stmt, 2, day);
     
     sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &errormessage);
     sql_error(&errormessage);
     for (i = 0; i < n_chr; ++i) {
          for (j = 0; j < chromosomes[i].ncds; ++j) {
               printf("   id: %d, beg: %d, end: %d\n",
                      chromosomes[i].cds_regions[j].id,
                      chromosomes[i].cds_regions[j].ref_beg,
                      chromosomes[i].cds_regions[j].ref_end);
          }
          fprintf(stderr, "%s", bamin->header->target_name[i]);
          sqlite3_bind_int64(stmt, 3, chromosomes[i].id);
          for (j = 0; j < bamin->header->target_len[i]; j += win_len) {
               printf("win: %d -- %d\n", j, j + win_len);
               for (int k = 0; k < chromosomes[i].ncds; ++k) {
                    if ((chromosomes[i].cds_regions[k].ref_beg >= j
                         && chromosomes[i].cds_regions[k].ref_beg < (j + win_len))
                        || (chromosomes[i].cds_regions[k].ref_end >= j
                            && chromosomes[i].cds_regions[k].ref_beg < (j + win_len))) {
                         printf("cds id: %d\n", chromosomes[i].cds_regions[k].id);
                    }
               }
               win.beg = j;
               bam_fetch(bamin->x.bam, bamidx,
                         i, j, j + win_len,
                         buf, &fetch_func);
               bam_plbuf_push(0, buf);
               window_to_db(db, stmt,
                            bamin->header->target_len[i], &win);
               bam_plbuf_reset(buf);
               memset(win.nucarrays, 0, sizeof(PosNucs)*win_len);
               fprintf(stderr, ".");
          }
          fprintf(stderr, "\n");
     }
     sqlite3_exec(db, "COMMIT TRANSACTION;", NULL, NULL, &errormessage);
     sql_error(&errormessage);
     sqlite3_finalize(stmt);

     free(win.nucarrays);
     bam_index_destroy(bamidx);
     samclose(bamin);
     free(chromosomes);
     sqlite3_close(db);

     return 0;
}
