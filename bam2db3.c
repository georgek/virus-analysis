#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <limits.h>

/* #define NDEBUG */
#include <assert.h>

#include <sqlite3.h>

#include <bam/sam.h>
#include <bam/bam.h>

#include "utils.h"
#include "errors.h"

char sql_select_chrid[] = "SELECT id FROM chromosomes WHERE name = ?;";
char sql_select_animalid[] = "SELECT id FROM animals WHERE name = ?;";
char sql_nuc_insert[] = "INSERT INTO "
     "nucleotides(animal, day, chromosome, position,"
     "Af, Cf, Gf, Tf, Ar, Cr, Gr, Tr, D) "
     "VALUES("
     "?1, ?2, ?3, ?4,"
     "?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)";
char sql_cod_insert[] = "INSERT INTO "
     "codons(animal, day, chromosome, position,"
     "AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT,"
     "AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT,"
     "CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT,"
     "CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT,"
     "GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT,"
     "GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT,"
     "TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT,"
     "TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT) "
     "VALUES("
     " ?1, ?2, ?3, ?4,"
     " ?5, ?6, ?7, ?8, ?9,?10,?11,?12,?13,?14,?15,?16,?17,?18,?19,?20,"
     "?21,?22,?23,?24,?25,?26,?27,?28,?29,?30,?31,?32,?33,?34,?35,?36,"
     "?37,?38,?39,?40,?41,?42,?43,?44,?45,?46,?47,?48,?49,?50,?51,?52,"
     "?53,?54,?55,?56,?57,?58,?59,?60,?61,?62,?63,?64,?65,?66,?67,?68);";

void memerror(void *ptr)
{
     if (!ptr) {
          printf("Couldn't allocate memory.\n");
          exit(MEM_ERROR);
     }
}

typedef sqlite3_int64 Nucs[4];
typedef sqlite3_int64 Cods[4*4*4];
typedef sqlite3_int64 Dels;

static size_t buffer_init_size = 64;
struct buffer {
     size_t beg, size;
     /* forward and reverse nuc arrays */
     Nucs *nucs[2];
     Cods *cods;
     Dels *dels;
};

int init_buffer(struct buffer *buf)
{
     buf->beg = 0;
     buf->size = buffer_init_size;
     buf->nucs[0] = calloc(buffer_init_size, sizeof(Nucs));
     if (!buf->nucs[0]) {
          return -1;
     }
     buf->nucs[1] = calloc(buffer_init_size, sizeof(Nucs));
     if (!buf->nucs[1]) {
          return -1;
     }
     buf->cods = calloc(buffer_init_size, sizeof(Cods));
     if (!buf->cods) {
          return -1;
     }
     buf->dels = calloc(buffer_init_size, sizeof(Dels));
     if (!buf->dels) {
          return -1;
     }
     return 0;
}

void free_buffer(struct buffer *buf)
{
     buf->beg = 0;
     buf->size = 0;
     free(buf->nucs[0]);
     free(buf->nucs[1]);
     free(buf->cods);
     free(buf->dels);
}

/* makes buffer at least as big as new_size */
int resize_buffer(struct buffer *buf, size_t new_size)
{
     size_t old_size = buf->size, diff;
     while (buf->size < new_size) {
          buf->size <<= 1;
     }
     buf->nucs[0] = realloc(buf->nucs[0], buf->size * sizeof(Nucs));
     if (!buf->nucs[0]) {
          return -1;
     }
     buf->nucs[1] = realloc(buf->nucs[1], buf->size * sizeof(Nucs));
     if (!buf->nucs[1]) {
          return -1;
     }
     buf->cods = realloc(buf->cods, buf->size * sizeof(Cods));
     if (!buf->cods) {
          return -1;
     }
     buf->dels = realloc(buf->dels, buf->size * sizeof(Dels));
     if (!buf->dels) {
          return -1;
     }
     diff = buf->size - old_size;
     if (diff > 0) {
          assert(diff >= old_size);
          /* can use memcpy because size is doubled */
          memcpy(buf->nucs[0] + buf->beg + diff,
                 buf->nucs[0] + buf->beg,
                 (old_size - buf->beg)*sizeof(Nucs));
          memcpy(buf->nucs[1] + buf->beg + diff,
                 buf->nucs[1] + buf->beg,
                 (old_size - buf->beg)*sizeof(Nucs));
          memcpy(buf->cods + buf->beg + diff,
                 buf->cods + buf->beg,
                 (old_size - buf->beg)*sizeof(Cods));
          memcpy(buf->dels + buf->beg + diff,
                 buf->dels + buf->beg,
                 (old_size - buf->beg)*sizeof(Dels));
          /* zero the new bit */
          memset(buf->nucs[0] + buf->beg,
                 0,
                 diff*sizeof(Nucs));
          memset(buf->nucs[1] + buf->beg,
                 0,
                 diff*sizeof(Nucs));
          memset(buf->cods + buf->beg,
                 0,
                 diff*sizeof(Cods));
          memset(buf->dels + buf->beg,
                 0,
                 diff*sizeof(Dels));
          buf->beg += diff;
     }

     return 0;
}

int buffer_index(struct buffer *buf, int index)
{
     return (buf->beg + index)%buf->size;
}

void buffer_next(struct buffer *buf) 
{
     memset(buf->nucs[0] + buf->beg, 0, sizeof(Nucs));
     memset(buf->nucs[1] + buf->beg, 0, sizeof(Nucs));
     memset(buf->cods + buf->beg, 0, sizeof(Cods));
     memset(buf->dels + buf->beg, 0, sizeof(Dels));
     buf->beg = (buf->beg + 1)%buf->size;
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

static int bambasetable[] = {-1, 0, 1, -1, 2, -1, -1, -1, 3,
                             -1, -1, -1, -1, -1, -1, -1};
void buffer_insert_read(struct buffer *buf, bam1_t *read)
{
     int buf_pos, seq_pos, cig_pos, cig_len;
     int bambase, bambase1, bambase2;
     Nucs *nucs;

     nucs = *(buf->nucs + ((read->core.flag & BAM_FREVERSE) >> 4));
     for(buf_pos = 0, seq_pos = 0, cig_pos = 0;
         cig_pos < read->core.n_cigar;
         cig_pos++) {
          cig_len = bam1_cigar(read)[cig_pos] >> 4;
          switch(bam1_cigar(read)[cig_pos] & 0xF) {
          case BAM_CMATCH:
               for (; cig_len > 0; cig_len--, buf_pos++, seq_pos++) {
                    bambase = bam1_seqi(bam1_seq(read), seq_pos);
                    bambase = bambasetable[bambase];
                    if (bambase > 0) {
                         nucs[buffer_index(buf, buf_pos)][bambasetable[bambase]]++;
                         if (cig_len >= 3) {
                              bambase1 = bam1_seqi(bam1_seq(read), seq_pos+1);
                              bambase1 = bambasetable[bambase1];
                              bambase2 = bam1_seqi(bam1_seq(read), seq_pos+2);
                              bambase2 = bambasetable[bambase2];
                              buf->cods[buffer_index(buf, buf_pos)][bambase*16+
                                                                    bambase1*4+
                                                                    bambase2]++;
                         }
                    }
               }
               break;
          case BAM_CINS:
               seq_pos += cig_len;
               break;
          case BAM_CREF_SKIP:
          case BAM_CDEL:
               for (; cig_len > 0; cig_len--, buf_pos++) {
                    buf->dels[buffer_index(buf, buf_pos)]++;
               }
               break;
          case BAM_CSOFT_CLIP:
               seq_pos += cig_len;
               break;
          case BAM_CHARD_CLIP:
          case BAM_CPAD:
               break;
          default:
               fprintf(stderr, "Unrecognised CIGAR operation.\n");
               break;
          }
     }
}

int beg_to_db(sqlite3 *db, sqlite3_stmt *nuc_stmt, sqlite3_stmt *cod_stmt,
              int32_t pos, struct buffer *buf)
{
     int i, res_code;
     size_t beg = buf->beg;

     /* do nucs */
     sqlite3_bind_int64(nuc_stmt, 4, pos + 1);
     sqlite3_bind_int64(nuc_stmt, 5, buf->nucs[0][beg][0]);
     sqlite3_bind_int64(nuc_stmt, 6, buf->nucs[0][beg][1]);
     sqlite3_bind_int64(nuc_stmt, 7, buf->nucs[0][beg][2]);
     sqlite3_bind_int64(nuc_stmt, 8, buf->nucs[0][beg][3]);
     sqlite3_bind_int64(nuc_stmt, 9, buf->nucs[1][beg][0]);
     sqlite3_bind_int64(nuc_stmt, 10, buf->nucs[1][beg][1]);
     sqlite3_bind_int64(nuc_stmt, 11, buf->nucs[1][beg][2]);
     sqlite3_bind_int64(nuc_stmt, 12, buf->nucs[1][beg][3]);
     sqlite3_bind_int64(nuc_stmt, 13, buf->dels[beg]);

     res_code = sqlite3_step(nuc_stmt);
     if (res_code != SQLITE_DONE) {
          fprintf(stderr, "SQLite3 error: %d!\n", res_code);
          exit(SQL_ERROR);
     }

     /* do cods */
     sqlite3_bind_int64(cod_stmt, 4, pos + 1);
     for (i = 0; i < 64; i++) {
          sqlite3_bind_int64(cod_stmt, i + 5, buf->cods[beg][i]);
     }
     res_code = sqlite3_step(cod_stmt);
     if (res_code != SQLITE_DONE) {
          fprintf(stderr, "SQLite3 error: %d!\n", res_code);
          exit(SQL_ERROR);
     }

     sqlite3_reset(nuc_stmt);
     sqlite3_reset(cod_stmt);
     return 0;
}

/* flush the next n things in buffer */
int flush_to_db(sqlite3 *db, sqlite3_stmt *nuc_stmt, sqlite3_stmt *cod_stmt,
                int32_t pos, struct buffer *buf, int32_t n) 
{
     for (; n > 0; n--, pos++) {
          beg_to_db(db, nuc_stmt, cod_stmt, pos, buf);
          buffer_next(buf);
     }
     return 0;
}

int main(int argc, char *argv[])
{
     int i;

     char *progname;
     char *dbname;
     char *bamfilename;
     char *animal;
     int day;

     sqlite3 *db = NULL;
     sqlite3_stmt *nuc_stmt, *cod_stmt;
     int dbopencode = 0;
     char *errormessage = NULL;
     int proper_pairs;

     samfile_t *samin;

     int32_t n_chr;
     sqlite3_int64 *db_chrids = NULL;
     sqlite3_int64 animal_id;

     bam1_t read;
     int32_t cur_tid, cur_pos, read_len;
     struct buffer buf;
     memset(&read, 0, sizeof(bam1_t));

     progname = *argv;
     argv++; argc--;

     proper_pairs = 0;
     while (argc > 0) {
          if ('-' == **argv) {
               switch (*(*argv + 1)) {
               case 'p':
                    proper_pairs = 1;
                    break;
               default:
                    fprintf(stderr, "Unknown option: %c\n", *(*argv + 1));
                    break;
               }
               argv++; argc--;
          }
          else {
               break;
          }
     }

     if (argc < 4) {
          printf("Usage: %s [-p] database bamfile animal day\n", progname);
          exit(ARG_ERROR);
     }
     else {
          fprintf(stderr, "DB name: %s... ", argv[0]);
          dbname = argv[0];
          bamfilename = argv[1];
          animal = argv[2];
          day = strtoul(argv[3], NULL, 10);
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
                        &nuc_stmt, NULL);
     for (i = 0; i < n_chr; i++) {
          sqlite3_bind_text(nuc_stmt, 1, samin->header->target_name[i], -1,
                            SQLITE_STATIC);
          sqlite3_step_onerow(nuc_stmt);
          db_chrids[i] = sqlite3_column_int64(nuc_stmt, 0);
          sqlite3_reset(nuc_stmt);
     }
     sqlite3_finalize(nuc_stmt);

     /* get id of animal */
     sqlite3_prepare_v2(db, sql_select_animalid, sizeof(sql_select_animalid),
                        &nuc_stmt, NULL);
     sqlite3_bind_text(nuc_stmt, 1, animal, -1, SQLITE_STATIC);
     sqlite3_step_onerow(nuc_stmt);
     animal_id = sqlite3_column_int64(nuc_stmt, 0);
     sqlite3_finalize(nuc_stmt);

     /* increases speed of insertion (means if the program crashes the db is */
     /* left invalid) */
     sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errormessage);

     /* start prepared statement for nucleotide insert */
     sqlite3_prepare_v2(db, sql_nuc_insert, sizeof(sql_nuc_insert),
                        &nuc_stmt, NULL);
     sqlite3_bind_int64(nuc_stmt, 1, animal_id);
     sqlite3_bind_int(nuc_stmt, 2, day);
     /* start prepared statement for codon insert */
     sqlite3_prepare_v2(db, sql_cod_insert, sizeof(sql_cod_insert),
                        &cod_stmt, NULL);
     sqlite3_bind_int64(cod_stmt, 1, animal_id);
     sqlite3_bind_int(cod_stmt, 2, day);

     sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &errormessage);
     sql_error(&errormessage);

     cur_tid = cur_pos = 0;
     sqlite3_bind_int64(nuc_stmt, 3, db_chrids[cur_tid]);
     sqlite3_bind_int64(cod_stmt, 3, db_chrids[cur_tid]);
     printf("chr %d\n", cur_tid);
     init_buffer(&buf);

     read_len = 0;
     while(samread(samin, &read) > 0) {
          if (read.core.flag & BAM_FUNMAP) {
               continue;
          }
          if (proper_pairs && !(read.core.flag & BAM_FPROPER_PAIR)) {
               continue;
          }
          if (cur_tid < read.core.tid) {
               /* flush buffer */
               flush_to_db(db, nuc_stmt, cod_stmt, cur_pos, &buf,
                           samin->header->target_len[cur_tid] - cur_pos);
               cur_tid = read.core.tid;
               sqlite3_bind_int64(nuc_stmt, 3, db_chrids[cur_tid]);
               sqlite3_bind_int64(cod_stmt, 3, db_chrids[cur_tid]);
               printf("chr %d\n", cur_tid);
               cur_pos = read.core.pos;
          }
          while (cur_pos < read.core.pos) {
               beg_to_db(db, nuc_stmt, cod_stmt, cur_pos, &buf);
               cur_pos++;
               buffer_next(&buf);
          }
          read_len = bam_cigar2qlen(&read.core, bam1_cigar(&read));
          if (read_len > buf.size) {
               resize_buffer(&buf, read_len);
          }
          buffer_insert_read(&buf, &read);
     }
     flush_to_db(db, nuc_stmt, cod_stmt, cur_pos, &buf,
                 samin->header->target_len[cur_tid] - cur_pos);

     sqlite3_exec(db, "COMMIT TRANSACTION;", NULL, NULL, &errormessage);
     sql_error(&errormessage);
     sqlite3_finalize(nuc_stmt);
     sqlite3_finalize(cod_stmt);

     free_buffer(&buf);
     free(db_chrids);
     samclose(samin);
     sqlite3_close(db);

     return 0;
}
