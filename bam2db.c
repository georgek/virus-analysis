/* Populates pileup_temp table from pileup, reads pileup from
 * stdin. Afterwards a proper pileup and animals etc. tables can be made from
 * it */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sqlite3.h>

#define PILEUP_CREATE "create table if not exists pileup_temp(animal, day int, \
chromosome int, position int,\
Af int, Cf int, Gf int, Tf int, \
Ar int, Cr int, Gr int, Tr int, D int)"

#define TO_NEXT_TAB while (getchar() != '\t')
#define TO_NEXT_NL while (getchar() != '\n')

static sqlite3 *db = NULL;
static sqlite3_stmt *stmt;
static char* dbname, *animal;
static int day;
#define chr_l 30
static char chr[chr_l+1];
static int nerrors = 0, maxerrors = 10;
static int nwarns = 0, maxwarns = 10;

static inline int char2base(char c)
{
     return ((c >> 1) & 0x3) | ((c & 0x20) >> 3);
}

static inline int numberp(char c)
{
     return c >= '0' && c <= '9';
}

void parse_chr()
{
     char c;
     int i;
     for (c = getchar(), i = 0; c != '\t' && i < chr_l; c = getchar(), ++i) {
          chr[i] = c;
     }
     chr[i] = '\0';
     if (c != '\t' && nwarns < maxwarns) {
          fprintf(stderr, "Warning, chromosome name truncated to %s.\n", chr);
          nwarns++;
     }
     ungetc(c, stdin);
}

int parse_number()
{
     char c;
     int number = 0;

     for (c = getchar(); numberp(c); c = getchar()) {
          number *= 10;
          number += c & 0xF;
     }
     ungetc(c, stdin);
     return number;
}

/* statics for parse_pileup_line */
static unsigned char pileup_c;
static int position;
static char ref_base;
static int nucs[8] = {0};
static int dels = 0;

/* just reads over indel and ignores it */
void do_indel()
{
     int length = parse_number();
     for(; length > 0; --length) {
          getchar();
     }
}

void do_del()
{
     dels++;
}

void do_quality()
{
     getchar();                 /* read over quality */
}

void do_base()
{
     nucs[char2base(pileup_c)]++;
}

void do_dot()
{
     if ('N' == ref_base && nwarns < maxwarns) {
          fprintf(stderr, "Warning: N ref base and dot\n");
     }
     else{
          nucs[(char2base(ref_base) & 0x3)]++;
     }
}

void do_comma()
{
     if ('N' == ref_base && nwarns < maxwarns) {
          fprintf(stderr, "Warning: N ref base and comma\n");
     }
     else {
          nucs[(char2base(ref_base) | 0x4)]++;
     }
}

void do_tab()
{
     TO_NEXT_NL;
     ungetc('\n', stdin);
}

void do_nothing()
{
     return;
}

void insert_row()
{
     int res_code;

     sqlite3_bind_text(stmt, 3, chr, strlen(chr), SQLITE_STATIC);
     sqlite3_bind_int(stmt, 4, position);
     sqlite3_bind_int(stmt, 5, nucs[char2base('A')]);
     sqlite3_bind_int(stmt, 6, nucs[char2base('C')]);
     sqlite3_bind_int(stmt, 7, nucs[char2base('G')]);
     sqlite3_bind_int(stmt, 8, nucs[char2base('T')]);
     sqlite3_bind_int(stmt, 9, nucs[char2base('a')]);
     sqlite3_bind_int(stmt, 10, nucs[char2base('c')]);
     sqlite3_bind_int(stmt, 11, nucs[char2base('g')]);
     sqlite3_bind_int(stmt, 12, nucs[char2base('t')]);
     sqlite3_bind_int(stmt, 13, dels);

     res_code = sqlite3_step(stmt);
     if (res_code != SQLITE_DONE) {
          fprintf(stderr, "SQLite3 error: %d!\n", res_code);
          exit(1);
     }

     sqlite3_reset(stmt);
}

typedef void (*Pfun)(void);
static Pfun fun_table[8] = {do_nothing, /* 0 */
                            do_indel,   /* 1 */
                            do_del,     /* 2 */
                            do_quality, /* 3 */
                            do_base,    /* 4 */
                            do_dot,     /* 5 */
                            do_comma,   /* 6 */
                            do_tab};    /* 7 */

static unsigned char a2f[256] = {
  /* 0     1     2     3     4     5     6     7           */
  /* 8     9     A     B     C     D     E     F           */
     0,    0,    0,    0,    0,    0,    0,    0,    /* 00 */
     0,    7,    0,    0,    0,    0,    0,    0,
     0,    0,    0,    0,    0,    0,    0,    0,    /* 10 */
     0,    0,    0,    0,    0,    0,    0,    0,
     0,    0,    0,    0,    0,    0,    0,    0,    /* 20 */
     0,    0,    2,    1,    6,    1,    5,    0,
     0,    0,    0,    0,    0,    0,    0,    0,    /* 30 */
     0,    0,    0,    0,    0,    0,    0,    0,
     0,    4,    0,    4,    0,    0,    0,    4,    /* 40 */
     0,    0,    0,    0,    0,    0,    0,    0,
     0,    0,    0,    0,    4,    0,    0,    0,    /* 50 */
     0,    0,    0,    0,    0,    0,    3,    0,
     0,    4,    0,    4,    0,    0,    0,    4,    /* 60 */
     0,    0,    0,    0,    0,    0,    0,    0,
     0,    0,    0,    0,    4,    0,    0,    0,    /* 70 */
     0,    0,    0,    0,    0,    0,    0,    0
};

void parse_pileup_line()
{
     int i;
     for(i = 0; i < 8; ++i){
          nucs[i] = 0;
     }
     dels = 0;

     parse_chr();
     TO_NEXT_TAB;
     position = parse_number();
     TO_NEXT_TAB;
     ref_base = getchar();
     TO_NEXT_TAB;
     TO_NEXT_TAB;               /* number of pileups */

     /* now read main pileup to tab (qualities follow) */
     while((pileup_c = getchar()) != '\n') {
          fun_table[a2f[pileup_c]]();
     }
     insert_row();
}

void parse_pileup()
{
     char nextchar;
     int nlines = 0;
     while((nextchar = getc(stdin)) != EOF) {
          nlines++;
          ungetc(nextchar, stdin);
          parse_pileup_line();
     }
}

void do_error(char **errormessage)
{
     if (*errormessage) {
          fprintf(stderr, "SQLite3 error: %s\n", *errormessage);
          sqlite3_free(*errormessage);
          *errormessage = NULL;
          nerrors++;
          if (nerrors >= maxerrors) {
               fprintf(stderr, "(Stopping after %d errors.)\n", maxerrors);
               exit(1);
          }
     }
}

int main(int argc, char *argv[])
{
     int dbopencode = 0;
     char* errormessage = NULL;
     char sql[] = "insert into pileup_temp values(?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)";

     if (argc < 4) {
          printf("Usage: bam2db database animal day\n");
          exit(1);
     }
     else {
          printf("DB name: %s\n", argv[1]);
          dbname = argv[1];
          animal = argv[2];
          day = strtoul(argv[3], NULL, 10);
     }

     dbopencode = sqlite3_open(argv[1], &db);
     if (SQLITE_OK == dbopencode) {
          printf("OK\n");
     }
     else {
          printf("SQLite3 error: %d\n", dbopencode);
          exit(1);
     }

     sqlite3_exec(db, PILEUP_CREATE, NULL, NULL, &errormessage);
     do_error(&errormessage);

     sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL, &errormessage);
     /* sqlite3_exec(db, "PRAGMA count_changes=OFF", NULL, NULL, &errormessage); */
     sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL, &errormessage);

     sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errormessage);
     do_error(&errormessage);

     /* init prepared statement */
     sqlite3_prepare_v2(db, sql, strlen(sql), &stmt, NULL);

     /* these are the same for entire transaction */
     sqlite3_bind_text(stmt, 1, animal, strlen(animal), SQLITE_STATIC);
     sqlite3_bind_int(stmt, 2, day);

     parse_pileup();

     sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errormessage);
     do_error(&errormessage);
     sqlite3_finalize(stmt);

     sqlite3_close(db);

     return 0;
}

