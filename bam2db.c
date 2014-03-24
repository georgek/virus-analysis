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
Ar int, Cr int, Gr int, Tr int, intD int)"

#define TO_NEXT_TAB while (getchar() != '\t')
#define TO_NEXT_NL while (getchar() != '\n')

static sqlite3 *db = NULL;
static sqlite3_stmt *stmt;
static char* dbname, *animal;
static int day;
static int chr_l = 21;
static char chr[21];
static int nerrors = 0, maxerrors = 10;

inline int char2base(char c)
{
     return ((c >> 1) & 0x3) | ((c & 0x20) >> 3);
}

inline int numberp(char c)
{
     return c >= '0' && c <= '9';
}

void parse_string()
{
     char c;
     int i;
     for (c = getchar(), i = 0; c != '\t' && i < chr_l-1; c = getchar(), ++i) {
          chr[i] = c;
     }
     chr[i] = '\0';
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

/* just reads over indel and ignores it */
void do_indel(char type)
{
     int length = parse_number();
     for(; length > 0; --length) {
          getchar();
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

void insert_row(int position, int (*nucs)[8], int dels)
{
     sqlite3_bind_text(stmt, 1, animal, strlen(animal), SQLITE_STATIC);
     sqlite3_bind_int(stmt, 2, day);
     sqlite3_bind_text(stmt, 3, chr, strlen(chr), SQLITE_STATIC);
     sqlite3_bind_int(stmt, 4, position);
     sqlite3_bind_int(stmt, 5, (*nucs)[char2base('A')]);
     sqlite3_bind_int(stmt, 6, (*nucs)[char2base('C')]);
     sqlite3_bind_int(stmt, 7, (*nucs)[char2base('G')]);
     sqlite3_bind_int(stmt, 8, (*nucs)[char2base('T')]);
     sqlite3_bind_int(stmt, 9, (*nucs)[char2base('a')]);
     sqlite3_bind_int(stmt, 10, (*nucs)[char2base('c')]);
     sqlite3_bind_int(stmt, 11, (*nucs)[char2base('g')]);
     sqlite3_bind_int(stmt, 12, (*nucs)[char2base('t')]);
     sqlite3_bind_int(stmt, 13, (*nucs)[char2base('D')]);

     if (sqlite3_step(stmt) != SQLITE_DONE) {
          fprintf(stderr, "Commit Failed!\n");
          exit(1);
     }

     sqlite3_reset(stmt);
}

void parse_pileup_line()
{
     char c;
     int position;
     char ref_base;
     int nucs[8] = {0};
     int dels = 0;

     parse_string();
     TO_NEXT_TAB;
     position = parse_number();
     TO_NEXT_TAB;
     ref_base = getchar();
     TO_NEXT_TAB;
     TO_NEXT_TAB;               /* number of pileups */

     /* now read main pileup to tab (qualities follow) */
     while((c = getchar())) {
          switch (c) {
          case '^':             /* read over quality */
               getchar();
               break;
          case 'a': case 'A':
          case 'c': case 'C':
          case 'g': case 'G':
          case 't': case 'T':
               nucs[char2base(c)]++;
               break;
          case '.':             /* forward */
               if ('N' == ref_base) {
                    fprintf(stderr, "Warning: N ref base and dot\n");
               }
               else{
                    nucs[(char2base(ref_base) & 0x3)]++;
               }
               break;
          case ',':             /* reverse */
               if ('N' == ref_base) {
                    fprintf(stderr, "Warning: N ref base and comma\n");
               }
               else {
                    nucs[(char2base(ref_base) | 0x4)]++;
               }
               break;
          case '+': case '-':
               do_indel(c);
               break;
          case '*':
               dels++;
               break;
          case '\t':
               TO_NEXT_NL;
               ungetc('\n', stdin);
               break;
          case '\n':
               insert_row(position, &nucs, dels);
               return;
          default:
               break;
          }
     }
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
     if (errormessage) {
          printf("SQLite3 error: %s\n", errormessage);
          sqlite3_free(errormessage);
          errormessage = NULL;
     }

     sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA count_changes=OFF", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errormessage);
     sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL, &errormessage);

     sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errormessage);
     do_error(&errormessage);

     /* init prepared statement */
     sqlite3_prepare_v2(db, sql, strlen(sql), &stmt, NULL);

     parse_pileup();

     sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errormessage);
     do_error(&errormessage);
     sqlite3_finalize(stmt);

     sqlite3_close(db);

     return 0;
}

