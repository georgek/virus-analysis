#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include <bam/bam.h>
#include <bam/sam.h>

#include "errors.h"

typedef struct pos_nucs {
     uint32_t forward[4];
     uint32_t reverse[4];
     float qforward[4];
     float qreverse[4];
     float avg_rl_f[4];
     float avg_rl_r[4];
} PosNucs;

typedef struct pos_data {
     size_t n_reads_f, n_reads_r;
     int32_t genome_position;
     int32_t read_length;
     PosNucs *read_positions;
} PosData;

static size_t n_fetched = 0;
static size_t n_pushed = 0;

/* callback for bam_fetch() */
static int fetch_func(const bam1_t *b, void *data)
{
     bam_plbuf_t *buf = data;
     if (b->core.flag & BAM_FPROPER_PAIR) {
          bam_plbuf_push(b, buf);
          n_pushed++;
     }
     n_fetched++;
     return 0;
}

/* callback for bam_plbuf_init() */
static int pileup_func(uint32_t tid, uint32_t pos, int n,
                       const bam_pileup1_t *pl, void *data)
{
     PosData *p = data;
     int i, j;
     uint32_t *nucs;
     float *quals, *rl;

     if (p->genome_position == pos) {
          for (i = 0; i < n; ++i) {
               if (!pl[i].is_del) {
                    if (pl[i].b->core.flag & BAM_FREVERSE) {
                         p->n_reads_r++;
                         nucs = p->read_positions[pl[i].qpos].reverse;
                         quals = p->read_positions[pl[i].qpos].qreverse;
                         rl = p->read_positions[pl[i].qpos].avg_rl_r;
                    }
                    else {
                         p->n_reads_f++;
                         nucs = p->read_positions[pl[i].qpos].forward;
                         quals = p->read_positions[pl[i].qpos].qforward;
                         rl = p->read_positions[pl[i].qpos].avg_rl_f;
                    }
                    switch(bam1_seqi(bam1_seq(pl[i].b), pl[i].qpos)) {
                    case 1:
                         j = 0;
                         break;
                    case 2:
                         j = 1;
                         break;
                    case 4:
                         j = 2;
                         break;
                    case 8:
                         j = 3;
                         break;
                    default:
                         return 0;
                    }
                    nucs[j]++;
                    quals[j] *= ((float)nucs[j]-1)/nucs[j];
                    quals[j] += ((float)((char *)bam1_qual(pl[i].b))[j]-33)/nucs[j];
                    rl[j] *= ((float)nucs[j]-1)/nucs[j];
                    rl[j] += ((float)pl[i].b->core.l_qseq)/nucs[j];
               }
          }
     }

     return 0;
}

int main(int argc, char *argv[])
{
     char *progname;
     int verbose, ggplot, quality;

     size_t i, j;
     char basenames[] = {'a','c','g','t'};
     char *bamfilename;
     int32_t tid;

     samfile_t *bamin;
     bam_index_t *bamidx;
     bam_plbuf_t *buf;
     bam1_t *bam_read;

     PosData pos;

     progname = *argv;
     argv++; argc--;
     verbose = 0;
     ggplot = 0;
     quality = 0;
     while (argc > 0) {
          if ('-' == **argv) {
               switch (*(*argv + 1)) {
               case 'v':
                    verbose = 1;
                    break;
               case 'g':
                    ggplot = 1;
                    break;
               case 'q':
                    quality = 1;
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
          printf("Usage: %s [-v] [-g] [-q] bam_file read_length chromosome_id position\n",
                 progname);
          exit(ARG_ERROR);
     }
     else {
          bamfilename = argv[0];
          pos.read_length = strtol(argv[1], NULL, 10);
          tid = strtol(argv[2], NULL, 10);
          pos.genome_position = strtol(argv[3], NULL, 10) - 1;
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
     bam_read = bam_init1();

     pos.n_reads_f = 0;
     pos.n_reads_r = 0;
     pos.read_positions = calloc(pos.read_length, sizeof(PosNucs));

     buf = bam_plbuf_init(&pileup_func, &pos);
     /* disable maximum pileup depth */
     bam_plp_set_maxcnt(buf->iter, INT_MAX);
     bam_fetch(bamin->x.bam, bamidx,
               tid, pos.genome_position, pos.genome_position+1,
               buf, &fetch_func);
     bam_plbuf_push(0, buf);    /* finish pileup */

     if (verbose) {
          printf("Total reads: %zu\n", n_fetched);
          printf("Proper pairs: %zu\n", n_pushed);
          printf("Reads forward: %zu, reads reverse: %zu\n",
                 pos.n_reads_f, pos.n_reads_r);
     }

     if (ggplot) {
          printf("%5s %5s %6s %6s %6s\n", "pos", "base", "count", "qual", "avgrl");
          for (i = 0; i < pos.read_length; ++i) {
               for (j = 0; j < 4; ++j) {
                    printf("%5zu %4cf %6d %6.2f %6.2f\n",
                           i+1,
                           basenames[j],
                           pos.read_positions[i].forward[j],
                           pos.read_positions[i].qforward[j],
                           pos.read_positions[i].avg_rl_f[j]);
               }
               for (j = 0; j < 4; ++j) {
                    printf("%5zu %4cr %6d %6.2f %6.2f\n",
                           i+1,
                           basenames[j],
                           pos.read_positions[i].reverse[j],
                           pos.read_positions[i].qreverse[j],
                           pos.read_positions[i].avg_rl_r[j]);
               }
          }
     }
     else if (quality) {
          printf("%5s %6s %14s %14s %14s %14s %14s %14s %14s\n",
                 "pos", "af", "cf", "gf", "tf", "ar", "cr", "gr", "tr");
          for (i = 0; i < pos.read_length; ++i) {
               printf("%5zu %6d (%5.2f) %6d (%5.2f) %6d (%5.2f) %6d (%5.2f) %6d (%5.2f) %6d (%5.2f) %6d (%5.2f) %6d (%5.2f)\n",
                      i + 1,
                      pos.read_positions[i].forward[0],
                      pos.read_positions[i].qforward[0],
                      pos.read_positions[i].forward[1],
                      pos.read_positions[i].qforward[1],
                      pos.read_positions[i].forward[2],
                      pos.read_positions[i].qforward[2],
                      pos.read_positions[i].forward[3],
                      pos.read_positions[i].qforward[3],
                      pos.read_positions[i].reverse[0],
                      pos.read_positions[i].qreverse[0],
                      pos.read_positions[i].reverse[1],
                      pos.read_positions[i].qreverse[1],
                      pos.read_positions[i].reverse[2],
                      pos.read_positions[i].qreverse[2],
                      pos.read_positions[i].reverse[3],
                      pos.read_positions[i].qreverse[3]);
          }
     }
     else {
          printf("%5s %6s %6s %6s %6s %6s %6s %6s %6s\n",
                 "pos", "af", "cf", "gf", "tf", "ar", "cr", "gr", "tr");
          for (i = 0; i < pos.read_length; ++i) {
               printf("%5zu %6d %6d %6d %6d %6d %6d %6d %6d\n",
                      i + 1,
                      pos.read_positions[i].forward[0],
                      pos.read_positions[i].forward[1],
                      pos.read_positions[i].forward[2],
                      pos.read_positions[i].forward[3],
                      pos.read_positions[i].reverse[0],
                      pos.read_positions[i].reverse[1],
                      pos.read_positions[i].reverse[2],
                      pos.read_positions[i].reverse[3]);
          }
     }
     bam_plbuf_destroy(buf);
     free(pos.read_positions);
     bam_destroy1(bam_read);
     bam_index_destroy(bamidx);
     samclose(bamin);
     return 0;
}
