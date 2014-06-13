#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include <unistd.h>

#include <bam/bam.h>
#include <bam/sam.h>

#include "errors.h"

/* int do_read(const bam1_t *read, void *d) */
/* { */
/*      /\* TODO: check CIGAR strings to find real read_pos *\/ */
/*      PosData *pos = d; */
/*      /\* int32_t read_pos = pos->genome_position - read->core.pos; *\/ */
/*      int32_t read_pos = 0, gen_pos = read->core.pos, diff; */
/*      int i; */
/*      int base = bam1_seqi(bam1_seq(read), read_pos); */
/*      uint32_t *nucs; */

/*      if (!(read->core.flag & BAM_FPROPER_PAIR)) { */
/*           return 0; */
/*      } */

/*      diff = pos->genome_position - gen_pos; */
/*      for (i = 0; (i < read->core.n_cigar) && diff > 0; ++i) { */
/*           switch (bam_cigar_op(bam1_cigar(read)[i])) { */
/*           case BAM_CMATCH: */
               
/*                break; */
/*           case BAM_CINS: */

/*                break; */
/*           case BAM_CDEL: */

/*                break; */
/*           default: */
/*                break; */
/*           } */
/*      } */

/*      nucs = (read->core.flag & BAM_FREVERSE) ? */
/*           pos->read_positions[read_pos].reverse : */
/*           pos->read_positions[read_pos].forward; */
/*      if (read->core.flag & BAM_FREVERSE) { */
/*           pos->n_reads_r++; */
/*      } */
/*      else { */
/*           pos->n_reads_f++; */
/*      } */

/*      if (read_pos < pos->read_length) { */
/*           switch (base) { */
/*           case 1: */
/*                nucs[0]++; */
/*                break; */
/*           case 2: */
/*                nucs[1]++; */
/*                break; */
/*           case 4: */
/*                nucs[2]++; */
/*                break; */
/*           case 8: */
/*                nucs[3]++; */
/*                break; */
/*           } */
/*      } */
/*      else { */
/*           printf("pos: %d, len: %d\n", read->core.pos, read->core.l_qseq); */
/*      } */

/*      return 0; */
/* } */

typedef struct pos_nucs {
     uint32_t forward[4];
     uint32_t reverse[4];
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
     int i;
     uint32_t *nucs;

     if (p->genome_position == pos) {
          for (i = 0; i < n; ++i) {
               if (!pl[i].is_del) {
                    if (pl[i].b->core.flag & BAM_FREVERSE) {
                         p->n_reads_r++;
                         nucs = p->read_positions[pl[i].qpos].reverse;
                    }
                    else {
                         p->n_reads_f++;
                         nucs = p->read_positions[pl[i].qpos].forward;
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

int main(int argc, char *argv[])
{
     size_t i;
     char *bamfilename;
     int32_t tid;

     samfile_t *bamin;
     bam_index_t *bamidx;
     bam_plbuf_t *buf;
     bam1_t *bam_read;

     PosData pos;

     if (argc < 5) {
          printf("Usage: %s bam_file read_length chromosome_id position\n",
                 argv[0]);
          exit(ARG_ERROR);
     }
     else {
          bamfilename = argv[1];
          pos.read_length = strtol(argv[2], NULL, 10);
          tid = strtol(argv[3], NULL, 10);
          pos.genome_position = strtol(argv[4], NULL, 10) - 1;
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
     /* for (i = 0; i < bamin->header->n_targets; ++i) { */
     /*      printf("%zu : %s\n", i, bamin->header->target_name[i]); */
     /* } */
     bam_read = bam_init1();

     pos.n_reads_f = 0;
     pos.n_reads_r = 0;
     pos.read_positions = calloc(pos.read_length, sizeof(PosNucs));

     buf = bam_plbuf_init(pileup_func, &pos);
     /* remove maximum pileup depth */
     bam_plp_set_maxcnt(buf->iter, INT_MAX);
     bam_fetch(bamin->x.bam, bamidx,
               tid, pos.genome_position, pos.genome_position+1,
               buf, &fetch_func);
     bam_plbuf_push(0, buf);    /* finish pileup */

     printf("Total reads: %zu\n", n_fetched);
     printf("Proper pairs: %zu\n", n_pushed);
     printf("Reads forward: %zu, reads reverse: %zu\n",
            pos.n_reads_f, pos.n_reads_r);
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

     bam_plbuf_destroy(buf);
     free(pos.read_positions);
     bam_destroy1(bam_read);
     bam_index_destroy(bamidx);
     samclose(bamin);
     return 0;
}
