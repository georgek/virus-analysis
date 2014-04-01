#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include <bam/bam.h>
#include <bam/sam.h>

#define USAGE "you need to specify an input file\n"

struct aln_control_s {
     size_t pos;
     size_t cigar_len;
     uint32_t * cigar;
     char * seq;
     size_t seq_idx;
     char regular_print;
     char ins_print;
     char out[0x10000];
};

char translate_seqenc[] = {
       [0] = '='
     , [1] = 'A'
     , [2] = 'C'
     , [3] = 'M'
     , [4] = 'G'
     , [5] = 'R'
     , [6] = 'S'
     , [7] = 'V'
     , [8] = 'T'
     , [9] = 'W'
     , [10] = 'Y'
     , [11] = 'H'
     , [12] = 'K'
     , [13] = 'D'
     , [14] = 'B'
     , [15] = 'N'
};

char translate_seqeni[] = {
       [1] = 0
     , [2] = 1
     , [4] = 2
     , [8] = 3
};



size_t print_bam_seq(bam1_t * read) {
     uint8_t * seqc = bam1_seq(read);
     size_t length = read->core.l_qseq;
     size_t i = 0;

     while (length) {
          putchar(translate_seqenc[bam1_seqi(seqc, i)]);
          i++;
          length--;
     }
     return read->core.l_qseq;
}

void print_ccs(int ccs[0x4000][4], size_t lastreflen, char * readname) {
     char seq[lastreflen+1];
     int ref_pos;

     for (ref_pos = 0; ref_pos < lastreflen; ref_pos++) {
          char nucleotide = 'A';
          int nc = ccs[ref_pos][0];
          if (nc < ccs[ref_pos][1]) {
               nucleotide = 'C';
               nc = ccs[ref_pos][1];
          }
          if (nc < ccs[ref_pos][2]) {
               nucleotide = 'G';
               nc = ccs[ref_pos][2];
          }
          if (nc < ccs[ref_pos][3]) {
               nucleotide = 'T';
               nc = ccs[ref_pos][3];
          }
          if (nc)
               seq[ref_pos] = nucleotide;
          else
               seq[ref_pos] = 'N';
     }
     seq[ref_pos] = '\0';

     // while(seq[ref_pos - 1] == 'N') ref_pos--;
     // while(seq[0] == 'N') seq++;

     printf(">%s/ccs\n%.*s\n", readname, ref_pos, seq);
}

char cigar_translate[] = {
       [BAM_CMATCH] = 'M'
     , [BAM_CINS] = 'I'
     , [BAM_CDEL] = 'D'
     , [BAM_CREF_SKIP] = 'N'
     , [BAM_CSOFT_CLIP] = 'S'
     , [BAM_CHARD_CLIP] = 'H'
     , [BAM_CPAD] = 'P'
};

void printcigar(uint32_t * cigar, size_t cigarlen) {
     while (cigarlen--) {
          printf("%d%c", (*cigar) >> 4, cigar_translate[*cigar & 0xF]);
          cigar++;
     }
     putchar('\t');
}

void print_ccs_sam(int ccs[0x4000][4], size_t lastreflen,
                   char * readname, char * refname) {
     char seq[lastreflen+1];
     uint32_t cigar[lastreflen];
     size_t cigarlen = 0;
     int seq_pos = 0;

     int started = 0;
     int pos = 0;

     int lastcigar = BAM_CMATCH;
     int lastcigarcount = 0;

     int deletions = 0;

     size_t ref_pos;
     for (ref_pos = 0; ref_pos < lastreflen; ref_pos++) {
          char nucleotide = 'A';
          int nc = ccs[ref_pos][0];
          if (nc < ccs[ref_pos][1]) {
               nucleotide = 'C';
               nc = ccs[ref_pos][1];
          }
          if (nc < ccs[ref_pos][2]) {
               nucleotide = 'G';
               nc = ccs[ref_pos][2];
          }
          if (nc < ccs[ref_pos][3]) {
               nucleotide = 'T';
               nc = ccs[ref_pos][3];
          }
          if (nc) {
               started = 1;
               seq[seq_pos++] = nucleotide;

               if (lastcigar != BAM_CMATCH) {
                    cigar[cigarlen++] = (lastcigarcount << 4) | lastcigar;
                    lastcigar = BAM_CMATCH;
                    lastcigarcount = 0;
               }
               lastcigarcount++;
          }
          else {
               if (!started) {
                    pos++;
               } else {
                    if (lastcigar != BAM_CDEL) {
                         cigar[cigarlen++] = (lastcigarcount << 4) | lastcigar;
                         lastcigar = BAM_CDEL;
                         lastcigarcount = 0;
                    }

                    lastcigarcount++;
                    deletions++;
               }
          }

     }

     if (lastcigar != BAM_CDEL) {
          cigar[cigarlen++] = (lastcigarcount << 4) | lastcigar;
     }

     if (!seq_pos || deletions > 15)
          return;

     printf(
          "%s/accs\t" // qname
          "%d\t" // flag
          "%s\t" // RNAME:
          "%d\t" // pos
          "255\t" // MAPQ
          , readname
          , 0x2
          , refname
          , pos + 1
          );

     printcigar(cigar, cigarlen);

     printf(
          "*\t" // RNEXT
          "*\t" // PNEXT
          "%d\t" // TLEN
          "%.*s\t" // seq
          "*\n" // QUAL
          , seq_pos
          , seq_pos
          , seq
          );

     // seq[ref_pos] = '\0';

     // while(seq[ref_pos - 1] == 'N') ref_pos--;
     // while(seq[0] == 'N') seq++;

     // printf(">%s/ccs\n%.*s\n", readname, ref_pos, seq);
}

#define QNAME_TAG "/accs"

void export_ccs_sam(samfile_t * samfile, int ccs[0x4000][4],
                    int32_t tid, char * qname) {
     char seq[samfile->header->target_len[tid] + 1];
     memset(seq, 0, sizeof(seq));
     uint32_t cigar[samfile->header->target_len[tid] + 1];
     size_t cigarlen = 0;
     int seq_pos = 0;
     int seq_offset = 0;
     int seq_suboffset = 4;

     int started = 0;
     int pos = 0;

     int lastcigar = BAM_CMATCH;
     int lastcigarcount = 0;

     int deletions = 0;

     size_t ref_pos;

     size_t longestdeletion = 0;

     for (ref_pos = 0; ref_pos <= samfile->header->target_len[tid]; ref_pos++) {
          char nucleotide = 1;
          int nc = ccs[ref_pos][0];
          if (nc < ccs[ref_pos][1]) {
               nucleotide = 2;
               nc = ccs[ref_pos][1];
          }
          if (nc < ccs[ref_pos][2]) {
               nucleotide = 4;
               nc = ccs[ref_pos][2];
          }
          if (nc < ccs[ref_pos][3]) {
               nucleotide = 8;
               nc = ccs[ref_pos][3];
          }
          if (nc) {
               started = 1;
               seq[seq_offset] |= nucleotide << seq_suboffset;
               if (seq_suboffset ^= 4)
                    seq_offset++;

               seq_pos++;

               if (lastcigar != BAM_CMATCH) {
                    cigar[cigarlen++] = (lastcigarcount << 4) | lastcigar;

                    if (lastcigar == BAM_CDEL) {
                         deletions += lastcigarcount;
                         if (longestdeletion < lastcigarcount)
                              longestdeletion = lastcigarcount;
                    }

                    lastcigar = BAM_CMATCH;
                    lastcigarcount = 0;
               }
               lastcigarcount++;
          }
          else {
               if (!started) {
                    pos++;
               } else {
                    if (lastcigar != BAM_CDEL) {
                         cigar[cigarlen++] = (lastcigarcount << 4) | lastcigar;
                         lastcigar = BAM_CDEL;
                         lastcigarcount = 0;
                    }

                    lastcigarcount++;
               }
          }

     }

     if (lastcigar != BAM_CDEL) {
          cigar[cigarlen++] = (lastcigarcount << 4) | lastcigar;
     }

     if (!seq_pos || deletions > (seq_pos >> 5) ||
         longestdeletion > (seq_pos >> 6))
          return;

     bam1_t out;

     memset(&out, 0, sizeof(out));

     out.core.tid = tid;
     out.core.pos = pos;
     out.core.bin = bam_reg2bin(pos, pos+seq_pos);
     out.core.qual = 0xFF;
     out.core.l_qname = strlen(qname) + sizeof(QNAME_TAG);
     out.core.flag = 0x2;
     out.core.n_cigar = cigarlen;
     out.core.l_qseq = seq_pos;
     out.core.mtid = -1;
     out.core.mtid = -1;
     // out.core.isize = 0;

     out.l_aux = 0;
     out.data_len = (cigarlen * sizeof(*cigar)) + strlen(qname) +
          1 + ((seq_pos+ 1) >> 1) + seq_pos;
     out.m_data = out.data_len;

     uint8_t data[out.m_data];
     memcpy(data, qname, strlen(qname));
     memcpy(data + strlen(qname), QNAME_TAG, sizeof(QNAME_TAG));
     memcpy(data + out.core.l_qname, cigar, cigarlen * sizeof(*cigar));
     memcpy(data + cigarlen * sizeof(*cigar) +
            out.core.l_qname, seq, ((seq_pos+ 1) >> 1));
     memset(data + cigarlen * sizeof(*cigar) +
          out.core.l_qname + ((seq_pos+ 1) >> 1), 0xFF, seq_pos);

     out.data = data;

     samwrite(samfile, &out);

     // seq[ref_pos] = '\0';

     // while(seq[ref_pos - 1] == 'N') ref_pos--;
     // while(seq[0] == 'N') seq++;

     // printf(">%s/ccs\n%.*s\n", readname, ref_pos, seq);
}

int main(int argc, char** argv) {

     /* int show_inserts = 0; */
     char * filename = "-";
     char * outfilename = "-";

     char ** arg = argv+1;

     char * input_mode = "rb";
     char * output_mode = "wh";

     while (*arg) {
          switch (**arg) {
          case '-':
               switch ((*arg)[1]) {
               case 'S':
                    input_mode = "r";
                    break;
               case 'b':
                    output_mode = "wb";
                    break;
               case 'u':
                    output_mode = "wbu";
                    break;
               case 'o':
                    outfilename = *(++arg);
                    break;
               case '\0':
                    filename = *arg;
                    break;
               default:
                    fprintf(stderr, "unknown option %s\n", *arg);
                    break;
               }
               break;
          default:
               filename = *arg;
               break;
          }
          arg++;
     }

     if (argc < 2) {
          fprintf(stderr, USAGE);
          return -1;
     }

     // tamFile fp = sam_open(argv[1]);

     samfile_t * samfile = samopen(filename, input_mode, NULL);

     if (!samfile) {
          fprintf(stderr, "Error in file %s\n", filename);
          perror("samopen");
          return -1;
     }
     bam_header_t * bam_header = samfile->header;
     samfile_t * outfile = samopen(outfilename, output_mode, bam_header);
     if (!samfile) {
          fprintf(stderr, "Error in file %s\n", outfilename);
          perror("samopen");
          return -1;
     }

     // int num_headerlines;
     // if (!(num_headerlines = sam_header_parse(bam_header))) {
     // 	fprintf(stderr, "error sam_header_parse\n");
     // 	return -1;
     // } else {
     // 	fprintf(stderr, "read %d header lines\n", num_headerlines);
     // }

     char readname[1024] = {'\0'};
     // bam1_t reads[0x100];
     bam1_t current_read;
     memset(&current_read, 0, sizeof(bam1_t));
     size_t num_reads;
     int ccs[0x4000][4];
     int lastrefid = -1;
     size_t lastreflen = 0;

     /* size_t idmaxlen = 0; */

     for (num_reads = 0; samread(samfile, &current_read) >= 0; num_reads++) {
          // fprintf(stderr, "read in read %s\n", bam1_qname(&current_read));

          // samwrite(outfile, &current_read);

          if (current_read.core.flag & BAM_FUNMAP) {
               // printf(">%s/u\n", bam1_qname(&current_read));
               // print_bam_seq(&current_read);
               // putchar('\n');
               continue;
          }

          if (strcmp(bam1_qname(&current_read), readname)) {
               if (*readname) {
                    // print_ccs(ccs, lastreflen, readname);
                    // print_ccs_sam(ccs, lastreflen, readname,
                    // bam_header->target_name[lastrefid]);

                    export_ccs_sam(outfile, ccs, lastrefid, readname);

               }

               strcpy(readname, bam1_qname(&current_read));
               memset(ccs, 0, sizeof(ccs));
               lastrefid = current_read.core.tid;
          }

          if (lastrefid != current_read.core.tid) {
               // printf(">%s/da\n", bam1_qname(&current_read));
               // print_bam_seq(&current_read);
               // putchar('\n');
               continue;
          }

          size_t pos = current_read.core.pos;
          size_t cigar_len = current_read.core.n_cigar;
          uint32_t * cigar_p = bam1_cigar(&current_read);
          size_t seq_idx = 0;

          if (cigar_len && ((cigar_p[0] & 0xF) == BAM_CSOFT_CLIP)) {
               seq_idx += cigar_p[0] >> 4;
               cigar_p++;
               cigar_len--;
          }

          size_t ref_pos = pos;
          size_t outlen;
          lastreflen = bam_header->target_len[current_read.core.tid];

          for (outlen = 0; ref_pos < lastreflen; outlen++) {

          repeat:


               if (!cigar_len)
                    break;

               uint32_t cigar = *cigar_p;

               switch(cigar & 0xF) {
               case BAM_CINS:
                    // aln_control->ins_print =
                    // translate_seqenc[bam1_seqi(aln_control->seq,
                    // aln_control->seq_idx)];
                    //
                    // insert = 1;
                    // break;
               case BAM_CSOFT_CLIP:
                    seq_idx += cigar >> 4;
                    cigar_p++;
                    cigar_len--;
                    goto repeat;
                    break;
               case BAM_CPAD:
               case BAM_CREF_SKIP:
               case BAM_CHARD_CLIP:
                    fprintf(stderr, "undefined CIGAR option %d\n", cigar & 0xF);
                    return -1;
                    break;
               case BAM_CDEL:
                    // aln_control->regular_print = ' ';
                    ref_pos += cigar >> 4;
                    cigar_p++;
                    cigar_len--;
                    break;
               case BAM_CMATCH:
               {
                    int i;
                    for (i = 0; i < (cigar >> 4); i++) {
                         switch (bam1_seqi(bam1_seq(&current_read),
                                           seq_idx + i)) {
                         case 1: ccs[ref_pos][0]++; break;
                         case 2: ccs[ref_pos][1]++; break;
                         case 4: ccs[ref_pos][2]++; break;
                         case 8: ccs[ref_pos][3]++; break;
                         default: break;
                         }
                         ref_pos++;
                    }

                    seq_idx += cigar >> 4;
                    cigar_p++;
                    cigar_len--;
               }

               break;
               default:
                    fprintf(stderr, "undefined CIGAR option %d with value %d\n",
                            cigar & 0xF, cigar >> 4);
                    return -1;
                    break;
               }


          }

          free(current_read.data);
          memset(&current_read, 0, sizeof(bam1_t));
     }
     if (*readname) {
          // print_ccs(ccs, lastreflen, readname);
          // print_ccs_sam(ccs, lastreflen, readname,
          // bam_header->target_name[lastrefid]);

          export_ccs_sam(outfile, ccs, lastrefid, readname);

     }

     samclose(samfile);
     samclose(outfile);

     fprintf(stderr, "read %zu reads\n", num_reads);

     return 0;

}
