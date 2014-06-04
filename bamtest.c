/* Counts the number of duplicates removed by looking at the DR field in each
 * read (if it exists). The DR field is updated with patched version of
 * bam_rmdup.c (samtools) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <bam/sam.h>
#include <bam/bam.h>

int main(int argc, char *argv[])
{
     char *filename = NULL;
     samfile_t *in;
     bam1_t *b;
     uint8_t *dr = 0;
     uint32_t dr_count = 0;

     b = bam_init1();

     if (argc < 2) {
          fprintf(stderr, "Usage: %s <bamfile>\n", argv[0]);
          exit(1);
     }
     else {
          filename = argv[1];
     }
     printf("filename: %s\n", filename);

     in = samopen(filename, "rb", NULL);
     if (!in) {
          fprintf(stderr, "Error opening file %s\n", filename);
     }

     while (samread(in, b) >= 0) {
          dr = bam_aux_get(b, "DR");
          if (dr) {
               dr_count += bam_aux2i(dr);
          }
     }

     samclose(in);
     bam_destroy1(b);

     printf("Duplicates removed: %u\n", dr_count);

     return 0;
}
