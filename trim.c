#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>


#define handle_error(msg) \
     do { perror(msg); exit(EXIT_FAILURE); } while (0)

#define MAX_SEQLEN 0x400

void help(char* name) {
     fprintf(stderr,
             "Usage: %s [-p prefix] [-m minumum_length (60)] "
             "[-q phred_min (33)] [-t q_threshold (2)] "
             "filename [filename_pair]\n",
             name);
     fprintf(stderr,"\nEither filename can be set to - (stdin).\n");
}

int main (int argc, char** argv) {
     char *filename1 = 0;
     char *filename2 = 0;
     char *prefix = 0;
     int minlength = 60;
     char baseq = '!';           /* 33 */
     char threshold = 2;

     for (++argv;*argv;argv++) {
          if (**argv == '-') {
               switch((*argv)[1]) {
               case 'p': prefix = *(++argv); break;
               case 'm': minlength = strtol(*(++argv), NULL, 10); break;
               case 'h': help("trim_reads"); exit(0); break;
               case 'q': baseq = strtol(*(++argv), NULL, 10); break;
               case 't': threshold = strtol(*(++argv), NULL, 10); break;
               case '-': if (!memcmp((*argv) + 2,"help",4)) {
                         help("trim_reads");
                         exit(0);
                    } break;
               case 0: {if (!filename1)
                              filename1 = *argv;
                         else if (!filename2)
                              filename2 = *argv;
               } break;
               }
          } else if (!filename1)
               filename1 = *argv;
          else if (!filename2)
               filename2 = *argv;
          else {
               fprintf(stderr,"Unknown option %s\n",*argv);
               help("trim_reads");
               exit(1);
          }
     }
     if (baseq < 33) {
          fprintf(stderr, "Base quality must be at least 33.\n");
          exit(1);
     }

     FILE *in1,*in2 = 0;

     if (filename1 && *filename1 != '-') {
          in1 = fopen(filename1,"r");
          if (!in1)
               handle_error("fopen");
     } else {
          in1 = stdin;
          filename1 = "-";
     }

     if (filename2 && *filename2 != '-') {
          in2 = fopen(filename2,"r");
          if (!in2)
               handle_error("fopen");
     } else if (filename2 && *filename2 == '-')
          in2 = stdin;

     int i,i1,i2;
     /* int j; */

     char keptnames[MAX_SEQLEN];
     char keptname1[MAX_SEQLEN];
     char keptname2[MAX_SEQLEN];
     char discardedname[MAX_SEQLEN];
     if (!prefix) {
          prefix = malloc(MAX_SEQLEN);
          if (*filename1 != '-')
               memcpy(prefix,filename1,strlen(filename1)+1);
          else if (filename2 && *filename2 != '-') {
               memcpy(prefix,filename2,strlen(filename2)+1);
          } else {
               fprintf(stderr,"You must specify an output prefix "
                       "when reading from the pipe.\n");
               exit(1);
          }
          for (i = strlen(prefix)-1; i >= 0; i--) {
               if (prefix[i] == '.') {
                    prefix[i] = 0;
                    break;
               }
               if (prefix[i] == '/') {
                    if (*filename1 != '-')
                         memcpy(prefix,filename1,strlen(filename1)+1);
                    else if (filename2 && *filename2 != '-') {
                         memcpy(prefix,filename2,strlen(filename2)+1);
                    }
                    break;
               }
          }
          if (i < 0) {
               if (*filename1 != '-')
                    memcpy(prefix,filename1,strlen(filename1)+1);
               else if (filename2 && *filename2 != '-') {
                    memcpy(prefix,filename2,strlen(filename2)+1);
               }
          }
     }

     sprintf(keptnames,"%s_filtered.fq",prefix);
     sprintf(discardedname,"%s_discarded.fq",prefix);
     if (in2) {
          sprintf(keptname1,"%s_1_filtered.fq",prefix);
          sprintf(keptname2,"%s_2_filtered.fq",prefix);
     }

     //fprintf(stderr,"%s\n",prefix);

     FILE *keptfp1 = NULL, *keptfp2 = NULL;
     FILE *keptfps = fopen(keptnames,"w");
     if (in2) {
          keptfp1 = fopen(keptname1,"w");
          keptfp2 = fopen(keptname2,"w");
     }
     FILE* discardedfp = fopen(discardedname,"w");

     char id1[MAX_SEQLEN],sequence1[MAX_SEQLEN],sequencecpy1[MAX_SEQLEN],
          qid1[MAX_SEQLEN],quality1[MAX_SEQLEN],qualitycpy1[MAX_SEQLEN];
     char id2[MAX_SEQLEN],sequence2[MAX_SEQLEN],sequencecpy2[MAX_SEQLEN],
          qid2[MAX_SEQLEN],quality2[MAX_SEQLEN],qualitycpy2[MAX_SEQLEN];
     char* line1 = id1;
     char* line2 = id2;
     int stage = 0; // 0 for read id, 1 for sequence, 2 for qid, 3 for quality
     /* long reads = 0; */
     long kept = 0;
     long discarded = 0;
     long orphaned = 0;
     long lengthcounts[MAX_SEQLEN];
     memset(lengthcounts, 0, sizeof(lengthcounts));
     i2 = -1;
     while (fgets(line1,MAX_SEQLEN,in1) &&
            (!in2 || fgets(line2,MAX_SEQLEN,in2))) {
          switch(stage) {
          case 0:
               line1 = sequence1;
               line2 = sequence2;
               break;
          case 1:
               line1 = qid1;
               line2 = qid2;
               break;
          case 2:
               line1 = quality1;
               line2 = quality2;
               break;
          case 3:
               line1 = id1;
               line2 = id2;
               memcpy(sequencecpy1,sequence1,MAX_SEQLEN);
               memcpy(qualitycpy1,quality1,MAX_SEQLEN);
               if (in2) {
                    memcpy(sequencecpy2,sequence2,MAX_SEQLEN);
                    memcpy(qualitycpy2,quality2,MAX_SEQLEN);
               }
               for (i1 = strlen(quality1)-1; i1 >= 0; i1--) {
                    if (quality1[i1] <= (baseq + threshold)) { /* '\n' is less than baseq */
                         quality1[i1] = '\0';
                         sequence1[i1] = '\0';
                    } else
                         break;
               }
               if (in2) {
                    for (i2 = strlen(quality2)-1; i2 >= 0; i2--) {
                         if (quality2[i2] <= (baseq + threshold)) {
                              quality2[i2] = '\0';
                              sequence2[i2] = '\0';
                         } else
                              break;
                    }
               }
               if (i1 >= (minlength-1) && i2 >= (minlength-1)) {
                    lengthcounts[i1]++;
                    lengthcounts[i2]++;
                    fprintf(keptfp1,"%s%s\n%s%s\n",id1,sequence1,qid1,quality1);
                    fprintf(keptfp2,"%s%s\n%s%s\n",id2,sequence2,qid2,quality2);
                    kept+=2;
               } else {
                    if (i1 >= (minlength-1)) {
                         lengthcounts[i1]++;
                         fprintf(keptfps,"%s%s\n%s%s\n",
                                 id1,sequence1,qid1,quality1);
                         kept++;
                         orphaned++;
                    } else if (i2 >= (minlength-1)) {
                         lengthcounts[i2]++;
                         fprintf(keptfps,"%s%s\n%s%s\n",
                                 id2,sequence2,qid2,quality2);
                         kept++;
                         orphaned++;
                    }
                    if (i1 < (minlength-1)) {
                         fprintf(discardedfp,"%s%s%s%s",
                                 id1,sequencecpy1,qid1,qualitycpy1);
                         discarded++;
                    }
                    if (in2 && i2 < (minlength-1)) {
                         fprintf(discardedfp,"%s%s%s%s",
                                 id2,sequencecpy2,qid2,qualitycpy2);
                         discarded++;
                    }
               }
               // reads++;
               break;
          }
          stage++;
          stage %= 4;
     }

     fclose(keptfps);
     fclose(discardedfp);
     if (in2) {
          fclose(keptfp1);
          fclose(keptfp2);
     }

     fprintf(stderr,"%ld reads processed.\n%ld reads filtered, "
             "%ld reads discarded, %ld reads orphaned.\n",
             kept+discarded,kept,discarded, orphaned);

     int longest = MAX_SEQLEN - 1;

     while (longest > 0 && !lengthcounts[longest]) longest--;

     double avglen = 0;
     for (i = 1; i < MAX_SEQLEN; i++)
          avglen += lengthcounts[i] * i;

     avglen /= kept;

     printf("%s,%ld,%ld,%ld,%ld,%f,%s\n",
            filename1, kept+discarded,kept,discarded,
            orphaned,avglen, filename2);

     exit(0);
}
