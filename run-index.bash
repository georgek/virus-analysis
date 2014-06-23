#!/bin/bash

source samtools-0.1.18

mkdir -p samtools_output

for bam in *.bam; do
    bsub -J "idx $(basename $bam .bam)" -o samtools_output/$(basename $bam .bam) -q Prod128 -n 1 "mv $bam ${bam}.dup; samtools rmdup ${bam}.dup $bam; rm ${bam}.dup; samtools index $bam"
done
