#!/bin/bash

# assumes the reads are paired and filenames for pairs end in _1.fq and _2.fq

source bowtie-2.0.6
source samtools-0.1.18

mkdir -p bowtie2_output

stubs=$(ls -1 ../reads_autofs/*.fq | \
  awk '{print substr($0, 0, length($0)-5)}' | uniq)

for stub in $stubs; do
  bsub -J $(basename $stub) -o bowtie2_output/$(basename $stub) \
    -q Test128 -n 8 -R"rusage[mem=1024]" \
    bowtie2 -p 8 --fr -x ../refbt2.2.0/A_Victoria_3_75 \
    -1 ${stub}_1.fq -2 ${stub}_2.fq \
    | samtools view -bS - | samtools sort - $(basename $stub)
done

mkdir -p samtools_output

for bam in *.bam; do
  bsub -J "i $(basename $bam .bam)" -o samtools_output/$(basename $bam .bam) \
    -q Test128 -n 1 \
    samtools index $bam
done
