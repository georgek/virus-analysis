#!/bin/bash
set -o nounset
set -o errexit

# assumes the reads are paired and filenames for pairs end in _1.fq and _2.fq

mkdir -p trimmed/trim-output

stubs=$(ls -1 reads_autofs/*.fq | awk '{print substr($0, 0, length($0)-5)}' | uniq)

for stub in $stubs; do
    bsub -J $(basename $stub) -o trimmed/trim-output/$(basename $stub) -q Test128 -R"rusage[mem=256]" "scripts/trim -p trimmed/$(basename $stub) reads_autofs/$(basename $stub)_1.fq reads_autofs/$(basename $stub)_2.fq"
done
