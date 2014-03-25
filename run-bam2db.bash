#!/bin/bash

# the bam files must be named like <animal>-d<day>-aln.bam

prog=$1                       # the bam2db program
db=$2                           # the database file
dir=$3                          # the diectory with bams in it

samtoolscmd="samtools mpileup -BQ0 -d10000000 "

main=$(ls -1 ${dir}/*.bam | xargs -n1 basename | sed -nr '/^F[0-9]+_d[0-9]+_aln\.bam$/p')
misc=$(ls -1 ${dir}/*.bam | xargs -n1 basename | sed -r  '/^F[0-9]+_d[0-9]+_aln\.bam$/d')

animaldays=$(echo "$main" | sed -r 's/^F([0-9]+)_d([0-9]+)_aln\.bam$/\1:\2/g' | uniq)
miscnames=$(echo "$misc" | sed -r 's/^(.*)_aln\.bam$/\1/g' | uniq)

for animalday in $animaldays; do
    a=${animalday%:*}
    d=${animalday#*:}
    echo "F${a}_d${d}"
    $samtoolscmd ${dir}/F${a}_d${d}_aln.bam | ${prog} ${db} F$a $d
done

for name in $miscnames; do
    echo "${name}"
    $samtoolscmd ${dir}/${name}_aln.bam | ${prog} ${db} $name 0
done

