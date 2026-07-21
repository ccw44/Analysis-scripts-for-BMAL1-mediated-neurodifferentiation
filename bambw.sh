#! /bin/bash

dir1=$1
# ~/CUTtag/map/
dir1=${dir1%/*}

ls ${dir1}/sam/*.sam | while read id
do
  samid=$(basename ${id})
  samtools view -@ 10 -bS ${id}| samtools sort -@ 10 -o ${dir1}/bam/${samid}.sorted.bam -
  samtools index ${dir1}/bam/${samid}.sorted.bam
done

ls ${dir1}/bam/*bam | while read id
do
  sam=${id##*/};sam1=${sam%%.*};
  bamCoverage -p 8 -b ${id} --centerReads --ignoreDuplicates -o ${dir1}/bam/${sam1}.bw
done