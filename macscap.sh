#! /bin/bash

dir1=$1
# ~/m6a/macscap/
dir1=${dir1%/*}

dir2=$2
# ~/m6a/map/
dir2=${dir2%/*}

cat group.txt | while read id1 id2 id3
do
  macs3 callpeak -t ${dir2}/${id2}.Hisat_aln.sorted.bam \
          -c ${dir2}/${id1}.Hisat_aln.sorted.bam \
          -n ${id3} \
	        --outdir ${dir1}\
          -B \
          -f BAMPE\
          --broad \
          -g hs \
          --broad-cutoff 0.1
          # --nomodel --extsize 150 \
          
done