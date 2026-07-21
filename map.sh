#! /bin/bash

dir1=$1
# ~/m6a/map/
# outdir=./map
dir1=${dir1%/*}

dir2=$2
# ~/m6a/bowtie2/
# indir=./bowtie2
dir2=${dir2%/*}

index=~/refseq/genome/human/GRCh38.111/index/Homo_sapiens.GRCh38.dna.primary_assembly

ls -1 $dir2/*gz |  rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | sort -n| uniq >$dir2/sampleid.txt

cat $dir2/sampleid.txt | while read id
do
  hisat2 -p 10 -x ${index} -1 ${dir2}/${id}.derRNA.fq.1.gz -2 ${dir2}/${id}.derRNA.fq.1.gz 2>$dir1/${id}.log | samtools sort -@ 10 -o ${dir1}/${id}.Hisat_aln.sorted.bam -  && samtools index ${dir1}/${id}.Hisat_aln.sorted.bam ${dir1}/${id}.Hisat_aln.sorted.bam.bai
done