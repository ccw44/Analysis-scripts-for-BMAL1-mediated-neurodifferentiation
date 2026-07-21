#! /bin/bash
 

dir1=$1
# outdir=./map
dir1=${dir1%/*}

dir2=$2
# ~/CUTtag/cleanfq/
dir2=${dir2%/*}

index=GRCh38

cores=8

ls -1 $dir2/*gz | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2 | sort -n | uniq >$dir2/sampleid.txt

cat $dir2/sampleid.txt | while read id
do
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${index} -S ${dir1}/sam/${id}_bowtie2.sam -1 ${dir2}/${id}_1_cut.fq.gz -2 ${dir2}/${id}_2_cut.fq.gz &> ${dir1}/sam/bowtie2_summary/${id}_bowtie2.txt
done

# # ls ${dir1}/bam | while read id
# do
#   sam=${id##*/};sam1=${sam%%.*};
#   bamCoverage -p 6 -b ${id} --centerReads --ignoreDuplicates -o ${dir1}/bam/${sam1}.bw
# done