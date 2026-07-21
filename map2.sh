#! /bin/bash

dir1=$1
dir1=${dir1%/*}

dir2=$2
dir2=${dir2%/*}

tool=$3

if [ ! -d "$dir1" ]; then
    mkdir -p $dir1 
else
    echo "exist $dir1"
fi
 
index=~/refseq/genome/human/GRCh38.111/index/Homo_sapiens.GRCh38.dna.primary_assembly
index2=~/refseq/genome/human/GRCh38.111/index_STAR/
 
# ls -1 $dir2/*gz | rev | cut -d'/' -f1 | rev | cut -d'_' -f1 | sort | uniq > $dir2/sampleid.txt
ls -1 $dir2/*gz |  rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | sort -n| uniq >$dir2/sampleid.txt
 
if [ ${tool} == "hisat2" ]; then
    cat ${dir2}/sampleid.txt | while read id
    do
        hisat2 -p 10 -x ${index} \
            -1 ${dir2}/${id}_R1_val_1.fq.gz \
            -2 ${dir2}/${id}_R2_val_2.fq.gz \
            2> ${dir1}/${id}.log | \
            samtools sort -@ 10 -o ${dir1}/${id}.Hisat_aln.sorted.bam - && \
            samtools index ${dir1}/${id}.Hisat_aln.sorted.bam
    done

elif [ ${tool} == "STAR" ]; then
    cat ${dir2}/sampleid.txt | while read id
    do
        STAR --runThreadN 30 \
            --genomeDir ${index2} \
            --readFilesCommand zcat \
            --readFilesIn ${dir2}/${id}.derRNA.fq.1.gz  ${dir2}/${id}.derRNA.fq.2.gz \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${dir1}/${id}.STAR.
            # --readFilesIn ${dir2}/${id}_L1_1_val_1.fq.gz ${dir2}/${id}_L1_2_val_2.fq.gz \
        samtools index ${dir1}/${id}.STAR.Aligned.sortedByCoord.out.bam
    done
fi