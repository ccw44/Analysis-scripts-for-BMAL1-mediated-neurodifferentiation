#! /bin/bash


dir_bam=$1
dir_bam=${dir_bam%/*}

dir_out=$2
dir_out=${dir_out%/*}

tool=$3

if [ ! -d "$dir_out" ]; then
    mkdir -p $dir_out 
    echo "Created output directory: $dir_out"
else
    echo "Directory exists: $dir_out"
fi
 
gtf=~/refseq/genome/human/GRCh38.111/Homo_sapiens.GRCh38.111.gtf
 
if [ ${tool} == "hisat2" ]; then
    ls -1 $dir_bam/*.Hisat_aln.sorted.bam | rev | cut -d'/' -f1 | rev | sed 's/.Hisat_aln.sorted.bam//g' > $dir_out/sampleid.txt
    bam_suffix=".Hisat_aln.sorted.bam"
elif [ ${tool} == "STAR" ]; then
    ls -1 $dir_bam/*.STAR.Aligned.sortedByCoord.out.bam | rev | cut -d'/' -f1 | rev | sed 's/.STAR.Aligned.sortedByCoord.out.bam//g' > $dir_out/sampleid.txt
    bam_suffix=".STAR.Aligned.sortedByCoord.out.bam"
else
    echo "Error: tool must be hisat2 or STAR"
    exit 1
fi

echo "========================================"
featureCounts -T 20 -p \
    -a ${gtf} \
    -t exon -g gene_id \
    --extraAttributes gene_name \
    -o ${dir_out}/gene_counts_raw.txt \
    ${dir_bam}/*${bam_suffix}

sed '1d' ${dir_out}/gene_counts_raw.txt | cut -f1,2,7- | tr '\t' ',' > ${dir_out}/Gene_Counts_Clean.csv
echo "Gene counts output to: ${dir_out}/Gene_Counts_Clean.csv"

featureCounts -T 20 -p \
    -a ${gtf} \
    -t exon -g transcript_id \
    --extraAttributes gene_id,gene_name \
    -o ${dir_out}/transcript_counts_raw.txt \
    ${dir_bam}/*${bam_suffix}

sed '1d' ${dir_out}/transcript_counts_raw.txt | cut -f1,2,7- | tr '\t' ',' > ${dir_out}/transcript_counts_clean.csv
echo "transcript counts output to: ${dir_out}/transcript_Counts_Clean.csv"