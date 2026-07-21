#! /bin/bash
 

dir1=$1
# ~/RNAseq/map/
# outdir=./map
dir1=${dir1%/*}
 
conda activate m6a

for bam_file in ${dir1}/*.STAR.Aligned.sortedByCoord.out.bam; do
    if [ -f "$bam_file" ]; then
        samtools index ${bam_file}
        echo "Finished processing $bam_file -> ${bam_file%.bam}"
    else
        echo "No bam files found in the current directory."
        break
    fi
done

wait
echo "All indexing jobs completed."