#! /bin/bash
 
 
index=~/refseq/genome/human/GRCh38.111/
index=${index%/*}
 
dir1=$1
# ~/CUTtag/map/
dir1=${dir1%/*}
 
computeMatrix reference-point -p 8 \
-b 3000 -a 3000  \
-R ${index}/Homo_sapiens.GRCh38.111.bed \
--referencePoint TSS \
-S $(ls ${dir1}/bam/*.bw | sort -V) \
-o ${dir1}/All_Sample.TSS-body-TES_distribute.gz
 
plotHeatmap --dpi 300 --heatmapHeight 20 \
  --matrixFile ${dir1}/All_Sample.TSS-body-TES_distribute.gz \
  --outFileName ${dir1}/All_Sample.TSS-body-TES_Heatmap.pdf
# --outFileName ${dir1}/All_Sample.TSS-body-TES_Heatmap.png