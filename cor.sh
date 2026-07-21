#! /bin/bash


dir1=$1
# ~/m6a/cor/
dir1=${dir1%/*}

dir2=$2
# ~/m6a/map/
dir2=${dir2%/*}



multiBamSummary bins -bs 1000 -p 10 --bamfiles $(ls ${dir2}/*.bam | sort -V) -o ${dir1}/results.npz

plotCorrelation \
    -in ${dir1}/results.npz \
    --corMethod pearson \
    --skipZeros  \
    --plotNumbers \
    --labels k1-07-29-Input k2-08-01-Input sha1-07-29-Input sha2-08-01-Input k1-07-29-IP k2-08-01-IP sha1-07-29-IP sha2-08-01-IP \
    --plotTitle "m6A Correlation across samples" \
    --whatToPlot heatmap --colorMap Blues \
    -o ${dir1}/heatmap_pearsonCor.pearson.pdf --removeOutliers  \
    --outFileCorMatrix ${dir1}/pearsonCorr.rmOut.tab