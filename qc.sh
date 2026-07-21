#!/bin/bash

dir1=$1
dir2=$2
dir1=${dir1%/*}
dir2=${dir2%/*}

ls $dir2/*gz | xargs fastqc -t 20 -o $dir1/qc
multiqc $dir1/qc -o $dir1/qc

conda deactivate