#!/bin/bash

dir=$1
dir=${dir%/*}
conda activate rna-seq
ls $dir/fq/*gz | xargs fastqc -t 8 -o $dir/qc
conda deactivate
multiqc $dir/qc -o $dir/qc
conda deactivate