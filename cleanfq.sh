#! /bin/bash

dir1=$1
dir2=$2
dir1=${dir1%/*}
dir2=${dir2%/*}

cd  $dir1/cleanfq

ls  $dir2/*.read1.fastq.gz > $dir1/cleanfq/1
ls  $dir2/*.read2.fastq.gz > $dir1/cleanfq/2

paste 1 2  > $dir1/cleanfq/config
rm $dir1/cleanfq/{1,2}
config=$dir1/cleanfq/config

cat $config | while read id
do
        arr=(${id})
        fq1=${arr[0]}
        fq2=${arr[1]} 
trim_galore -q 25 --phred33 --length 15 --stringency 3 --fastqc --max_n 3 --paired -o $dir1/cleanfq $fq1 $fq2 
done 
rm $dir/cleanfq/config

multiqc $dir1/cleanfq -o $dir1/cleanfq