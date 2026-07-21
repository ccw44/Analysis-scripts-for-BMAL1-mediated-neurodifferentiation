#! /bin/bash

dir1=$1
# ~/CUTtag/
dir2=$2
# ~/rawdata/cuttag/
dir1=${dir1%/*}
dir2=${dir2%/*}


cd  $dir1/cleanfq
ls  $dir2/*_1.fq.gz > $dir1/cleanfq/1
ls  $dir2/*_2.fq.gz > $dir1/cleanfq/2

paste 1 2  > $dir1/cleanfq/config
rm $dir1/cleanfq/{1,2}
config=$dir1/cleanfq/config

cat $config | while read id
do
        arr=(${id})
        fq1=${arr[0]}
        fq2=${arr[1]} 
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 -m 75 -o $dir1/cleanfq/$(basename $fq1 .fq.gz)_cut.fq.gz -p $dir1/cleanfq/$(basename $fq2 .fq.gz)_cut.fq.gz $fq1 $fq2
done

rm $dir1/cleanfq/config

multiqc $dir1/cleanfq -o $dir1/cleanfq
