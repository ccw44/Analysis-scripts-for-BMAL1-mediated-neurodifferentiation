#!/bin/bash

dir1=$1
dir2=$2
dir1=${dir1%/*}
dir2=${dir2%/*}

conda activate m6a
if  [ ! -d "$dir2/cleanfq" ]; then
	# [space]before"!"
	mkdir $dir2/cleanfq 
else
	echo "exist $dir2/cleanfq"
fi
cd  $dir2/cleanfq 
#ls  $dir1/*.raw_1.fastq.gz > $dir2/cleanfq/1
#ls  $dir1/*.raw_2.fastq.gz > $dir2/cleanfq/2
# ls  $dir1/*_1.fastq.gz > $dir2/cleanfq/1
# ls  $dir1/*_2.fastq.gz > $dir2/cleanfq/2
#ls  $dir1/*.read1.fastq.gz > $dir2/cleanfq/1
#ls  $dir1/*.read2.fastq.gz > $dir2/cleanfq/2
ls $dir1/*_1.fq.gz > $dir2/cleanfq/1
ls $dir1/*_2.fq.gz > $dir2/cleanfq/2
#ls $dir1/*_1.fastq.gz > $dir2/cleanfq/1
#ls $dir1/*_2.fastq.gz > $dir2/cleanfq/2
paste 1 2  > $dir2/cleanfq/config
rm $dir2/cleanfq/{1,2}
config=$dir2/cleanfq/config

cat $config | while read id
do
        arr=(${id})
        fq1=${arr[0]}
        fq2=${arr[1]} 
trim_galore -q 25 --phred33 --length 36 --stringency 3 --paired -o $dir2/cleanfq $fq1 $fq2 
done 
rm $dir2/cleanfq/config

conda deactivate
#nohup bash ~/TEtranscripts/tecount.sh -d ~/RNAseq/map-STAR -o ~/TEtranscripts/result/tecount -s "C1 C2 C3 K1 K2 K3" > ~/TEtranscripts/DOR-tecount.log &
#nohup bash ~/RNAseq/cleanfq.sh ~/rawdata/ZY/00_Rawdata/SI/ ~/rawdata/ZY/ > ~/rawdata/ZY/clean1.log 2>&1 &