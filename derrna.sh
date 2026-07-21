#! /bin/bash

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/data/t020507/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/data/t020507/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/data/t020507/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/data/t020507/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate m6a

dir1=$1
# ~/m6a/bowtie2/
dir2=$2
# ~/m6a/cleanfq/

dir1=${dir1%/*}
dir2=${dir2%/*}
od=$dir1
index=~/refseq/rRNA_index/hg38_rRNA.rRNA

ls $dir2/*gz|cut -d"_" -f 1 |sort -u |while read id;
# ls $dir2/*gz|cut -d"." -f 1 |sort -u |while read id;
do
   id2=${id##*/};
   # bowtie2 -x $index --un-conc-gz ${od}/${id2}.derRNA.fq.gz -1 ${id}_1_val_1.fq.gz -2 ${id}_2_val_2.fq.gz -p 8 -S ${od}/${id2}.rRNA.sam 2> ${od}/${id2}.log;
   bowtie2 -x $index --un-conc-gz ${od}/${id2}.derRNA.fq.gz -1 ${id}_L1_1_val_1.fq.gz -2 ${id}_L1_2_val_2.fq.gz -p 8 -S ${od}/${id2}.rRNA.sam 2> ${od}/${id2}.log;
   # bowtie2 -x $index --un-conc-gz ${od}/${id2}.derRNA.fq.gz -1 ${id}.read1_Clean.fastq.gz -2 ${id}.read2_Clean.fastq.gz -p 8 -S ${od}/${id2}.rRNA.sam 2> ${od}/${id2}.log;
   rm ${od}/${id2}.rRNA.sam
done

conda deactivate 