#! /bin/bash


dir1=$1
# OUTPUT
# ~/CUTtag/motif/
dir1=${dir1%/*}
 
dir2=$2
# INPUT
dir2=${dir2%/*}

index=$3

ref_hsa=~/refseq/genome/human/GRCh38.111/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ref_mmu=~/refseq/genome/mouse/GRCm39.111/Mus_musculus.GRCm39.dna.primary_assembly.fa

if [ "$index" == "human" ]; then
    ref=$ref_hsa
elif [ "$index" == "mouse" ]; then
    ref=$ref_mmu
else
    echo "Invalid index. Please specify 'human' or 'mouse'."
    exit 1
fi

# ls ${dir2}/*.narrowPeak | while read id
ls ${dir2}/*.bed | while read id
do
	bedtools getfasta -fi ${ref} -bed ${id} -name -fo $dir1/$(basename $id).fa
done