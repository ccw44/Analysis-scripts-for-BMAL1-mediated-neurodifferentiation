#! /bin/bash


dir1=$1
# ~/CUTtag/motif/meme_res/
dir1=${dir1%/*}
 
dir2=$2
# ~/CUTtag/motif/
# ~/CUTtag/motif/region_fa/
dir2=${dir2%/*}

index=~/refseq/meme_ref/motif_databases/
index=${index%/*}

ref1=$3

ref_hsa=HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
ref_mmu=MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme

if [ "$ref1" == "human" ]; then
    ref=$ref_hsa
elif [ "$ref1" == "mouse" ]; then
    ref=$ref_mmu
else
    echo "Invalid index. Please specify 'human' or 'mouse'."
    exit 1
fi

conda activate motif

ls ${dir2}/*fa |while read id;
do 
	meme-chip -meme-p 10 \
	-oc ${dir1}/$(basename ${id}).results/ \
	-db ${index}/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme \
	-db ${index}/${ref} \
	$id -meme-nmotifs 3
done