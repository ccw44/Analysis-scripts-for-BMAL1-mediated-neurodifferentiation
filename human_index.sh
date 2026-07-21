#! /bin/bash


tool=$1

conda activate m6a

cd ~/refseq/genome/human/GRCh38.111
if [ ${tool} == "hisat" ];then
  hisat2-build -p 8 Homo_sapiens.GRCh38.dna.primary_assembly.fa  Homo_sapiens.GRCh38.dna.primary_assembly
  mkdir index
  mv Homo_sapiens.GRCh38.dna.primary_assembly*ht2 index
  # mv GRCh38.primary_assembly* index
elif [ ${tool} == "STAR" ];then
  mkdir index_STAR
	STAR  --runMode genomeGenerate \
	--runThreadN 20 \
	--genomeDir  index_STAR \
	--genomeFastaFiles GRCh38.primary_assembly.genome.fa \
	--sjdbGTFfile gencode.v47.annotation.gtf \
	--sjdbOverhang 149
fi