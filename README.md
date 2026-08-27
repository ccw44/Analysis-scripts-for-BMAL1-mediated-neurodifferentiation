# Analysis-scripts-for-BMAL1-mediated-neurodifferentiation
This release contains scripts used for bioinformatics analysis in the study of BMAL1 regulation during neural differentiation.
# Transcriptomic, CUT&Tag, and m6A-seq Analysis of BMAL1-Mediated Neural Differentiation

This repository contains the computational scripts used for transcriptomic, CUT&Tag-seq, and m6A-seq analyses in the study investigating the role of the circadian regulator **BMAL1** in human embryonic stem cell (ESC) differentiation toward neural stem cells (NSCs).

The scripts enable reproducible analysis of:

- RNA sequencing (RNA-seq)
- BMAL1 CUT&Tag sequencing
- m6A methylation sequencing (MeRIP-seq/m6A-seq)
- Differential expression analysis
- Peak identification and annotation
- Functional enrichment analysis
- Multi-omics integration analyses

---

## Overview

BMAL1-mutant human ESCs and derived NSCs were profiled using transcriptomic and epitranscriptomic approaches to characterize BMAL1-dependent regulation during neural lineage commitment.

The sequencing datasets analyzed in this study include:

| Dataset | Cell type | Purpose |
|---|---|---|
| RNA-seq | ESCs and NSCs | Transcriptomic profiling and differential gene expression analysis |
| CUT&Tag-seq | BMAL1-3×Flag ESCs and NSCs | Identification of BMAL1 genomic binding sites |
| m6A-seq (MeRIP-seq) | NSCs | Characterization of m6A modification patterns |

---

# Analysis Workflow

## 1. RNA-seq Analysis

### Sequencing and preprocessing

RNA libraries were generated using the VAHTS mRNA-seq V2 Library Prep Kit and sequenced using Illumina NovaSeq platforms with paired-end reads.

Raw sequencing reads were processed using:

- **FastQC** for sequencing quality assessment
- **Trimmomatic** for adapter and low-quality read removal
- **HISAT2** for alignment to the human reference genome (GRCh38/hg38)

### Transcript quantification and differential expression

Gene expression quantification was performed using:

- **StringTie**

Sample relationships were evaluated using:

- Principal component analysis (PCA)

Differentially expressed genes (DEGs) were identified using:

- **DESeq2**

Criteria:

- False discovery rate (FDR) < 0.05
- Absolute fold change ≥ 2

Functional enrichment analyses were performed using:

- Gene Ontology (GO)
- Kyoto Encyclopedia of Genes and Genomes (KEGG)
- Reactome pathway analysis

---

## 2. CUT&Tag-seq Analysis

### Experimental design

BMAL1 genomic occupancy was profiled using CUT&Tag sequencing in:

- 3×Flag-BMAL1 ESCs
- Differentiated NSCs

Libraries were sequenced using Illumina NovaSeq with 150-bp paired-end reads.

### Bioinformatic analysis

Raw reads were processed using:

- **bcl2fastq**
- **Cutadapt**

Reads were aligned to GRCh38 using:

- **Bowtie2 v2.5.3**

Downstream analyses included:

### TSS enrichment analysis

- **deepTools computeMatrix**

### Peak calling

- **MACS3 v3.0.1**

### Peak annotation

- **ChIPseeker v1.40.0**

### Motif analysis

- **MEME v5.5.5**

---

## 3. m6A-seq (MeRIP-seq) Analysis

### Experimental design

Low-input m6A sequencing was performed using immunoprecipitation with an anti-m6A antibody.

Both:

- Input RNA libraries
- m6A-immunoprecipitated libraries

were prepared and sequenced using Illumina HiSeq X Ten with 150-bp paired-end reads.

### Bioinformatic analysis

Raw sequencing reads were processed using:

- **MultiQC v1.21**
- **Trim Galore v0.6.10**

Ribosomal RNA removal and alignment were performed using:

- **Bowtie2**
- **HISAT2 v2.2.1**

m6A peak identification was performed using:

- **exomePeak2 v1.16.2**

Peak annotation was performed using:

- **ChIPseeker**

Differential m6A modification analysis and integrated analyses were performed together with transcriptomic datasets.

---

# Software Requirements

The following software and versions were used for analysis:

| Software | Version |
|---|---|
| FastQC | - |
| Trimmomatic | - |
| HISAT2 | 2.2.1 |
| StringTie | - |
| DESeq2 | - |
| Bowtie2 | 2.5.3 |
| Cutadapt | - |
| MACS3 | 3.0.1 |
| ChIPseeker | 1.40.0 |
| exomePeak2 | 1.16.2 |
| MultiQC | 1.21 |
| Trim Galore | 0.6.10 |
| MEME | 5.5.5 |

---

# Repository Structure

This repository is licensed under the MIT License. See the LICENSE file for details.
