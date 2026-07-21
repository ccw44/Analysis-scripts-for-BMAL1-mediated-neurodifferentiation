library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicFeatures)

# annotation file gtf
txdb <- makeTxDbFromGFF(file = "Homo_sapiens.GRCh38.111.gtf.gz", format="gtf", 
                        dataSource="Ensembl", organism="Homo sapiens")


txdb <- makeTxDbFromGFF(file = "hg38.knownGene.gtf.gz", 
                        format="gtf", 
                        dataSource="Knowngene", 
                        organism="Homo sapiens")
txdb <- makeTxDbFromGFF(file = "hg38.ncbiRefSeq.gtf.gz", 
                        format="gtf", 
                        dataSource="Knowngene", 
                        organism="Homo sapiens")



peak_file <- readPeakFile("peaks.bed")

library(GenomicRanges)
peak_file@seqinfo@seqnames
seqlevels(peak_file) <- sub("^chr", "", seqlevels(peak_file))


## peak annotation
#region select："Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"
peak_anno<- annotatePeak(peak_file,
                           tssRegion = c(-3000, 3000),
                           TxDb = txdb,
                           assignGenomicAnnotation = TRUE,
                           # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Intergenic"),
                           # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Promoter"),
                           genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron"),
                           addFlankGeneInfo = TRUE,
                           flankDistance = 5000)
peak_anno<- annotatePeak(peak_file,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         assignGenomicAnnotation = TRUE,
                         # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Intergenic"),
                         # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Promoter"),
                         # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron"),
                         addFlankGeneInfo = TRUE,
                         flankDistance = 5000)

