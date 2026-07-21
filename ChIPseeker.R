library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicFeatures)

# annotation file gtf
txdb <- makeTxDbFromGFF(file = "GRCh38.111/Homo_sapiens.GRCh38.111.gtf.gz", format="gtf", 
                        dataSource="Ensembl", organism="Homo sapiens")
# annotation file gff3
txdb2 <- makeTxDbFromGFF(file = "GRCh38.111/Homo_sapiens.GRCh38.105.gff3.gz", format="gff", 
                        dataSource="Ensembl", organism="Homo sapiens")


# CUTtag
peak_file<-readPeakFile("peaks2.xls")

## peak annotation
# region select："Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"

## CUTtag
# broad
# gtf
peak_anno<- annotatePeak(peak_file,
                           tssRegion = c(-3000, 3000),
                           TxDb = txdb,
                           assignGenomicAnnotation = TRUE,
                           # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Intergenic"),
                           genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Promoter"),
                           # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron"),
                           addFlankGeneInfo = TRUE,
                           flankDistance = 5000,
                         annoDb="org.Hs.eg.db")
# broad
# gff3
peak_anno2<- annotatePeak(peak_file,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb2,
                         assignGenomicAnnotation = TRUE,
                         # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Intergenic"),
                         genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Promoter"),
                         # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron"),
                         addFlankGeneInfo = TRUE,
                         flankDistance = 5000,
                         annoDb="org.Hs.eg.db"
                         )
# narrow
# gtf
peak_na_anno<- annotatePeak(peak_na_file,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         assignGenomicAnnotation = TRUE,
                         # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Intergenic"),
                         genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Promoter"),
                         # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron"),
                         addFlankGeneInfo = TRUE,
                         flankDistance = 5000,
                         annoDb="org.Hs.eg.db"
                         )
# narrow
# gff3
peak_na_anno2<- annotatePeak(peak_na_file,
                            tssRegion = c(-3000, 3000),
                            TxDb = txdb2,
                            assignGenomicAnnotation = TRUE,
                            # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Intergenic"),
                            genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Promoter"),
                            # genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron"),
                            addFlankGeneInfo = TRUE,
                            flankDistance = 5000,
                            annoDb="org.Hs.eg.db")



# CUTtag
pdf(file = "Anno.Pie.pdf",width =12,height = 8)
plotAnnoPie(peak_anno)
dev.off()