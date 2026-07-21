library(Gviz)
library(GenomeInfoDb)


library(rtracklayer)

gtf_file <- "Homo_sapiens.GRCh38.111.gtf.gz"

peak_na_file<-readPeakFile("peaks2.xls")

#
axis_track <- GenomeAxisTrack(name="hg38")
i_track <- IdeogramTrack(genome = "hg38", chromosome = 'chr1')

txdb <- makeTxDbFromGFF(file = gtf_file, 
                        format="gtf", 
                        dataSource="Ensembl", organism="Homo sapiens")
seqlevels(txdb) <- paste0("chr", seqlevels(txdb))
k <- keys(txdb, keytype = "GENEID")
gene_symbols <- mapIds(org.Hs.eg.db, keys = k, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
options(ucscChromosomeNames = T)
gene_track<- GeneRegionTrack(txdb, genome = "hg38", name = "YTHDF2",
                             symbol = gene_symbols,
                             showId = TRUE,
                             geneSymbol = TRUE)


peak_track0 <- DataTrack(range = "B-1_L1_bowtie2_chr.bw", genome = "hg38", type = "histogram",
                        name = "H1 BMAL1 Peaks", col.histogram = "#2fe1e9", fill.histogram = "#2fe1e9",ylim = c(0, 25))
peak_track <- DataTrack(range = "B-3_L1_bowtie2_chr.bw", genome = "hg38", type = "histogram",
                        name = "NPC BMAL1 Peaks", col.histogram = "red", fill.histogram = "red",ylim = c(0, 25))
peak_track2 <- DataTrack(range = "G_L1_bowtie2_chr.bw", genome = "hg38", type = "histogram",
                        name = "NPC IgG Peaks", col.histogram = "grey", fill.histogram = "grey",ylim = c(0, 25))

bed_data <- read.table("narrow_summits.bed", header = FALSE, sep = "\t")
bed_data$V1 <- paste0("chr", bed_data$V1)
write.table(bed_data, "summits_chr.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
bed_data<-"narrow_summits_chr.bed"
a_track <- AnnotationTrack(bed_data, genome = "hg38", name = "Region", col = "purple", fill = "purple")

library(Biostrings)
# YTHDF2 peak region
seq<-'CCCTAGGGGTCGGGTTACCTCCCAGATGGGGGTGGGGTATCACTTGCAGCTGGAAGGCAGGTTCATTTTAAGGTCTCGCGGGCTTCGGAGGCCCCTCCCTCACTGAGAAACCGGCTTTTGTCGTTCTTCTGAAGCCAGAGCTAGTCTTTCCAGGTGTTAGTCGAAACCTCGTGGTGCGACCCTGGTCGTCCCAAACCCCCTAGGCCTTAATCC'
dna <- DNAStringSet(seq)
names(dna) <- "BMAL1-binded seq"
options(ucscChromosomeNames = FALSE)
seqTrack <- SequenceTrack(dna, name = "Sequence",genome="hg38")

png(filename=paste0("Gviz_peaks_","YTHDF2",".png") ,width=2000, height=1500, res=100)
pdf(file=paste0("Gviz_peaks_","YTHDF2",".pdf"), width = 20, height = 15)


plotTracks(list(
                i_track,
                axis_track,
                gene_track,
                peak_track0,
                peak_track,
                peak_track2,
                a_track
                ),
           from = 28736000, to = 28741000,
           # from =28734624, to = 28771775,
           # from =45845118, to = 45885244,
           chromosome = "chr1",
           showTitle = TRUE,
           col.title = "black",
           col.axis = "black",
           cex.title = 1.2,    # 调整标题字体大小
           cex.axis = 1   # 调整坐标轴字体大小
           )
dev.off()