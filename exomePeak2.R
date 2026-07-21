rm(list=ls())
options(stringsAsFactors = F)
library(getopt)
spec <- matrix(c(    
  'help',   'h',  0, 'logical', 
  'gtf',    'g',  1,  'character', 
  'bamdir', 'b',  1, 'character',
  'ip' ,    'i',  1, 'character',
  'input' , 'I',  1, 'character',
  'lib',    'l',  2, 'character',
  'mode',   'm',  2, 'character',
  'pvalue', 'P',  2, 'numeric',
  'fc',     'f',  2, 'numeric',
  'od',     'o',  1, 'character'), 
  byrow = TRUE, ncol = 4)
opt <- getopt(spec)

print_usage <- function(spec=NULL){
  cat(' Rscript exomePeak2.R  --gtf Homo_sapiens.gtf --bamdir Hisat2 --ip ip1,ip2 --input input1,input2 --od ./
Options: 
    --help   -h   NULL
    --gtf   character    gtf file, genome annotation [forced]
    --bamdir character   mapping file, bam dir [forced]
    --ip     character    IP Sample, used for peak calling only [optional]
    --input  character    INPUT Sample, used for peak calling only [optional]
                          
    --lib    character    the protocal type of the RNA-seq library, can be 
                          one in c("unstranded", "1st_strand", "2nd_strand"); 
                          default = "unstrand"
                          
    --mode   character    a character specifies the scope of peak calling on genome, can be 
                          one of c("exon", "full_transcript", "whole_genome"); 
                          Default = "whole_genome"
    
    --pvalue numeric      the cutoff on p values in peak calling, default = 1e-05
    --fc     numeric      the cutoff on IP over input fold changes in peak calling, default = 2
    --od     character    outdir [forced]
      \n')
  q(status=1)
}

if ( !is.null(opt$help) | is.null(opt$gtf) | is.null(opt$bamdir) ) { print_usage(spec) }
if ( is.null(opt$ip) |is.null(opt$input)) { print_usage(spec) }
# if ( is.null(opt$paired) )  { opt$paired <- TRUE }
if ( is.null(opt$lib) )     { opt$lib <- "unstrand" }
if ( is.null(opt$mode) )    { opt$mode <- "full_transcript" }
if ( is.null(opt$pvalue) )  { opt$pvalue <-  1e-05 }
if ( is.null(opt$fc) )      { opt$fc <- 2 }
if (!file.exists(opt$od))   { dir.create(opt$od) }

library(exomePeak2)

## peak calling for every sample or groups
ip <- unlist(strsplit(opt$ip,split=','))
input <- unlist(strsplit(opt$input,split=','))

ip_bam <- paste0(opt$bamdir,ip,'.Hisat_aln.sorted.bam')
input_bam <- paste0(opt$bamdir,input,'.Hisat_aln.sorted.bam')

print(paste('ip_bam:',ip_bam))
print(paste('input_bam',input_bam))

library(GenomicFeatures)
# txdb <- GenomicFeatures::makeTxDbFromGFF(file = opt$gtf, format="gtf", dataSource="Ensembl")
txdb <- makeTxDbFromGFF(file = "~/refseq/genome/human/GRCh38.111/Homo_sapiens.GRCh38.111.gtf.gz", format="gtf", 
                        dataSource="Ensembl", organism="Homo sapiens")


result <- exomePeak2(bam_ip = ip_bam, 
                     bam_input = input_bam,
                     
                     # data type
                     # paired_end = opt$paired,
                     strandness  = opt$lib,
                     # txdb = txdb, 
                     gff=opt$gtf,
                     # genome = "hg38", no
                     # genome="GRCh38",
                     # cut off
                     p_cutoff = opt$pvalue,
                     # log2FC_cutoff = log2(as.numeric(opt$fc)),
                     mode = opt$mode, 
                     
                     # save
                     save_dir = opt$od,
)

## result write out default