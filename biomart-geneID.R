#!/usr/bin/Rscript

library(biomaRt)

rerun_script <- function() {
  # 运行脚本
  system("Rscript your_script.R your_csv_file.csv -o /path/to/output/folder > output.log")
}

repeat {
  tryCatch({
    # 运行脚本
    rerun_script()
    # 如果成功运行脚本，则跳出重试循环
    break
  }, error = function(e) {
    cat("Error occurred while running the script:\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("Retrying in 5 seconds...\n")
    Sys.sleep(5)  
  })
}

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

args <- commandArgs(trailingOnly = TRUE)

if ("-o" %in% args) {
  output_index <- which(args == "-o")
  
  if (output_index < length(args)) {
    output_folder <- args[output_index + 1]
  } else {
    cat("Error: No output folder specified.\n")
    quit(status = 1)
  }
}

get_gene_sequence <- function(gene_id) {
  seq<-getSequence(id = gene_id, type = "ensembl_gene_id", seqType = "gene_exon_intron",  mart = ensembl)
  seq_up<- getSequence(id = gene_id, type = "ensembl_gene_id", seqType = "gene_flank", upstream  = 400, mart = ensembl)
  seq_down<- getSequence(id = gene_id, type = "ensembl_gene_id", seqType = "gene_flank", downstream   = 400, mart = ensembl)
  seq_start<-paste(seq_up,substr(seq, 1, 400),sep='')
  seq_end<-paste(substr(seq, nchar(seq) - 399, nchar(seq)),seq_down,sep='')
  return(list(start=seq_start,end=seq_end))
}

for (csv_file in args) {
  gene_list <- read.csv(csv_file, header = TRUE)
  gene_ids <- gene_list$ENSEMBL  
  
  sequences <- list()
  
  for (gene_id in gene_ids) {
    repeat {
      tryCatch({
        seq_data <- get_gene_sequence(gene_id)
        sequences[[gene_id]] <- seq_data
        break
      }, error = function(e) {
        cat("Error occurred for gene_id:", gene_id, "\n")
        cat("Error message:", conditionMessage(e), "\n")
        cat("Retrying...\n")
        Sys.sleep(5)  
      })
    }
  }
  
  result <- data.frame(gene_id = gene_ids,
                       seq_start= sapply(sequences, function(x) x$start[1]),                     
                       seq_end = sapply(sequences, function(x) x$end[1])
  )
  
  output_file <- paste(output_folder, sub(".csv", "_result.csv", basename(csv_file)), sep = "/")
  write.csv(result, file = output_file, row.names = FALSE)
  
  cat("Processed", length(gene_ids), "gene IDs from", csv_file, "\n")
}