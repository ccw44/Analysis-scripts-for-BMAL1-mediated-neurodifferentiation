library(biomaRt)

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

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

get_transcript_sequence <- function(transcript_id) {
  seq<-getSequence(id = transcript_id, type = "ensembl_transcript_id", seqType = "transcript_exon_intron",  mart = ensembl)
  seq_up<- getSequence(id = transcript_id, type = "ensembl_transcript_id", seqType = "transcript_flank", upstream  = 400, mart = ensembl)
  seq_down<- getSequence(id = transcript_id, type = "ensembl_transcript_id", seqType = "transcript_flank", downstream   = 400, mart = ensembl)
  seq_start<-paste(seq_up,substr(seq, 1, 400),sep='')
  seq_end<-paste(substr(seq, nchar(seq) - 399, nchar(seq)),seq_down,sep='')
  return(list(start=seq_start,end=seq_end))
}

for (csv_file in args) {
  if (file.exists(csv_file)) {
    gene_list <- read.csv(csv_file, header = TRUE)
    transcript_ids <- gene_list$ENSEMBL  
    sequences <- list()
    for (transcript_id in transcript_ids) {
      repeat {
        tryCatch({
          seq_data <- get_transcript_sequence(transcript_id)
          sequences[[transcript_id]] <- seq_data
          break
        }, error = function(e) {
          cat("Error occurred for transcript_id:", transcript_id, "\n")
          cat("Error message:", conditionMessage(e), "\n")
          cat("Retrying...\n")
          Sys.sleep(5)
        })
      }
    }
    
    result <- data.frame(transcript_id = transcript_ids,
                         seq_start = sapply(sequences, function(x) x$start[1]),                     
                         seq_end = sapply(sequences, function(x) x$end[1]))
    
    output_file <- file.path(output_folder, paste0(sub(".csv", "", basename(csv_file)), "_result.csv"))
    write.csv(result, file = output_file, row.names = FALSE)
    
    cat("Processed", length(transcript_ids), "gene IDs from", csv_file, "\n")
  } else {
    cat("Error: Input file does not exist -", csv_file, "\n")
  }
}