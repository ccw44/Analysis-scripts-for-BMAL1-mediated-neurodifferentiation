library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(reshape2)
library(ChIPseeker)
library(ggpubr)

m6acontent<-function(path1, type,diff=F){

  
  seqlevelsStyle(peak_file) <- "UCSC" 
  peak_file <- keepStandardChromosomes(peak_file, pruning.mode = "coarse")
  
  peak_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peak_file)
  
  # Count adenosines in each peak sequence
  a_counts <- vapply(peak_sequences, function(seq) {
    alphabetFrequency(seq)["A"]
  }, numeric(1))
  # Calculate the Total Adenosine Counts in Peaks
  a_peaks <- sum(a_counts)
  
  # Calculate the Total Peak Width #1
  w_peak <- width(peak_file)
  ap_peaks <- a_peaks / sum(w_peak)
  
  # Calculate the Total Peak Width #2
  # tw_peak <- sum(width(peak_file))
  # ap_peaks <- a_peaks / tw_peak
  return(ap_peaks)
}

m6acontent2 <- function(path1,path2, type, 
                       bam = NULL, 
                       mapping_stats = NULL) {
  library(GenomicAlignments)
  library(simsem) 
  library(Biostrings)
  
  peak_file <- readPeakFile(file.path(path1,type,"exomePeak2_output/peaks.bed"))
  seqlevelsStyle(peak_file) <- "UCSC" 
  peak_file <- keepStandardChromosomes(peak_file, pruning.mode = "coarse")
  
  if(!is.null(bam)) {
    ip_bam<-file.path(path2,paste0(bam,'IP.Hisat_aln.sorted.bam'))
    ip <- readGAlignments(ip_bam)
    input_bam<-file.path(path2,paste0(bam,'Input.Hisat_aln.sorted.bam'))
    input <- readGAlignments(input_bam)
    
    ip_coverage <- getCoverage(peak_file, ip)
    input_coverage <- getCoverage(peak_file, input)
    
    enrichment_scores <- ip_coverage / input_coverage
    mean_enrichment <- mean(enrichment_scores, na.rm = TRUE)
    
    ip_quality <- mean(mcols(ip)$mapq, na.rm = TRUE)
    input_quality <- mean(mcols(input)$mapq, na.rm = TRUE)
  }
  
  peak_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peak_file)
  
  a_counts <- vapply(peak_sequences, function(seq) {
    alphabetFrequency(seq)["A"]
  }, numeric(1))
  
  a_peaks <- sum(a_counts)
  total_width <- sum(width(peak_file))
  raw_ratio <- a_peaks / total_width
  
  adjusted_ratio <- raw_ratio
  
  if(!is.null(mapping_stats)) {
    ip_mapping_rate <- mapping_stats$ip_mapping_rate
    input_mapping_rate <- mapping_stats$input_mapping_rate
    
    mapping_ratio <- ip_mapping_rate / input_mapping_rate
    adjusted_ratio <- adjusted_ratio * mapping_ratio
  }
  
  if(!is.null(ip_bam) && !is.null(input_bam)) {
    depth_ratio <- mean(ip_coverage) / mean(input_coverage)
    
    enrichment_weight <- mean_enrichment / (1 + mean_enrichment)  
    
    quality_weight <- ((ip_quality + input_quality) / 2) / 60
    
    adjusted_ratio <- adjusted_ratio * depth_ratio * enrichment_weight * quality_weight
  }
  
  return(list(
    raw_m6a_ratio = raw_ratio,
    adjusted_m6a_ratio = adjusted_ratio,
    total_peaks = length(peak_file),
    total_adenines = a_peaks,
    total_width = total_width,
    sequencing_metrics = if(!is.null(ip_bam) && !is.null(input_bam)) {
      list(
        ip_depth = mean(ip_coverage),
        input_depth = mean(input_coverage),
        depth_ratio = depth_ratio,
        enrichment_score = mean_enrichment,
        ip_quality = ip_quality,
        input_quality = input_quality,
        enrichment_weight = enrichment_weight,
        quality_weight = quality_weight
      )
    } else NULL,
    mapping_stats = if(!is.null(mapping_stats)) {
      list(
        ip_mapping_rate = mapping_stats$ip_mapping_rate,
        input_mapping_rate = mapping_stats$input_mapping_rate,
        mapping_ratio = mapping_ratio
      )
    } else NULL
  ))
}



run_m6acount<-function(
    ap1,ap2,
    con,case,
    path,pro){
  df_ap <- data.frame(
    group = rep(c(case, con), each = length(ap1)),  
    ap = c(ap1, ap2) 
  )
  mdf_ap<- reshape2::melt(df_ap)
  mdf_ap$group <- factor(mdf_ap$group, levels =c(con,case))
  pbar <- ggplot(mdf_ap,aes(group,100*value))+
    stat_summary(mapping=aes(fill = group),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.5)+
    scale_fill_manual(values = c("#ff7f0e","#1f77b4"))+
    # stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",linewidth=0.2)+
    stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
                 width = 0.3, 
                 linewidth=0.5, 
                 position = position_dodge(1)
    )+
    stat_compare_means(method = "t.test",
                       label = "p.signif",
                       comparisons = list(c(con,case)),
                       size = 10
    )+
    labs(x = "",y = "m6A Content(%)")+
    scale_y_continuous(expand = c(0,0),limit=c(0,50),
                       breaks = seq(0,50,by = 10))+
    theme_classic()+
    theme(panel.background=element_rect(fill="white",colour="black",size=0.25), 
          axis.line=element_line(colour="black",size=0.25), 
          axis.title=element_text(size=16,color="black"), 
          axis.text = element_text(size=16,color="black"), 
          axis.ticks.x = element_blank(),
          legend.position="none"
    )
  ggsave(pbar, file=file.path(path,pro,'m6acontent.png'),scale=1,height=6,width=6)
  ggsave(pbar, file=file.path(path,pro,'m6acontent.pdf'),height=6,width=6,units = "in", dpi = 300)
}