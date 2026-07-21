run_CBN <- function(enrichlist,check,cluster=NULL,
                            path,pro,
                            dat, grouplist,control,case,
                            enrichtype
){
  library(BiocFileCache)
  library(CBNplot)
  library(org.Hs.eg.db)
  library(DESeq2)
  library(enrichplot)
  library(parallel)
  colData <- data.frame(row.names=grouplist$sample, 
                        grouplist1=grouplist$group,
                        batch=grouplist$batch
  )
  colData<-subset(colData,grouplist1 %in% c(control,case))
  grouplist<-subset(grouplist,group %in% c(control,case))
  grouplist1 <- grouplist$group

  dat<-dat[match(rownames(colData), names(dat))]
  dds<- DESeqDataSetFromMatrix(countData = dat,
                               colData = colData,
                               design = ~ grouplist1+batch)
  
  filt <- rowSums(counts(dds) < 10) > dim(colData)[1]*0.9 
  dds <- dds[!filt,]
  dds$grouplist1<-relevel(factor(grouplist1),ref=control)
  dds2 <- DESeq(dds)
  
  v <- vst(dds2, blind = FALSE)
  vsted <- assay(v)
  
  # pway<-DOSE::setReadable(enrichlist, OrgDb='org.Hs.eg.db',
  #                         keyType='ENTREZID')

  pway<-enrichlist
  pway <- enrichplot::pairwise_termsim(pway, method = "JC")
  pway<-pway %>% arrange(desc(EnrichmentFactor),pvalue)
  p1<-pway@result
  incSample <- rownames(subset(colData, grouplist1 == control))
  
  enrichlist2<-enrichlist@result
  enrichlist2<-subset(enrichlist2,pvalue<0.05)
  pathID<-enrichlist2[grepl(check,enrichlist2$Description),'ID']
  

  cache_dir <- if (interactive()) {
  } else {
    message("Non-interactive mode detected, using a temporary cache directory.")
    tempdir()
  }
  if (!interactive() && cache_dir == tempdir()) {
    message("Cache usage skipped in non-interactive mode.")
    bfc <- NULL
  } else {
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
      message("Cache directory created at: ", cache_dir)
    }
    bfc <- BiocFileCache(cache_dir)
  }
  
  
  if (grepl('rea',enrichtype)){
    tryCatch({
    bg<-bngeneplot(pway,
               exp = vsted,
               pathNum = which(pway@result$ID %in% pathID),
               labelSize = 3, 
               shadowText = TRUE,
               strThresh = 0.5, 
               R = 10, 
               showDir = T, 
               pathDb = "reactome", 
               expRow = "ENSEMBL",
               convertSymbol = T,
               strengthPlot = T, 
               nStrength = 10, 
    )+  theme(text = element_text(family = "sans"))
    ggsave(bg,filename = file.path(path,pro,paste0(enrichtype,"_",check,'_cbngeneplot.png')),height=10,width=10)
    }, error = function(e) {
      cat("Cannot paint bngeneplot: ", e$message, "\n")
    })
    tryCatch({
    bp<-bnpathplot(results = pway, 
               exp = vsted, 
               # expSample = incSample,
               R = 10,
               # nCategory = 50,
               labelSize=5, 
               shadowText=T,
               compareRef=T,
               expRow = "ENSEMBL",
               showDir=T,	
               qvalueCutOff=0.05)+  
      theme(text = element_text(family = "sans")) 
    ggsave(bp,filename = file.path(path,pro,paste0(enrichtype,"_",check,'_cbnpathplot.png')),scale=1,height=12,width=12)
    ggsave(bp,filename = file.path(path,pro,paste0(enrichtype,"_",check,'_cbnpathplot.pdf')),height=12,width=12,units = "in", dpi = 300)
    }, error = function(e) {
      cat("Cannot paint bnpathplot: ", e$message, "\n")
    })
    
  } else if(grepl('go',enrichtype)) {
    tryCatch({
    bg<-bngeneplot(results = pway,
             exp = vsted,
             pathNum =which(pway@result$ID %in% pathID), 
             orgDb = org.Hs.eg.db,
             R = 10,
             showDir = T,
             convertSymbol = T,
             expRow = "ENSEMBL",
             )+  
      theme(text = element_text(family = "sans"))
    ggsave(bg,filename = file.path(path,pro,paste0(enrichtype,"_",check,'_cbngeneplot.png')),scale=1,height=12,width=12)
    }, error = function(e) {
      cat("Cannot paint bngeneplot: ", e$message, "\n")
      cat("Change R from 10 to 5 and try")
      bg<-bngeneplot(results = pway,
                     exp = vsted,
                     pathNum =which(pway@result$ID %in% pathID), 
                     orgDb = org.Hs.eg.db,
                     R = 5,
                     showDir = T,
                     convertSymbol = T,
                     expRow = "ENSEMBL",
      )+  theme(text = element_text(family = "sans"))
      ggsave(bg,filename = file.path(path,pro,paste0(enrichtype,"_",check,'_cbngeneplot.png')),scale=1,height=10,width=10)
    })
  }
  
  print('CBNplot is finished!')
}
