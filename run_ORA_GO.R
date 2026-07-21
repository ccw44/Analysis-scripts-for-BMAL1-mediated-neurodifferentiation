run_ORA_GO <- function(deg, path,pro,
                       type=NULL,ct=F,m6a1=F){
  # type: used for merely up or dw genes enrichment for m6a/CT, or up+dw= diff2 to make loliplot2 for CT&RNA data 
  # ct：whether use CUTtag data
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  library(myenrichplot)
  library(cols4all)
  
  # deg<-df_sig
  if (!is.null(type) & isTRUE(ct)){
    deg<-deg %>% filter(!is.na(CT)) 
  }

  df <- bitr(unique(deg$name), fromType = "ENSEMBL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  DEG<-deg
  DEG<-merge(DEG,df,by.y='ENSEMBL',by.x='name',all.x=T)

  ## |logFC|>1
  if (is.null(type)|(!is.null(type) && grepl('up',type))){gene_up= unique(DEG[DEG$g == 'up','ENTREZID'])
                                                          gene_up2= unique(DEG[DEG$g == 'up','SYMBOL_1'])   # m6a
                                                          }
  if (is.null(type) |(!is.null(type) && grepl('dw',type))){gene_down=unique(DEG[DEG$g == 'down','ENTREZID'])
                                                            gene_dw2= unique(DEG[DEG$g == 'down','SYMBOL_1']) # m6a
                                                            }
  if ((!is.null(type) && grepl('diff2',type))){
                gene_up3= unique(DEG[DEG$g == 'up','ENTREZID'])
                gene_dw3=unique(DEG[DEG$g == 'down','ENTREZID'])
  }
  if (is.null(type)&& isFALSE(m6a1) ){
                    gene_diff=unique(c(gene_up,gene_down))
                    uptop<-DEG %>% filter(g=='up') %>% arrange(desc(logFC)) %>% head(100) 
                    uptop<-uptop$ENTREZID
                    dwtop<-DEG %>% filter(g=='down') %>% arrange(logFC) %>% head(100) 
                    dwtop<-dwtop$ENTREZID
  }
  
  
  goplot2 <- function(goresult,pattern){
    goresult<-goresult@result %>% arrange(pvalue)
    if (nrow(goresult)>=30){
      goresult<- goresult %>% head(30) 
    }
    goresult$Description <- factor(goresult$Description,
                                    levels = rev(goresult$Description))
    mytheme <- theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      axis.text.y = element_blank(),
      plot.title = element_text(size = 14,
                                hjust = 0.5,
                                face = "bold"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      plot.margin = margin(t = 5.5,
                           r = 10,
                           l = 5.5,
                           b = 5.5),
      axis.ticks.y = element_blank()
    )
    
    if ( 'ONTOLOGY' %in%  names(goresult)) {
    ggplot(data = goresult, aes(x = -log10(pvalue), y = Description,fill = RichFactor)) +
      scale_fill_distiller(palette = pattern,direction = 1)+ 
      facet_grid(ONTOLOGY~., scale="free") +
      geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
      scale_x_continuous(expand=c(0,0))+
      labs(x = "-log10(pvalue)", y = "", title = "GO pathways") +
      geom_text(size=4,aes(x = 0.1,
                    label = Description),
                hjust = 0)+
      # geom_text(size=3,aes(x=0.1,label=geneID))+
      theme_classic() + 
      mytheme
    } else {
      ggplot(data = goresult, aes(x = -log10(pvalue), y = Description,fill = RichFactor)) +
        scale_fill_distiller(palette = pattern,direction = 1)+ 
        geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
        scale_x_continuous(expand=c(0,0))+
        labs(x = "-log10(pvalue)", y = "", title = "GO BP pathways") +
        geom_text(size=4,aes(x = 0.1, 
                      label = Description),
                  hjust = 0)+
        # geom_text(size=3,aes(x=0.1,label=geneID))+
        theme_classic() + 
        mytheme
    }
  }

  
  ifelse(!dir.exists(file.path(path,pro)),
         dir.create(file.path(path,pro)),
         "Directory Exists")
  
  if (exists('gene_up')){
  ## GO-all_up
  go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all",
                 pAdjustMethod = "BH",readable= TRUE,
                 pvalueCutoff = 0.3,
                 qvalueCutoff =0.5) 
  go <- clusterProfiler::simplify(go, cutoff=0.7, by="p.adjust", select_fun=min)

  go<-DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(go@result,file = file.path(path,pro,'gene_up_GO_all.csv'))
  up_go<-go
  go@result <- go@result %>% filter(p.adjust<0.05 & qvalue<0.2)
  
  tryCatch({

  bar<-goplot2(go,pattern='YlOrRd')
  ggsave(bar,filename = file.path(path,pro,'gene_up_GO_all_barplot.png'),scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint barplot-goall-up: ", e$message, "\n")
  })
  
  tryCatch({
  tp<-myenrichplot::treeplot(pairwise_termsim(go),
               cluster.params=list(method='average',n=3),
               nCluster = 3)
  ggsave(filename=file.path(path,pro,"gene_up_GO_all_treeplot.png"),plot=tp,scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint treeplot-goall-up: ", e$message, "\n")
  })
  
  ## GO-BP_up
  go_BP<- enrichGO(gene_up, keyType = "ENTREZID",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    readable      = TRUE,
    pvalueCutoff =0.5,
    qvalueCutoff =0.5)
  go_BP <- clusterProfiler::simplify(go_BP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  go_BP<-DOSE::setReadable(go_BP, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(go_BP@result,file = file.path(path,pro,'gene_up_go_BP.csv'))
  up_gobp<-go_BP
  go_BP@result <- go_BP@result %>% filter(p.adjust<0.05 & qvalue<0.2)
  tryCatch({

  barbp<-goplot2(go_BP,pattern='YlOrRd')
  ggsave(barbp,filename = file.path(path,pro,'gene_up_go_BP_barplot.png'),scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint barplot-gobp-up: ", e$message, "\n")
  })
  
  tryCatch({
  tp<-myenrichplot::treeplot(pairwise_termsim(go_BP),
               cluster.params=list(method='average',n=3),
               nCluster = 3)
  ggsave(filename=file.path(path,pro,"gene_up_GO_BP_treeplot.png"),plot=tp,scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint treeplot-gobp-up: ", e$message, "\n")
  })
  
  }
  
  if (exists('gene_down')){
  ## GO-all_down
  go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all",pAdjustMethod = "BH",readable= TRUE,
                 pvalueCutoff= 0.5,
                 qvalueCutoff =0.5) 
  go <- clusterProfiler::simplify(go, cutoff=0.7, by="p.adjust", select_fun=min)
  
  go<-DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(go@result,file = file.path(path,pro,'gene_down_GO_all.csv'))
  down_go<-go
  go@result <- go@result %>% filter(p.adjust<0.05 & qvalue<0.2)
  
  tryCatch({

  bar<-goplot2(go,pattern='PuBu')
  ggsave(bar,filename = file.path(path,pro,'gene_down_GO_all_barplot.png'),scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint barplot-goall-dw: ", e$message, "\n")
  })
  
  tryCatch({
  tp<-myenrichplot::treeplot(pairwise_termsim(go),
               cluster.params=list(method='average',n=3),
               nCluster = 3)
  ggsave(filename=file.path(path,pro,"gene_down_GO_all_treeplot.png"),plot=tp,scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint treeplot-goall-dw: ", e$message, "\n")
  })
  
  ## GO-BP_down
  go_BP<- enrichGO(gene_down, keyType = "ENTREZID",
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   readable      = TRUE,
                   pvalueCutoff = 0.3,
                   qvalueCutoff =0.5)
  go_BP <- clusterProfiler::simplify(go_BP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  go_BP<-DOSE::setReadable(go_BP, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(go_BP@result,file = file.path(path,pro,'gene_down_go_BP.csv'))
  down_gobp<-go_BP
  go_BP@result <- go_BP@result %>% filter(p.adjust<0.05 & qvalue<0.2)
  
  tryCatch({

  barbp<-goplot2(go_BP,pattern='PuBu')
  ggsave(barbp,filename = file.path(path,pro,'gene_down_go_BP_barplot.png'),scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint barplot-gobp-dw: ", e$message, "\n")
  })

  
  tryCatch({
    tp <- myenrichplot::treeplot(pairwise_termsim(go),
                         cluster.params=list(method='average', n=3),
                         nCluster = 3)
    ggsave(filename = file.path(path, pro, "gene_down_GO_all_treeplot.png"), plot = tp, scale = 1, height = 8, width = 14)
  }, error = function(e) {
    cat("Cannot paint treeplot-gobp-dw: ", e$message, "\n")
  })
  
  }
  
  if (exists('gene_diff')){
  ## GO-all_diff
  go <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="all",pAdjustMethod = "BH",readable= TRUE,
                 pvalueCutoff = 0.3,
                 qvalueCutoff =0.5) 
  go <- clusterProfiler::simplify(go, cutoff=0.7, by="p.adjust", select_fun=min)
  
  go<-DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(go@result,file = file.path(path,pro,'gene_diff_GO_all.csv'))
  diff_go<-go
  go@result <- go@result %>% filter(p.adjust<0.05 & qvalue<0.2)
  
  tryCatch({
  bar<-barplot(go, split="ONTOLOGY",font.size =10)+ 
    facet_grid(ONTOLOGY~., scale="free") + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(bar,filename = file.path(path,pro,'gene_diff_GO_all_barplot.png'),scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint barplot-godiff: ", e$message, "\n")
  })
  tryCatch({
  tp<-myenrichplot::treeplot(pairwise_termsim(go),
               cluster.params=list(method='average',n=3),
               nCluster = 3)
  ggsave(filename=file.path(path,pro,"gene_diff_GO_all_treeplot.png"),plot=tp,scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint treeplot-godiff: ", e$message, "\n")
  })
  ## GO-BP_diff
  go_BP<- enrichGO(gene_diff, keyType = "ENTREZID",
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   readable      = TRUE,
                   pvalueCutoff = 0.3,
                   qvalueCutoff =0.5)
  go_BP <- clusterProfiler::simplify(go_BP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  go_BP<-DOSE::setReadable(go_BP, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(go_BP@result,file = file.path(path,pro,'gene_diff_go_BP.csv'))
  diff_gobp<-go_BP
  go_BP@result <- go_BP@result %>% filter(p.adjust<0.05 & qvalue<0.2)
  tryCatch({
  bar<-barplot(go_BP, font.size =10)+ 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(bar,filename = file.path(path,pro,'gene_diff_go_BP_barplot.png'),scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint barplot-gobp-diff: ", e$message, "\n")
  })
  tryCatch({
  tp<-myenrichplot::treeplot(pairwise_termsim(go_BP),
               cluster.params=list(method='average',n=3),
               nCluster = 3)
  ggsave(filename=file.path(path,pro,"gene_diff_GO_BP_treeplot.png"),plot=tp,scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint treeplot-gobp-diff: ", e$message, "\n")
  })
  # GO-BP_uptop
  go_BP<- enrichGO(uptop, keyType = "ENTREZID",
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   readable      = TRUE,
                   pvalueCutoff = 0.3,
                   qvalueCutoff =0.5)
  go_BP <- clusterProfiler::simplify(go_BP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  go_BP<-DOSE::setReadable(go_BP, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(go_BP@result,file = file.path(path,pro,'uptop_go_BP.csv'))
  uptop_gobp<-go_BP
  go_BP@result <- go_BP@result %>% filter(p.adjust<0.05 & qvalue<0.2)
  
  tryCatch({
  # bar<-barplot(go_BP, font.size =10)+ 
  #   scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  barbp<-goplot2(go_BP,pattern='YlOrRd')
  ggsave(barbp,filename = file.path(path,pro,'uptop_go_BP_barplot.png'),scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint barplot: ", e$message, "\n")
  })
  
  tryCatch({
  tp<-myenrichplot::treeplot(pairwise_termsim(go_BP),
               cluster.params=list(method='average',n=2),
               nCluster = 2)
  ggsave(filename=file.path(path,pro,"uptop_GO_BP_treeplot.png"),plot=tp,scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint treeplot-gobp-uptop: ", e$message, "\n")
  })
  
  # GO-BP_dwtop
  go_BP<- enrichGO(dwtop, keyType = "ENTREZID",
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   readable      = TRUE,
                   pvalueCutoff = 0.3,
                   qvalueCutoff =0.5)
  go_BP <- clusterProfiler::simplify(go_BP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  go_BP<-DOSE::setReadable(go_BP, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(go_BP@result,file = file.path(path,pro,'dwtop_go_BP.csv'))
  dwtop_gobp<-go_BP
  go_BP@result <- go_BP@result %>% filter(p.adjust<0.05 & qvalue<0.2)
  
  tryCatch({
    barbp<-goplot2(go_BP,pattern='PuBu')
  ggsave(barbp,filename = file.path(path,pro,'dwtop_go_BP_barplot.png'),scale=1,height=8,width=14)
  }, error = function(e) {
    cat("Cannot paint barplot-gobp-dwtop: ", e$message, "\n")
  })
  
  tryCatch({
  tp<-myenrichplot::treeplot(pairwise_termsim(go_BP),
               cluster.params=list(method='average',n=2),
               showCategory=12,
               fontsize=4,
               label_format=12,
               extend=0.2,
               nCluster = 3,
               nWords =6,
               color='p.adjust',
               group_color=brewer.pal(3, "Set2")[1:3]
               )
  ggsave(filename=file.path(path,pro,"dwtop_GO_BP_treeplot.png"),plot=tp,scale=1,height=7,width=12)
  ggsave(filename=file.path(path,pro,"dwtop_GO_BP_treeplot.pdf"),plot=tp,height=7,width=12,units = "in", dpi = 300)
  }, error = function(e) {
    cat("Cannot paint treeplot-gobp-dwtop: ", e$message, "\n")
  })
  
  save(diff_go,down_go,up_go,
       diff_gobp,down_gobp,up_gobp,
       uptop_gobp,dwtop_gobp,
       file =file.path(path,pro,'go_results.Rdata'))
  }
  
  if (!is.null(type) && grepl('up',type)){
    save(up_go,up_gobp, file =file.path(path,pro,paste(type,'go_results.Rdata',sep='-')))}
  
  tryCatch({
  if (!is.null(type) &&grepl('dw',type)){
    save(down_go,down_gobp,file =file.path(path,pro,paste(type,'go_results.Rdata',sep='-')))}
  }, error = function(e) {
    cat("No Go results: ", e$message, "\n")
  })
  
  
  if (exists('gene_up')& exists('gene_down')){
    goup2 <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="BP",
                     pAdjustMethod = "BH",readable= TRUE,
                     pvalueCutoff = 0.3,
                     qvalueCutoff =0.5) 
    goup2 <- clusterProfiler::simplify(goup2, cutoff=0.7, by="p.adjust", select_fun=min)
    goup2<-DOSE::setReadable(goup2, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    goup2 <- goup2@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    if(nrow(goup2)>=1){goup2$group=1}
    
    godw2 <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="BP",
                     pAdjustMethod = "BH",readable= TRUE,
                     pvalueCutoff = 0.3,
                     qvalueCutoff =0.5) 
    godw2 <- clusterProfiler::simplify(godw2, cutoff=0.7, by="p.adjust", select_fun=min)
    godw2<-DOSE::setReadable(godw2, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    godw2<- godw2@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    if(nrow(godw2)>=1){godw2$group=-1}
    
    
    dat=rbind(goup2,godw2)
    dat<-dat %>% subset(qvalue<0.2) %>% arrange(desc(RichFactor))
    # if (nrow(dat)>=20){
    #   dat<- dat %>% head(20) 
    # }
    
    if ((table(dat$group)[1]/table(dat$group)[2]) < (1/20)){
      dat1<-dat %>% filter(group==1) %>% head(19)
      dat2<-dat %>% filter(group==(-1)) %>% head(1)
    } else if ( (table(dat$group)[2]/table(dat$group)[1]) < (1/20) ){
      dat1<-dat %>% filter(group==1) %>% head(1)
      dat2<-dat %>% filter(group==(-1)) %>% head(19)
    }  else{
      dat1<-dat %>% filter(group==1) %>% head(floor(20*(table(dat$group)[2])/length(dat$group)))
      dat2<-dat %>% filter(group==(-1)) %>% head(20-floor(20*(table(dat$group)[2])/length(dat$group)))
    }
    dat<-rbind(dat1,dat2)
    
    dat$logPvalue = -log10(dat$pvalue)
    dat$logPvalue=dat$logPvalue*dat$group 
    dat=dat[order(dat$logPvalue,decreasing = F),]
    
    dat$FoldEnrichment <- dat$FoldEnrichment + runif(nrow(dat), -0.1, 0.1)
    if(nrow(goup2)>=1){datup<-dat %>% filter(group=='1')}
    if(nrow(godw2)>=1){datdw<-dat %>% filter(group=='-1')}
    
    
    dat$FoldEnrichment<-as.character(round(dat$FoldEnrichment,3))
    datup$FoldEnrichment<-as.character(round(datup$FoldEnrichment,3))
    datdw$FoldEnrichment<-as.character(round(datdw$FoldEnrichment,3))
    tryCatch({
    pe <- ggplot(dat, aes(x=`logPvalue`, y=reorder(FoldEnrichment,order(logPvalue, decreasing = F)), fill=`logPvalue`)) +
      geom_col(aes(fill = `logPvalue`), width = 0.1)+
      geom_point(aes(size = RichFactor,
                     color = `logPvalue`))+
      scale_size_continuous(range = c(2, 7)) +
      geom_text(data = datup,
                aes(x = -0.2, y = reorder(FoldEnrichment,order(logPvalue, decreasing = F)), label = reorder(Description,order(logPvalue, decreasing = F))),
                size =5,
                hjust = 1)+ 
      geom_text(data = datdw,
                aes(x = 0.2, y = reorder(FoldEnrichment,order(logPvalue, decreasing = F)), label = reorder(Description,order(logPvalue, decreasing = F))),
                size = 5,
                hjust = 0)+ 
      scale_color_continuous_c4a_div('sunset', mid = 0, reverse = F) +
      scale_fill_continuous_c4a_div('sunset', mid = 0, reverse = F) +
      scale_x_continuous(limits = c(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue))),
                         breaks = seq(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue)),by = 2),
                         labels = abs(seq(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue)),by = 2))
      ) + 
      labs(x = '-log(Pvalue)', y = '')+
      theme_classic()+
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 13), 
        axis.title.x = element_text(size = 13), 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 12), 
        axis.ticks.y=element_blank()
      )+
      guides(
        fill = "none",  
        color = "none"  
      ) 
    ggsave(pe,filename=file.path(path,pro,'gene_up&down_go_BP_loliplot.png'),scale=1,height=12,width=12)
    ggsave(pe,filename=file.path(path,pro,'gene_up&down_go_BP_loliplot.pdf'),height=12,width=12,units = "in", dpi = 300)
    }, error = function(e) {
      cat("Cannot paint loliplot-gobp-updw: ", e$message, "\n")
    })
    
    
    godiff <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="BP",
                      pAdjustMethod = "BH",readable= TRUE,
                      pvalueCutoff = 0.3,
                      qvalueCutoff =0.5) 
    godiff <- clusterProfiler::simplify(godiff, cutoff=0.7, by="p.adjust", select_fun=min)
    godiff<-DOSE::setReadable(godiff, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    got<-godiff
    godiff <- godiff@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    
    godiff <- godiff %>%
      rowwise() %>%
      mutate(
        upn = length(intersect(unlist(strsplit(geneID, "/")), gene_up2)),
        dwn = length(intersect(unlist(strsplit(geneID, "/")), gene_dw2))
      )
    
    
    
    
    
    tryCatch({
      pb <- godiff %>%
        mutate(
          tg = upn + dwn,
          log10_pval = -log10(pvalue)
        ) %>%
        top_n(20, wt = tg) %>% 
        arrange(desc(log10_pval)) %>%
        mutate(Description = fct_reorder(Description, tg)) 
      
     
      tp<-myenrichplot::treeplot(pairwise_termsim(got),
                                 cluster.params=list(method='average',n=3),
                                 nCluster =3)
      clusters <- as.data.frame(tp$data)
      
      pb2 <- pb %>%
        left_join(clusters[, c("label", "group")], by = c("Description" = "label")) %>%
        mutate(Category = group)
      
      
      pb2 <- pb2 %>% filter(!is.na(group))
      pb3 <- pb2 %>%
        pivot_longer(cols = c(upn, dwn), names_to = "regulation", values_to = "count") %>%
        mutate(regulation = ifelse(regulation == "upn", "up", "down")) %>% 
        arrange(group,log10_pval)
      pb3$Description <- fct_inorder(pb3$Description)
      
      pbs<-ggplot(pb3, aes(x = Description, y = count, fill = regulation)) +
        geom_bar(stat = "identity") +  
        scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#009ec5")) + 
        geom_text(data = pb3 %>% filter(regulation == "down"),
                  aes(label = round(log10_pval, 2), y = max(pb3$tg)), 
                  size = 3, 
                  hjust = 1,  
                  vjust = 0.5
        ) +  
        coord_flip() +  
        facet_wrap(~Category, scales = "free_y", ncol = 1) + 
        labs(
          x = "",
          y = "Gene Count",
          fill = pro
        ) +
        theme_classic()+
        theme(
          strip.text = element_blank(),
          # strip.text = element_text(size = 10, face = "bold"), 
          axis.text.y = element_text(size = 8),  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          axis.line.x = element_line(color = "black"), 
          axis.text.x = element_text(size = 10),  
          axis.title.x = element_text(size = 12),  
          axis.ticks.y=element_blank(),
          # strip.placement = "outside",
        ) +
        scale_y_continuous(position = "right",expand = c(0, 0)) 
      
      
      pb3 <- pb2 %>%
        pivot_longer(cols = c(upn, dwn), names_to = "regulation", values_to = "count") %>%
        mutate(regulation = ifelse(regulation == "upn", "up", "down")) %>% 
        arrange(desc(group),log10_pval)
      pb3$Description <- fct_inorder(pb3$Description)
      
      pbs2<-ggplot(pb3, aes(x = Description, y = count, fill = regulation)) +
        geom_bar(stat = "identity") +  
        scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#009ec5")) + 
        geom_text(data = pb3 %>% filter(regulation == "down"),
                  aes(label = round(log10_pval, 2), y = max(pb3$tg)),  
                  size = 7, 
                  hjust = 1,  
                  vjust = 0.5
        ) +  
        scale_x_discrete(labels=function(x) str_wrap(x, width=30)) +
        coord_flip() +  
        labs(
          x = "",
          y = "Gene Count",
          fill = pro
        ) +
        theme_classic()+
        theme(
          strip.text = element_blank(),
          # strip.text = element_text(size = 10, face = "bold"), 
          axis.text.y = element_text(size =13),  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(), 
          axis.line.x = element_line(color = "black"),  
          axis.text.x = element_text(size =13),  
          axis.title.x = element_text(size = 15),  
          axis.ticks.y=element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)
        ) +
        scale_y_continuous(position = "right",expand = c(0, 0))  
      
      
      
      ggsave(pbs,filename=file.path(path,pro,'gene_up&down_go_BP_treebarplot.png'),scale=1,height=12,width=12)
      ggsave(pbs2,filename=file.path(path,pro,'gene_up&down_go_BP_treebarplot2.png'),scale=1,height=12,width=12)
      ggsave(pbs2,filename=file.path(path,pro,'gene_up&down_go_BP_treebarplot2.pdf'),height=12,width=12,units = "in", dpi = 300)
    }, error = function(e) {
      cat("Cannot paint treebarplot-gobp-updw: ", e$message, "\n")
    })
    
  }
  
  if (exists('gene_up3')& exists('gene_dw3')){
    goup <- enrichGO(gene_up3, OrgDb = "org.Hs.eg.db", ont="BP",
                   pAdjustMethod = "BH",readable= TRUE,
                   pvalueCutoff = 0.3,
                   qvalueCutoff =0.5) 
    goup <- clusterProfiler::simplify(goup, cutoff=0.7, by="p.adjust", select_fun=min)
    goup<-DOSE::setReadable(goup, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    goup <- goup@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    if(nrow(goup)>=1){goup$group=1}
    
    godw <- enrichGO(gene_dw3, OrgDb = "org.Hs.eg.db", ont="BP",
                     pAdjustMethod = "BH",readable= TRUE,
                     pvalueCutoff = 0.3,
                     qvalueCutoff =0.5) 
    godw <- clusterProfiler::simplify(godw, cutoff=0.7, by="p.adjust", select_fun=min)
    godw<-DOSE::setReadable(godw, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    godw<- godw@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    if(nrow(godw)>=1){godw$group=-1}
    
    
    dat=rbind(goup,godw)
    dat<-dat %>% subset(qvalue<0.2) %>% arrange(desc(RichFactor))
    if (nrow(dat)>=20){
      dat<- dat %>% head(20) 
    }
    
    dat$logPvalue = -log10(dat$pvalue)
    dat$logPvalue=dat$logPvalue*dat$group 
    dat=dat[order(dat$logPvalue,decreasing = F),]
    if(nrow(goup)>=1){datup<-dat %>% filter(group=='1')}
    if(nrow(godw)>=1){datdw<-dat %>% filter(group=='-1')}
    
    tryCatch({
    pe <- ggplot(dat, aes(x=`logPvalue`, y=reorder(Description,order(logPvalue, decreasing = F)), fill=`logPvalue`)) +
      geom_col(aes(fill = `logPvalue`), width = 0.1)+
      geom_point(aes(size = RichFactor,
                     color = `logPvalue`))+
      scale_size_continuous(range = c(2, 7)) +
      geom_text(data = datup,
                aes(x = -0.2, y = reorder(Description,order(logPvalue, decreasing = F)), label = reorder(geneID,order(logPvalue, decreasing = F))),
                size = 3.5,
                hjust = 1)+ 
      geom_text(data = datdw,
                aes(x = 0.2, y = reorder(Description,order(logPvalue, decreasing = F)), label = reorder(geneID,order(logPvalue, decreasing = F))),
                size = 3.5,
                hjust = 0)+ 
      scale_color_continuous_c4a_div('sunset', mid = 0, reverse = F) +
      scale_fill_continuous_c4a_div('sunset', mid = 0, reverse = F) +
      scale_x_continuous(limits = c(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue))),
                         breaks = seq(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue)),by = 2),
                         labels = abs(seq(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue)),by = 2))
      ) + 
      labs(x = '-log(Pvalue)', y = 'GO terms')+
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 13),
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 12) 
      )+
      guides(
        fill = "none", 
        color = "none"  
      ) +
      theme_classic()
    ggsave(pe,filename=file.path(path,pro,'gene_up&down_go_BP_loliplot2.png'),scale=1,height=8,width=14)
    }, error = function(e) {
      cat("Cannot paint loliplot-gobp-updw: ", e$message, "\n")
    })
  }
  print('GO analysis is finished!')
   
}