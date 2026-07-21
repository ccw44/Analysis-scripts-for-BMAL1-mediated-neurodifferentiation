run_ORA_KEGG <- function(deg, path,pro,
                         type=NULL,ct=F,m6a1=F){
  # type: used for merely up or dw genes enrichment
  # ct：whether use CUTtag data
  
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(scales)
  library(stringr)
  library(cols4all)
  

  if (!is.null(type)& isTRUE(ct)){
    deg<-deg %>% filter(!is.na(CT)) 
  }
  
  df <- bitr(unique(deg$name), fromType = "ENSEMBL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  DEG=deg
  DEG=merge(DEG,df,by.y='ENSEMBL',by.x='name')
  
  if (is.null(type)|(!is.null(type) && grepl('up',type))){gene_up= unique(DEG[DEG$g == 'up','ENTREZID'])}
  if (is.null(type) |(!is.null(type) && grepl('dw',type))){gene_down=unique(DEG[DEG$g == 'down','ENTREZID'])}
  if ((!is.null(type) && grepl('diff2',type))){
    gene_up2= unique(DEG[DEG$g == 'up','ENTREZID'])
    gene_dw2=unique(DEG[DEG$g == 'down','ENTREZID'])
  }
  if (is.null(type)&&isFALSE(m6a1)){gene_diff=unique(c(gene_up,gene_down))
                    uptop<-DEG %>% filter(g=='up') %>% arrange(desc(logFC)) %>% head(100) 
                    uptop<-uptop$ENTREZID
                    dwtop<-DEG %>% filter(g=='down') %>% arrange(logFC) %>% head(100) 
                    dwtop<-dwtop$ENTREZID
  }
  
  ifelse(!dir.exists(file.path(path,pro)),
         dir.create(file.path(path,pro)),
         "Directory Exists")
  
  if (exists('gene_up')){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.5,
                      qvalueCutoff =0.9)
  
  kk.up<-DOSE::setReadable(kk.up, OrgDb='org.Hs.eg.db',
                       keyType='ENTREZID')
  write.csv(kk.up@result, file = file.path(path,pro,'gene_up_KEGG.csv'))
  }
  
  if (exists('gene_down')){
  kk.down <- enrichKEGG(gene      = gene_down,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.5,
                        qvalueCutoff =0.9)
  
  
  kk.down=DOSE::setReadable(kk.down, OrgDb='org.Hs.eg.db',
                          keyType='ENTREZID')
  write.csv(kk.down@result, file = file.path(path,pro,'gene_down_KEGG.csv'))
  }
  
  if (exists('gene_diff')){
  kk.diff <- enrichKEGG(gene      = gene_diff,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.5,
                        qvalueCutoff =0.9)

  
  kk.diff=DOSE::setReadable(kk.diff, OrgDb='org.Hs.eg.db',
                          keyType='ENTREZID')
  write.csv(kk.diff@result, file =file.path(path,pro,'gene_diff_KEGG.csv'))
  
  kk.uptop <- enrichKEGG(gene      = uptop,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.5,
                        qvalueCutoff =0.9)
  
  
  kk.uptop=DOSE::setReadable(kk.uptop, OrgDb='org.Hs.eg.db',
                            keyType='ENTREZID')
  write.csv(kk.uptop@result, file =file.path(path,pro,'gene_uptop_KEGG.csv'))
  
  kk.dwtop <- enrichKEGG(gene      = dwtop,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.5,
                        qvalueCutoff =0.9)
  
  
  kk.dwtop=DOSE::setReadable(kk.dwtop, OrgDb='org.Hs.eg.db',
                            keyType='ENTREZID')
  write.csv(kk.dwtop@result, file =file.path(path,pro,'gene_dwtop_KEGG.csv'))
  }
  
  if (exists('gene_diff')){
    kegg_diff_dt <- as.data.frame(kk.diff)
    kegg_uptop_dt <- as.data.frame(kk.uptop)
    uptop_kegg<-kegg_uptop_dt[kegg_uptop_dt$pvalue<0.05,]
    kegg_dwtop_dt <- as.data.frame(kk.dwtop)
    dwtop_kegg<-kegg_dwtop_dt[kegg_dwtop_dt$pvalue<0.05,]
    }
  if (exists('gene_diff') && exists('gene_down')){
    kegg_down_dt <- as.data.frame(kk.down)
    down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.01,]
    if(nrow(down_kegg)>=1){down_kegg$group=-1}
  } else if (!is.null(type) && grepl('dw',type)){
    kegg_down_dt <- as.data.frame(kk.down)
    down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,]
    }
  if (exists('gene_diff') && exists('gene_up')){
    kegg_up_dt <- as.data.frame(kk.up)
    up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.01,]
    if(nrow(up_kegg)>=1){up_kegg$group=1}
  }else if (!is.null(type) && grepl('up',type)){
    kegg_up_dt <- as.data.frame(kk.up)
    up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,]
  }
  
  
  kegg_plot <- function(up_kegg,down_kegg){
    dat=rbind(up_kegg,down_kegg)
    colnames(dat)
    dat$pvalue = -log10(dat$pvalue)
    dat$pvalue=dat$pvalue*dat$group 
    
    dat=dat[order(dat$pvalue,decreasing = F),]
    
    g_kegg <- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 

      geom_bar(stat="identity") + 
      scale_fill_gradient(low="#34bfb5",high="#ff6633",guide = 'none') + 
      scale_x_discrete(name ="Pathway names") +
      scale_y_continuous(name ="log10P-value") +
      ggtitle("KEGG Enrichment") +
      coord_flip() + theme_bw()+
      theme(plot.title = element_text(size = 15,face = 'bold',hjust = 0.5),  
            axis.text = element_text(size = 12,face = 'bold'),
            panel.grid = element_blank())
    
  }
  kegg_plot2 <- function(keggresult){
    if (nrow(keggresult)>=30){
      keggresult<- keggresult %>% head(30) 
    }
    keggresult$category<-factor(keggresult$category,levels=c('Metabolism',
                                                             'Organismal Systems',
                                                             "Human Diseases",
                                                             'Genetic Information Processing',
                                                             'Cellular Processes',
                                                             'Environmental Information Processing'
    ))
    keggresult<-keggresult %>% filter(pvalue<0.01)%>% arrange(category, pvalue)
    keggresult$Description <- factor(keggresult$Description,
                                     levels = rev(keggresult$Description))

    
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
    keggresult$RichFactor <- scales::rescale(keggresult$RichFactor, to = c(0, 1))
    keggresult$RichFactor <- (rank(keggresult$RichFactor) - 1) / (nrow(keggresult) - 1)
    ggplot(data = keggresult, aes(x = -log10(pvalue), y = Description)) +
      geom_bar(aes(fill = category, alpha = RichFactor), stat = "identity", width = 0.8) +
      scale_fill_manual(values = c("Metabolism" = "steelblue", 
                                   "Organismal Systems" = "red", 
                                   "Human Diseases" = "green", 
                                   "Genetic Information Processing" = "purple", 
                                   "Cellular Processes" = "cyan",
                                   "Environmental Information Processing" = "orange")) +
      scale_alpha_continuous(range = c(0.1, 1),guide = "none") + # Adjust alpha range as needed
      scale_x_continuous(expand=c(0,0))+
      labs(x = "-log10(pvalue)", y = "KEGG terms") +
      geom_text(aes(x = -log10(pvalue), y = Description, label = round(FoldEnrichment, 1)),
                vjust = 0.2, hjust=-1, size = 6) +
      geom_text(aes(x = 0.03, label = Description), hjust = 0) +
      # geom_text(size=3,aes(x=0.1,label=geneID))+
      theme_classic() + 
      mytheme
  }
  if (exists('gene_diff')){
    g_kegg<-kegg_plot(up_kegg,down_kegg)
    ggsave(file.path(path,pro,"gene_diff_KEGG_barplot.png"),plot=g_kegg,scale=1,height=8,width=14)
    
    tryCatch({
    g_kegg<-kegg_plot2(kegg_uptop_dt)
    ggsave(file.path(path,pro,"gene_uptop_KEGG_barplot2.png"),plot=g_kegg,scale=1,height=8,width=14)
    
    g_kegg<-kegg_plot2(kegg_dwtop_dt)
    ggsave(file.path(path,pro,"gene_dwtop_KEGG_barplot2.png"),plot=g_kegg,scale=1,height=8,width=14)
    }, error = function(e) {
      cat("Cannot paint kegg_plot2: ", e$message, "\n")
    })
    
  }
  if (exists('gene_up')){
    g_kegg<-kegg_plot2(kegg_up_dt)
    ggsave(file.path(path,pro,"gene_up_KEGG_barplot2.png"),plot=g_kegg,scale=1,height=8,width=14)
    ggsave(file.path(path,pro,"gene_up_KEGG_barplot2.pdf"),plot=g_kegg,height=8,width=14,units = "in", dpi = 300)
    
  }
  if (exists('gene_down')){
    g_kegg<-kegg_plot2(kegg_down_dt)
    ggsave(file.path(path,pro,"gene_down_KEGG_barplot2.png"),plot=g_kegg,scale=1,height=8,width=14)
    ggsave(file.path(path,pro,"gene_down_KEGG_barplot2.pdf"),plot=g_kegg,height=8,width=14,units = "in", dpi = 300)
  }
  if (exists('gene_diff')){
    save(kk.up,kk.down,kk.diff,
         kk.uptop,kk.dwtop,
         file =file.path(path,pro,'kegg_results.Rdata'))
  }
  
  if (!is.null(type) && grepl('up',type)){
    save(kk.up,file =file.path(path,pro,paste(type,'kegg_results.Rdata',sep='-')))}
  if (!is.null(type) &&grepl('dw',type)){
    save(kk.down,file =file.path(path,pro,paste(type,'kegg_results.Rdata',sep='-')))}
  
  
  if (exists('gene_up')& exists('gene_down')){
    kkup2 <- enrichKEGG(gene      = gene_up,
                       organism     = 'hsa',
                       #universe     = gene_all,
                       pvalueCutoff = 0.5,
                       qvalueCutoff =0.9)
    kkup2<-DOSE::setReadable(kkup2, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    kkup2 <- kkup2@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    if(nrow(kkup2)>=1){kkup2$group=1}
    
    kkdw2 <- enrichKEGG(gene      = gene_down,
                       organism     = 'hsa',
                       #universe     = gene_all,
                       pvalueCutoff = 0.5,
                       qvalueCutoff =0.9)
    kkdw2<-DOSE::setReadable(kkdw2, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    kkdw2 <- kkdw2@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    if(nrow(kkdw2)>=1){kkdw2$group=-1}
    
    dat=rbind(kkup2,kkdw2)
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
    datup<-dat %>% filter(group=='1')
    datdw<-dat %>% filter(group=='-1')
    
    dat$FoldEnrichment<-as.character(round(dat$FoldEnrichment,3))
    datup$FoldEnrichment<-as.character(round(datup$FoldEnrichment,3))
    datdw$FoldEnrichment<-as.character(round(datdw$FoldEnrichment,3))
    pe <- ggplot(dat, aes(x=`logPvalue`, y=reorder(FoldEnrichment,order(logPvalue, decreasing = F)), fill=`logPvalue`)) +
      geom_col(aes(fill = `logPvalue`), width = 0.1)+
      geom_point(aes(size = RichFactor,
                     color = `logPvalue`))+
      scale_size_continuous(range = c(2, 7)) +
      geom_text(data = datup,
                aes(x = -0.2, y = reorder(FoldEnrichment,order(logPvalue, decreasing = F)), label = reorder(Description,order(logPvalue, decreasing = F))),
                size = 3.5,
                hjust = 1)+ 
      geom_text(data = datdw,
                aes(x = 0.2, y = reorder(FoldEnrichment,order(logPvalue, decreasing = F)), label = reorder(Description,order(logPvalue, decreasing = F))),
                size = 3.5,
                hjust = 0)+
      scale_color_continuous_c4a_div('sunset', mid = 0, reverse = F) +
      scale_fill_continuous_c4a_div('sunset', mid = 0, reverse = F) +
      scale_x_continuous(limits = c(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue))),
                         breaks = seq(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue)),by = 2),
                         labels = abs(seq(floor(min(dat$logPvalue)), ceiling(max(dat$logPvalue)),by = 2))
      ) + 
      labs(x = '-log(Pvalue)', y = 'KEGG terms')+
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
    ggsave(pe,filename=file.path(path,pro,'gene_up&down_KEGG_loliplot.png'),scale=1,height=8,width=14)
    ggsave(pe,filename=file.path(path,pro,'gene_up&down_KEGG_loliplot.pdf'),height=8,width=14,units = "in", dpi = 300)
    
  }
  
  if (exists('gene_up2')& exists('gene_dw2')){
    kkup <- enrichKEGG(gene      = gene_up2,
                     organism     = 'hsa',
                     #universe     = gene_all,
                     pvalueCutoff = 0.5,
                     qvalueCutoff =0.9)
    kkup<-DOSE::setReadable(kkup, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    kkup <- kkup@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    if(nrow(kkup)>=1){kkup$group=1}
    
    kkdw <- enrichKEGG(gene      = gene_dw2,
                       organism     = 'hsa',
                       #universe     = gene_all,
                       pvalueCutoff = 0.5,
                       qvalueCutoff =0.9)
    kkdw<-DOSE::setReadable(kkdw, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    kkdw <- kkdw@result %>% filter(p.adjust<0.2 & qvalue<0.2)
    if(nrow(kkdw)>=1){kkdw$group=-1}
    
    dat=rbind(kkup,kkdw)
    dat<-dat %>% subset(qvalue<0.2) %>% arrange(desc(RichFactor))
    if (nrow(dat)>=20){
      dat<- dat %>% head(20) 
    }
    
    dat$logPvalue = -log10(dat$pvalue)
    dat$logPvalue=dat$logPvalue*dat$group 
    dat=dat[order(dat$logPvalue,decreasing = F),]
    datup<-dat %>% filter(group=='1')
    datdw<-dat %>% filter(group=='-1')
    
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
      labs(x = '-log(Pvalue)', y = 'KEGG terms')+
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
    ggsave(pe,filename=file.path(path,pro,'gene_up&down_KEGG_loliplot.png'),scale=1,height=8,width=14)
  }
  
  print('KEGG analysis is finished!')
}