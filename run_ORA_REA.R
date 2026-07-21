run_ORA_REA <- function(deg, path,pro,
                        type=NULL,ct=F,
                        type2=NULL,check=NULL,m6a1=F){
  # type: used for merely up or dw genes enrichment
  # ct：whether use CUTtag data
  
  library(BioEnricher)
  
  
  reaplot2 <- function(rearesult,pattern){
    rearesult<-rearesult@result %>% filter(pvalue<0.05 & p.adjust<0.1)%>% arrange(pvalue)
    if (nrow(rearesult)>=30){
      rearesult<- rearesult %>% head(30) 
    }
    rearesult$Description <- factor(rearesult$Description,
                                     levels = rev(unique(rearesult$Description)))
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
      axis.ticks.y = element_blank() # 去除 y 轴刻度线
    )

    ggplot(data = rearesult, aes(x = -log10(pvalue), y = Description,fill = EnrichmentFactor)) +
      scale_fill_distiller(palette = pattern,direction = 1)+ 

      geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
      scale_x_continuous(expand=c(0,0))+
      labs(x = "-log10(pvalue)", y = "", title = "Reactome pathways") +
      geom_text(aes(x = 0.1,
                    label = Description),
                hjust = 0)+
      theme_classic() + 
      mytheme
  }
  if (!is.null(type)& isTRUE(ct)){
    deg<-deg %>% filter(!is.na(CT)) 
  }
  df <- bitr(unique(deg$name), fromType = "ENSEMBL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  # head(df)
  DEG<-deg
  # head(DEG)
  DEG<-merge(DEG,df,by.y='ENSEMBL',by.x='name')

  
  ## |logFC|>1
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
  bora_up<- lzq_ORA(
    genes = gene_up,
    gene.type = "ENTREZID",
    enrich.type = 'Reactome',
    organism = "Human",
    pvalue.cutoff  = 0.1,
    qvalue.cutoff  = 0.5
  )

  ra_up<-reaplot2(bora_up,pattern='RdPu')

  ggsave(ra_up,filename=file.path(path,pro,'gene_up_Rea_barplot.png'),scale=1,width=8,height=15)
  
  write.csv(bora_up,file = file.path(path,pro,'gene_up_Rea.csv'),row.names = F)
  }
  
  if (exists('gene_down')){
  bora_down<- lzq_ORA(
    genes = gene_down,
    gene.type = "ENTREZID",
    enrich.type = 'Reactome',
    organism = "Human",
    pvalue.cutoff  = 0.1,
    qvalue.cutoff  = 0.5
  )
  

  tryCatch({
  ra_down<-reaplot2(bora_down,pattern='PuBuGn')
  ggsave(ra_down,filename=file.path(path,pro,'gene_down_Rea_barplot.png'),scale=1,width=8,height=15)
  }, error = function(e) {
    cat("Cannot paint barplot-radw: ", e$message, "\n")
  })
  
  write.csv(bora_down,file = file.path(path,pro,'gene_down_Rea.csv'),row.names = F)
  
  }
  
  if (exists('gene_diff')){
    
    bora_uptop<- lzq_ORA(
      genes = uptop,
      gene.type = "ENTREZID",
      enrich.type = 'Reactome',
      organism = "Human",
      pvalue.cutoff  = 0.05,
      qvalue.cutoff  = 0.05
    )
    tryCatch({

      ra_uptop<-reaplot2(bora_uptop,pattern='RdPu')
      ggsave(ra_uptop,filename=file.path(path,pro,'gene_uptop_Rea_barplot.png'),scale=1,width=8,height=15)
    }, error = function(e) {
      cat("Cannot paint treeplot: ", e$message, "\n")
    })
    
    write.csv(bora_uptop,file = file.path(path,pro,'gene_uptop_Rea.csv'),row.names = F)
    
    bora_dwtop<- lzq_ORA(
      genes = dwtop,
      gene.type = "ENTREZID",
      enrich.type = 'Reactome',
      organism = "Human",
      pvalue.cutoff  = 0.05,
      qvalue.cutoff  = 0.05
    )
    tryCatch({

      ra_dwtop<-reaplot2(bora_dwtop,pattern='PuBuGn')
      ggsave(ra_dwtop,filename=file.path(path,pro,'gene_dwtop_Rea_barplot.png'),scale=1,width=8,height=15)
    }, error = function(e) {
      cat("Cannot paint treeplot: ", e$message, "\n")
    })
    write.csv(bora_dwtop,file = file.path(path,pro,'gene_dwtop_Rea.csv'),row.names = F)
    
  }

  
  if (exists('gene_diff') && (nrow(bora_up)>0 & nrow(bora_down)>0) ){
  ra_updown<-lzq_ORA.barplot2(
    enrich.obj1 = bora_up,
    enrich.obj2 = bora_down,
    obj.types = c('Up','Down'),
  )+scale_y_discrete(labels=function(x)str_wrap(x,width=50))
  ggsave(ra_updown,filename=file.path(path,pro,'gene_up&down_Rea_barplot.png'),width=15,height=8)
  
  
  
  
  bodw<-bora_down@result
  if(nrow(bodw)>=1){bodw$group=-1 }
  boup<-bora_up@result
  if(nrow(boup)>=1){boup$group=1 }
  
  dat=rbind(boup,bodw)
  
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
  
  
  dat$logPadj = -log10(dat$p.adjust)
  dat$logPadj=dat$logPadj*dat$group 
  dat=dat[order(dat$logPadj,decreasing = F),]
  if(nrow(boup)>=1){datup<-dat %>% filter(group=='1')} 
  if(nrow(bodw)>=1){datdw<-dat %>% filter(group=='-1')}
  
  tryCatch({
  pe <- ggplot(dat, aes(x=`logPadj`, y=reorder(Description,order(logPadj, decreasing = F)), fill=`logPadj`)) +
    geom_col(aes(fill = `logPadj`), width = 0.1)+
    geom_point(aes(size = EnrichmentFactor,
                   color = `logPadj`))+
    scale_size_continuous(range = c(2, 7)) +
    geom_text(data = datup,
              aes(x = -0.2, y = reorder(Description,order(logPadj, decreasing = F)), label = reorder(str_wrap(Description, width=30),order(logPadj, decreasing = F))),
              size = 4.5,
              hjust = 1)+ 
    geom_text(data = datdw,
              aes(x = 0.2, y = reorder(Description,order(logPadj, decreasing = F)), label = reorder(str_wrap(Description, width=30),order(logPadj, decreasing = F))),
              size = 4.5,
              hjust = 0)+
    scale_color_continuous_c4a_div('sunset', mid = 0, reverse = F) +
    scale_fill_continuous_c4a_div('sunset', mid = 0, reverse = F) +
    
    scale_x_continuous(limits = c(-ceiling(max(abs(dat$logPadj))), ceiling(max(abs(dat$logPadj)))),
                       breaks = seq(-ceiling(max(abs(dat$logPadj))), ceiling(max(abs(dat$logPadj))),by = 1),
                       labels = abs(seq(-ceiling(max(abs(dat$logPadj))), ceiling(max(abs(dat$logPadj))),by = 1))
    ) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
    labs(x = '-log(Padj)', y = 'Reactome terms')+
    theme_classic()+
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(size = 14),
      axis.title = element_text(size = 15), 
      legend.title = element_text(size = 13), 
      legend.text = element_text(size = 12), 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank() 
    )+
    guides(
      fill = "none",  
      color = "none"  
    )
  ggsave(pe,filename=file.path(path,pro,'gene_up&down_Rea_loliplot.png'),scale=1,height=14,width=13)
  
  ggsave(pe,filename=file.path(path,pro,'gene_up&down_Rea_loliplot.pdf'),units = "in", dpi = 300,height=14,width=13)
    }, error = function(e) {
      cat("Cannot paint Rea_loliplot: ", e$message, "\n")
    })

  }
  if (is.null(type) && exists('gene_diff')){
    save(bora_up,bora_down,bora_uptop,bora_dwtop,
         file =file.path(path,pro,'Rea_results.Rdata'))}
  
  if (is.null(type) && !exists('gene_diff')){
    save(bora_up,bora_down,
         file =file.path(path,pro,'Rea_results.Rdata'))}
  
  if (!is.null(type) && grepl('up',type)){
    save(bora_up,file =file.path(path,pro,paste(type,'Rea_results.Rdata',sep='-')))}
  if (!is.null(type) &&grepl('dw',type)){
    save(bora_down,file =file.path(path,pro,paste(type,'Rea_results.Rdata',sep='-')))}
  
  
  if (exists('gene_up')& exists('gene_down')){
    
    boup2<-lzq_ORA(
      genes = gene_up,
      gene.type = "ENTREZID",
      enrich.type = 'Reactome',
      organism = "Human",
      pvalue.cutoff  = 0.05,
      qvalue.cutoff  = 0.5
    )
    
    boup2 <- boup2@result %>% filter(p.adjust<0.1 & qvalue<0.2)
    if(nrow(boup2)>=1){boup2$group=1 }
    
    bodw2<-lzq_ORA(
      genes = gene_down,
      gene.type = "ENTREZID",
      enrich.type = 'Reactome',
      organism = "Human",
      pvalue.cutoff  = 0.05,
      qvalue.cutoff  = 0.5
    )
    bodw2 <- bodw2@result %>% filter(p.adjust<0.1 & qvalue<0.2)
    if(nrow(bodw2)>=1){bodw2$group=-1 }
    
    dat=rbind(boup2,bodw2)
    dat<-dat %>% arrange(desc(EnrichmentFactor))
    
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
    if(nrow(boup)>=1){datup<-dat %>% filter(group=='1')} 
    if(nrow(bodw)>=1){datdw<-dat %>% filter(group=='-1')}
    
    if(!is.null(type)){dat<-dat[grepl(check,dat$Description),]}
    
    tryCatch({
      pe <- ggplot(dat, aes(x=`logPvalue`, y=reorder(Description,order(logPvalue, decreasing = F)), fill=`logPvalue`)) +
        geom_col(aes(fill = `logPvalue`), width = 0.1)+
        geom_point(aes(size = EnrichmentFactor,
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
        labs(x = '-log(Pvalue)', y = 'Rea terms')+
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
      ggsave(pe,filename=file.path(path,pro,'gene_up&down_Rea_loliplot2.png'),scale=1,height=8,width=14)
    }, error = function(e) {
      cat("Cannot paint Rea_loliplot2: ", e$message, "\n")
    })
  }
  if (exists('gene_up2')& exists('gene_dw2')){
    
    boup<-lzq_ORA(
      genes = gene_up2,
      gene.type = "ENTREZID",
      enrich.type = 'Reactome',
      organism = "Human",
      pvalue.cutoff  = 0.2,
      qvalue.cutoff  = 0.5
    )
    
    boup <- boup@result %>% filter(p.adjust<0.2)
    if(nrow(boup)>=1){boup$group=1 }
    
    bodw<-lzq_ORA(
      genes = gene_dw2,
      gene.type = "ENTREZID",
      enrich.type = 'Reactome',
      organism = "Human",
      pvalue.cutoff  = 0.2,
      qvalue.cutoff  = 0.5
    )
    bodw <- bodw@result %>% filter(p.adjust<0.2)
    if(nrow(bodw)>=1){bodw$group=-1 }

    dat=rbind(boup,bodw)
    dat<-dat %>% arrange(desc(EnrichmentFactor))
    if (nrow(dat)>=20){
      dat<- dat %>% head(20) 
    }
    
    dat$logPvalue = -log10(dat$pvalue)
    dat$logPvalue=dat$logPvalue*dat$group 
    dat=dat[order(dat$logPvalue,decreasing = F),]
    if(nrow(boup)>=1){datup<-dat %>% filter(group=='1')} 
    if(nrow(bodw)>=1){datdw<-dat %>% filter(group=='-1')}
    
    dat <- dat %>%
      mutate(
        Description = case_when(
          group == 1 & duplicated(Description) ~ paste0(Description, "_up"),
          group == -1 & duplicated(Description) ~ paste0(Description, "_down"),
          TRUE ~ Description  
        )
      )
    tryCatch({
    pe <- ggplot(dat, aes(x=`logPvalue`, y=reorder(Description,order(logPvalue, decreasing = F)), fill=`logPvalue`)) +
      geom_col(aes(fill = `logPvalue`), width = 0.1)+
      geom_point(aes(size = EnrichmentFactor,
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
      labs(x = '-log(Pvalue)', y = 'Rea terms')+
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
    ggsave(pe,filename=file.path(path,pro,'gene_up&down_Rea_loliplot2.png'),scale=1,height=8,width=14)
    }, error = function(e) {
      cat("Cannot paint Rea_loliplot2: ", e$message, "\n")
    })
  }
  
  print("Reactome analysis is finished!")
}