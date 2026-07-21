run_GSEA_REA <- function(deg, path,pro,check=NULL){
  library(myenrichplot)
  library(ggstatsplot)
  library(clusterProfiler)
  library(patchwork)
  library(ReactomePA)
  df <- bitr(unique(deg$name), fromType = "ENSEMBL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  DEG<-deg
  DEG<-merge(DEG,df,by.y='ENSEMBL',by.x='name')
  geneList=DEG$logFC
  names(geneList)=DEG$ENTREZID
  geneList=sort(geneList,decreasing = T)
  
  ra_gse <- gsePathway(geneList     = geneList,
                       organism = "human",
                  minGSSize    = 10,
                  verbose      = FALSE,
                  pvalueCutoff = 0.5)
  ra<-DOSE::setReadable(ra_gse, OrgDb='org.Hs.eg.db',
                        keyType='ENTREZID')
  write.csv(ra@result,file = file.path(path,pro,'gsea_rea.csv'))
  save(ra,file = file.path(path,pro,'gsea_rea.Rdata'))

  
  up_ra <- ra[head(order(ra$enrichmentScore,decreasing = T)),];up_ra$group=1
  dw_ra <- ra[tail(order(ra$enrichmentScore,decreasing = T)),];dw_ra$group=-1
  
  dat=rbind(up_ra,dw_ra)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]

  
  pe <- ggplot(dat, aes(x=`pvalue`, y=reorder(Description,order(pvalue, decreasing = F)), fill=`pvalue`)) +
    geom_col(aes(fill = `pvalue`), width = 0.1)+
    geom_point(aes(size = abs(NES),
                   color = `pvalue`))+
    scale_size_continuous(range = c(2, 7)) +
    scale_color_continuous_c4a_div('sunset', mid = 0, reverse = F) +
    scale_fill_continuous_c4a_div('sunset', mid = 0, reverse = F) +
    scale_x_continuous(limits = c(floor(min(dat$pvalue)), ceiling(max(dat$pvalue))),
                       breaks = seq(floor(min(dat$pvalue)), ceiling(max(dat$pvalue)),by = 2),
                       labels = abs(seq(floor(min(dat$pvalue)), ceiling(max(dat$pvalue)),by = 2))
    ) + 
    labs(x = 'log(Pvalue)', y = 'REA terms')+
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
  
  ggsave(pe,filename = file.path(path,pro,'gsea_rea_loliplot.png'),scale=1,height=8,width=14)
  
  peup <-myenrichplot::gseaplot2(ra, geneSetID = head(rownames(up_ra)), 
                                 title = 'Rea activated pathways',pvalue_table = T,
                                 color = ggsci::pal_lancet()(length(head(rownames(up_ra))))
  )  
  png(file=file.path(path,pro,'ra_up_gseaplot.png'),width=8,height = 6,units = 'in',res=300) 
  print(peup)
  dev.off()
  pedw <-myenrichplot::gseaplot2(ra, geneSetID = head(rownames(dw_ra)), 
                                 title = 'Rea inhibited pathways',pvalue_table = T,
                                 color = ggsci::pal_lancet()(length(head(rownames(dw_ra))))
  )
  png(file=file.path(path,pro,'ra_down_gseaplot.png'),width=8,height = 6,units = 'in',res=300) 
  print(pedw)
  dev.off()
  if (is.null(check)){
    ped <-myenrichplot::gseaplot2(ra, geneSetID = c(head(rownames(up_ra)),head(rownames(dw_ra))), 
                                  title = 'Rea affected pathways',pvalue_table = T,
                                  color = ggsci::pal_lancet()(length(c(head(rownames(up_ra)),head(rownames(dw_ra)))))
    )
    png(file=file.path(path,pro,'rea_updown_gseaplot.png'),width=15,height =15,units = 'in',res=300)
    print(ped)
    dev.off()
  }
  else{
    tryCatch({
      ped <-myenrichplot::gseaplot2(ra, geneSetID = ra[grepl(check,ra@result$Description),'ID'], 
                                    title = paste0('Rea affected_', check,'_pathways'),pvalue_table = T,
                                    # color = ggsci::pal_lancet()(length(c(head(rownames(up_g)),head(rownames(down_g)))))
      )
      png(file=file.path(path,pro,paste0('rea_updown_',check,'_gseaplot.png')),width=15,height =15,units = 'in',res=300)
      print(ped)
      dev.off()
    }, error = function(e) {
      cat("Cannot paint gseaplot2 for check term: ", e$message, "\n")
    })
  }
  
  print('GSEA REA analysis is finished!')
}