run_GSEA_KEGG <- function(deg, path,pro){
  library(myenrichplot)
  library(ggstatsplot)
  library(clusterProfiler)
  library(patchwork)
  df <- bitr(unique(deg$name), fromType = "ENSEMBL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  DEG<-deg
  DEG<-merge(DEG,df,by.y='ENSEMBL',by.x='name')
  geneList=DEG$logFC
  names(geneList)=DEG$ENTREZID
  geneList=sort(geneList,decreasing = T)
  head(geneList)
  #mmu
  #hsa
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 10,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  kk<-DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',
                       keyType='ENTREZID')
  write.csv(kk@result,file = file.path(path,pro,'gsea_kegg.csv'))
  save(kk,file = file.path(path,pro,'gsea_kk.Rdata'))
  
  up_k <- kk[head(order(kk$enrichmentScore,decreasing = T)),];up_k$group=1
  down_k <- kk[tail(order(kk$enrichmentScore,decreasing = T)),];down_k$group=-1
  
  dat=rbind(up_k,down_k)
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
    labs(x = 'log(Pvalue)', y = 'KEGG terms')+
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
  
  ggsave(pe,filename = file.path(path,pro,'kegg_gsea_loliplot.png'),scale=1,height=8,width=14)

  peup <-myenrichplot::gseaplot2(kk, geneSetID = head(rownames(up_k)), 
                          title = 'KEGG activated pathways',pvalue_table = T,
                          color = ggsci::pal_lancet()(length(head(rownames(up_k))))
                          )  
  png(file=file.path(path,pro,'kegg_up_gseaplot.png'),width=8,height = 6,units = 'in',res=300) 
  print(peup)
  dev.off()
  pedw <-myenrichplot::gseaplot2(kk, geneSetID = head(rownames(down_k)), 
                          title = 'KEGG inhibited pathways',pvalue_table = T,
                          color = ggsci::pal_lancet()(length(head(rownames(down_k))))
                          )
  png(file=file.path(path,pro,'kegg_down_gseaplot.png'),width=8,height = 6,units = 'in',res=300) 
  print(pedw)
  dev.off()
  
  ped <-myenrichplot::gseaplot2(kk, geneSetID = c(head(rownames(up_k)),head(rownames(down_k))), 
                                 title = 'KEGG affected pathways',pvalue_table = T,
                                 color = ggsci::pal_lancet()(length(c(head(rownames(up_k)),head(rownames(down_k)))))
                                )
  png(file=file.path(path,pro,'kegg_updown_gseaplot.png'),width=15,height =15,units = 'in',res=300) 
  print(ped)
  dev.off()

  print('GSEA KEGG analysis is finished!')
}