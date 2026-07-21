run_GSEA_GO <- function(deg, path, pro, ont, check=NULL){
  library(myenrichplot)
  library(ggstatsplot)
  library(clusterProfiler)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  
  df <- bitr(unique(deg$name), fromType = "ENSEMBL",
             toType = c("ENTREZID"),
             OrgDb = org.Mm.eg.db)
  
  DEG <- merge(deg, df, by.y='ENSEMBL', by.x='name')
  
  DEG <- DEG %>%
    dplyr::filter(!is.na(ENTREZID)) %>%          
    dplyr::group_by(ENTREZID) %>%                
    dplyr::summarize(logFC = mean(logFC)) %>%    
    dplyr::ungroup()
  
  geneList = DEG$logFC
  names(geneList) = DEG$ENTREZID
  geneList = sort(geneList, decreasing = T)
  
  message(sprintf("Running GSEA GO (%s)...", ont))
  go_gse <- gseGO(geneList     = geneList,
                  OrgDb        = 'org.Mm.eg.db',
                  ont          = ont,
                  minGSSize    = 10,
                  pvalueCutoff = 1.0,  # 阈值放宽至 1.0 以获取背景
                  verbose      = FALSE)
  
  if (is.null(go_gse) || nrow(go_gse@result) == 0) {
    message(sprintf("\n[Excellent Safety]: No GO %s terms enriched at all!", ont))
    message("The expression profiles are nearly identical. Skipping plots to prevent crash.\n")
    
    dir.create(file.path(path, pro), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(Message="Extremely identical expression. No pathways enriched."), 
              file = file.path(path, pro, paste0('gsea_go_EMPTY_', ont, '.csv')), row.names = FALSE)
    
    return(invisible(NULL))
  }
  # ====================================================
  
  go <- DOSE::setReadable(go_gse, OrgDb='org.Mm.eg.db', keyType='ENTREZID')
  
  dir.create(file.path(path, pro), recursive = TRUE, showWarnings = FALSE)
  write.csv(go@result, file = file.path(path, pro, 'gsea_go.csv'))
  save(go, file = file.path(path, pro, paste0('gsea_go_', ont, '.Rdata')))
  
  res_df <- go@result
  
  up_g <- res_df %>% 
    dplyr::filter(enrichmentScore > 0) %>% 
    dplyr::arrange(desc(enrichmentScore)) %>% 
    head(6) %>% 
    dplyr::mutate(group = 1)
  
  down_g <- res_df %>% 
    dplyr::filter(enrichmentScore < 0) %>% 
    dplyr::arrange(enrichmentScore) %>% 
    head(6) %>% 
    dplyr::mutate(group = -1)
  
  dat = rbind(up_g, down_g)
  
  if(nrow(dat) == 0) {
    message("No pathways passed the positive/negative enrichment score split. Skipping plots.")
    return(invisible(NULL))
  }
  # ====================================================
  
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue = dat$pvalue * dat$group 
  dat = dat[order(dat$pvalue, decreasing = F),]
  
  pe <- ggplot(dat, aes(x=reorder(Description, order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="#34bfb5", high="#ff6633", guide = "none") + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + 
    theme_ggstatsplot()+
    theme(plot.title = element_text(size = 15,hjust = 0.5),  
          axis.text = element_text(size = 12,face = 'bold'),
          panel.grid = element_blank())+
    ggtitle("GO Pathway Enrichment") 
  
  ggsave(pe, filename = file.path(path, pro, paste0('gsea_go_', ont, '_barplot.png')), scale=1, height=8, width=14)
  
  if(nrow(up_g) > 0) {
    peup <- myenrichplot::gseaplot2(go, geneSetID = head(rownames(up_g)), 
                                    title = 'GO activated pathways', pvalue_table = T,
                                    color = ggsci::pal_lancet()(nrow(up_g)))
    png(file=file.path(path, pro, paste0('go_up_', ont, '_gseaplot.png')), width=8, height = 6, units = 'in', res=300) 
    print(peup)
    dev.off()
  }
  
  if(nrow(down_g) > 0) {
    pedw <- myenrichplot::gseaplot2(go, geneSetID = head(rownames(down_g)), 
                                    title = 'GO inhibited pathways', pvalue_table = T,
                                    color = ggsci::pal_lancet()(nrow(down_g)))
    png(file=file.path(path, pro, paste0('go_down_', ont, '_gseaplot.png')), width=8, height = 6, units = 'in', res=300) 
    print(pedw)
    dev.off()
  }
  
  print('GSEA GO analysis is finished or elegantly skipped due to high similarity!')
}