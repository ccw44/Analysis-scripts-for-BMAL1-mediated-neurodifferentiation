run_updwclsuter <- function(path,pro){
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  library(myenrichplot)
  
  load(file.path(path,pro,'DEG_FC1.Rdata'))
  
  deg<-df
  
  df <- bitr(unique(deg$name), fromType = "ENSEMBL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  DEG<-deg
  DEG<-merge(DEG,df,by.y='ENSEMBL',by.x='name')
  
  ## |logFC|>1
  gene_up= unique(DEG[DEG$g == 'up','ENTREZID'])
  gene_down=unique(DEG[DEG$g == 'down','ENTREZID'])
  gene_diff=unique(c(gene_up,gene_down))
  
  ifelse(!dir.exists(file.path(path,pro)),
         dir.create(file.path(path,pro)),
         "Directory Exists")

  for (i in c("enrichKEGG", 'enrichGO', 'enrichPathway')) {
    if (i == 'enrichKEGG') {
      xx <- compareCluster(geneCluster =list(up=gene_up, down=gene_down), 
                           fun=i, organism="hsa", pvalueCutoff=0.3)
    } else if (i == 'enrichGO') {
      xx <- compareCluster(geneCluster =list(up=gene_up, down=gene_down), 
                           fun=i, OrgDb="org.Hs.eg.db", ont='BP', pvalueCutoff=0.3)
    } else {
      xx <- compareCluster(geneCluster =list(up=gene_up, down=gene_down),
                           fun=i, organism='human', pvalueCutoff=0.3)
    }
    x1 <- dotplot(xx) + 
          theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
                axis.text.y=element_text(size=11),
                legend.title = element_text(size = 11),
                legend.text = element_text(size = 10),) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=25))+
          labs(x = "")
    ggsave(x1, filename=file.path(path,pro, paste0('updown_cluster_', i, '_dotplot.png')), scale=1, width=7, height=7)
    ggsave(x1,filename=file.path(path,pro, paste0('updown_cluster_', i, '_dotplot.pdf')),width = 7, height =7, units = "in", dpi = 300)
    
    print(paste0(i,' is finished!'))
  }
  
  
  for (i in c("KEGG", 'GO', 'Rea')) {
    if (i == 'KEGG') {
      kk<- enrichKEGG(gene      = gene_diff,
                            organism     = 'hsa',
                            #universe     = gene_all,
                            pvalueCutoff = 0.9,
                            qvalueCutoff =0.9)
      kk<-DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
      kk@result <- kk@result %>% filter(p.adjust<0.05 & qvalue<0.05)
      
      # kk1<-kk@result %>% filter(subcategory =='Nervous system')
      
      kk@result$up <- 0
      kk@result$down <- 0
      for (i in 1:nrow(kk@result)) {
        genes_term <- unlist(strsplit(kk@result[i, "geneID"], "/")) 
        
        up_genes <- genes_term[genes_term %in% DEG$SYMBOL[DEG$g == "up"]]
        down_genes <- genes_term[genes_term %in% DEG$SYMBOL[DEG$g == "down"]]
        
        kk@result$up[i] <- length(up_genes)
        kk@result$down[i] <- length(down_genes)
      }
      
      kk_long <- kk@result %>%
        dplyr::arrange(desc(RichFactor)) %>%
        dplyr::slice_head(n = 20) %>%  
        dplyr::select(category, subcategory,Description, up, down) %>%
        pivot_longer(cols = c(up, down), names_to = "Updown", values_to = "GeneCount")
      kk_long$category<-factor(kk_long$category,levels=c('Metabolism',
                                                         'Organismal Systems',
                                                         "Human Diseases",
                                                         'Genetic Information Processing',
                                                         'Cellular Processes',
                                                         'Environmental Information Processing'
      ))
      kk_long$Description <- factor(kk_long$Description,level=unique(kk_long$Description))
      kk_long <-kk_long %>% arrange(category)
      kk_long$Description <- factor(kk_long$Description,
                                       levels = unique(kk_long$Description))
      barl2<-ggplot(kk_long, aes(x = Description, y = GeneCount, fill = Updown,color = category)) +
        geom_bar(stat = "identity", position = "dodge") + 
        scale_fill_manual(values = c("up" = "orange", "down" = "skyblue")) +  
               scale_x_discrete(labels=function(x) str_wrap(x, width=20))+ 
        scale_color_manual(values = c("Metabolism" = "steelblue", 
                                      "Organismal Systems" = "red", 
                                      "Human Diseases" = "green", 
                                      "Genetic Information Processing" = "purple", 
                                      "Cellular Processes" = "cyan",
                                      "Environmental Information Processing" = "orange"),
                           guide = guide_legend(override.aes = list(fill = NA)))+
        scale_y_continuous(expand = c(0, 0),limits = c(0, max(kk_long$GeneCount) * 1.1)) + 
        labs(title = "Level2 KEGG term", x = "KEGG Term", y = "Number of Genes") +
        theme_bw() + 
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
              axis.text = element_text(size = 12), 
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x =element_blank(),
              legend.position = "right",
              legend.key = element_blank(),
              panel.grid = element_blank())
      ggsave(barl2,filename = file.path(path,pro,'gene_diff_kegg_level2barplot(top20).png'),scale=1,height=8,width=16)
      ggsave(barl2,filename = file.path(path,pro,'gene_diff_kegg_level2barplot(top20).pdf'),height=8,width=16,units = "in", dpi = 300)
    } else if (i == 'GO') {
        go <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="all",pAdjustMethod = "BH",readable= TRUE,
                       pvalueCutoff = 0.2,
                       qvalueCutoff =0.5) 
        go <- clusterProfiler::simplify(go, cutoff=0.3, by="p.adjust", select_fun=min)
        
        go<-DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
        go@result <- go@result %>% filter(p.adjust<0.05 & qvalue<0.05)
        go@result$up <- 0
        go@result$down <- 0
      
        for (i in 1:nrow(go@result)) {
          genes_term <- unlist(strsplit(go@result[i, "geneID"], "/")) 
          
          up_genes <- genes_term[genes_term %in% DEG$SYMBOL[DEG$g == "up"]]
          down_genes <- genes_term[genes_term %in% DEG$SYMBOL[DEG$g == "down"]]
          
          go@result$up[i] <- length(up_genes)
          go@result$down[i] <- length(down_genes)
        }
      
        go_long <- go@result %>%
          dplyr::select(ONTOLOGY, Description, up, down) %>%
          pivot_longer(cols = c(up, down), names_to = "Updown", values_to = "GeneCount")
        go_long$ONTOLOGY<-factor(go_long$ONTOLOGY,level=c('BP','MF','CC'))
        go_long$Description <- factor(go_long$Description,level=unique(go_long$Description))
        go_long <-go_long %>% arrange(ONTOLOGY)
        barl2<-ggplot(go_long, aes(x = Description, y = GeneCount, fill = Updown,color = ONTOLOGY)) +
          geom_bar(stat = "identity", position = "dodge") +
          scale_fill_manual(values = c("up" = "orange", "down" = "skyblue")) + 
          scale_x_discrete(labels=function(x) str_wrap(x, width=20))+ 
          scale_color_manual(values = c("BP" = "black", 
                                        "MF" = "green", 
                                        "CC" = "pink"),
                             guide = guide_legend(override.aes = list(fill = NA)))+
          labs(title = "Level2 GO term", x = "GO Term", y = "Number of Genes") +
          theme_bw() + 
          theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
                axis.text = element_text(size = 10),  
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "right",
                legend.key = element_blank(),
                panel.grid = element_blank())
        ggsave(barl2,filename = file.path(path,pro,'gene_diff_go_level2barplot.png'),scale=1,height=8,width=14)
    } else {
      
    }
  }
}