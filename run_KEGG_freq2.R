run_KE_highfreq2 <- function(kelist,
                            kelist2=NULL,
                            check,df,
                            path,pro,
                            KEGGtype,
                            limit=2, # screen data
                            f1=7,f2=4, # barplot_break
                            y1=15,y2=2 # y_limit & step
){
  # kelist2: only for down kelist to make updown barplot
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(tidyverse)
  library(cols4all)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  
  
  if(!is.null(kelist2)){
    gl1<-kelist@result %>% subset(qvalue<0.2)
    gl2<-kelist2@result %>% subset(qvalue<0.2)
    gl1$logP<--log10(gl1$pvalue)
    gl2$logP<--log10(gl2$pvalue)
    
    gl1<-gl1[grepl(check,gl1$Description),c('Description','logP','Count','geneID')]
    gl2<-gl2[grepl(check,gl2$Description),c('Description','logP','Count','geneID')]
    gl2$logP<--gl2$logP
    gl<-rbind(gl1,gl2)
    gl$group <- case_when(gl$logP > 0 ~ 'up',
                          gl$logP < 0 ~ 'down')
    
    if (length(gl$Description[duplicated(gl$Description)]) >0){
      gl <- gl %>%
        mutate(Description = if_else(duplicated(Description) | Description %in% duplicated(Description),
                                     case_when(
                                       group == 'up' ~ paste0(Description, '(up)'),
                                       group == 'down' ~ paste0(Description, '(dw)'),
                                       TRUE ~ Description
                                     ),
                                     Description
        ))
    }
    gl$Description <- factor(gl$Description,
                             levels = rev(gl$Description))
    up<- gl[which(gl$logP > 0),]
    down<- gl[which(gl$logP < 0),]
    pb <- ggplot(gl,aes(x = `logP`,y = Description,fill = `logP`)) +
      geom_col(aes(fill = `logP`), width = 0.1)+
      geom_point(aes(size = Count,
                     color = `logP`))+
      scale_size_continuous(range = c(2, 7)) +
      geom_text(data = up,
                aes(x = -0.2, y = Description, label = geneID),
                size = 3.5,
                hjust = 1)+ 
      geom_text(data = down,
                aes(x = 0.2, y = Description, label = geneID),
                size = 3.5,
                hjust = 0)+ 
      scale_color_continuous_c4a_div('sunset', mid = 0, reverse = F) +
      scale_fill_continuous_c4a_div('sunset', mid = 0, reverse = F) +
      scale_x_continuous(limits = c(floor(min(gl$logP)), ceiling(max(gl$logP))),
                         breaks = seq(floor(min(gl$logP)), ceiling(max(gl$logP)),by = 2),
                         labels = abs(seq(floor(min(gl$logP)), ceiling(max(gl$logP)),by = 2))
      ) + 
      labs(x = 'log(Pvalue)', y = 'KEGG terms')+
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 13), 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 12) 
      )+
      theme_classic()
    ggsave(filename=file.path(path,pro,paste(KEGGtype,'with godown',check,'loliplot.png',sep='_')),plot=pb,scale=1,height=8,width=14)
    
  }
  
  
  
  
  kelist<-kelist@result
  
  kelist<-kelist %>% dplyr::select(category,ID,Description,pvalue,p.adjust,FoldEnrichment,geneID) %>% 
    rename('category'='category','term'='Description','adj_pval'='p.adjust','genes'='geneID')
  
  checkterm<-kelist[grepl(check,kelist$genes),]
  checkgenes<-stringr::str_split(string = checkterm$genes,
                                 pattern = "/",simplify = TRUE) %>%as.vector()
  checkgenes_tb<-sort(table(checkgenes),decreasing = T)
  checkgenes_tb<-checkgenes_tb[-1] # 排序后NULL值位于1，故删去
  checkgenes_tb<-as.data.frame(checkgenes_tb)
  checkgenes_tbar<-subset(checkgenes_tb,Freq>limit)
  hc<-ggplot(checkgenes_tbar,aes(x=checkgenes,y=Freq))+
    geom_bar(stat='identity', position='identity', 
             fill=ifelse(checkgenes_tbar$Freq>=f1,'#e474ab',
                         ifelse(checkgenes_tbar$Freq>=f2,'#64e259','grey80')),
             width=0.5)+
    scale_y_continuous(expand = c(0,0),
                       # all
                       # breaks=seq(0,25,2),
                       # limits=c(0,26)
                       # up
                       # breaks=seq(0,20,2),
                       # limits=c(0,20)
                       # down
                       breaks=seq(0,max(checkgenes_tbar$Freq),y2),
                       limits=c(0,max(checkgenes_tbar$Freq)+0.5)
    ) +
    labs(x = paste0('Genes in KEGG terms_',check), y = 'Freq')+
    theme_classic()+
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(color = "black",
                                     linewidth = 0.6),
          axis.title =  element_text(size=15,face = "bold"),
          axis.text.x = element_text(angle=45,  
                                     hjust =1,
                                     size=10),
          axis.text.y= element_text(size=12)
    )

  ## logFC top + freq + barplot term genes
  checkgenes_tb2<-checkgenes_tb
  names(checkgenes_tb2)[1]<-'SYMBOL'
  checkgenes_tb2<-merge(checkgenes_tb2,df,by.x='SYMBOL')
  
  ff<-ggplot(checkgenes_tb2, aes(x = Freq, y = abs(logFC), label = SYMBOL)) +
    geom_point(aes(color = ifelse(Freq > 4 & abs(logFC) > 2,"#26e701",
                                  ifelse(Freq > 4 | abs(logFC) > 2, "red", "blue"))), size = 3)+  # Conditional point color
    geom_text_repel(data = subset(checkgenes_tb2, Freq > 4 | abs(logFC) > 2), 
                    vjust = -0.5, hjust = 1.5, size = 4, family = 'serif',
                    min.segment.length =Inf) +  # Conditional text labels
    geom_vline(xintercept = 4, linetype = "dashed", color = "grey") +  # Vertical dashed line
    geom_hline(yintercept = 2, linetype = "dashed", color = "grey") +  # Horizontal dashed line
    scale_color_identity() +  # Use the specified colors directly
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(checkgenes_tbar$Freq), 1), limits = c(0, max(checkgenes_tbar$Freq)+0.5)) +
    
    labs(x = "Frequency", y =ifelse(grepl('down','KEGG'),"-logFC","logFC") , 
         # title = "Invasion-related Genes"
    ) +
    theme_classic() +
    theme(
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(hjust = 1, size = 10),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)  # Center the title
    )

  ggsave(filename=file.path(path,pro,paste0(KEGGtype,"_",check,'_scatterplot.png')),plot=ff,scale=1,height=8,width=14)
  
  library(ggplot2)
  library(ggsankey)
  checkgenes2<-unique(checkgenes)[grepl(check,unique(checkgenes))]
  checkgenes2<-sort(checkgenes2)
  checkterm$Enrichindex <- cut(
    -log10(checkterm$adj_pval),
    breaks = quantile(-log10(checkterm$adj_pval), probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)+ c(0, 1e-10, 2e-10),
    labels = c("low", "medium", "high"),
    include.lowest = TRUE
  )
  checkterm$key<-'Specific function'
  checkterm<- checkterm %>% arrange(desc(pvalue))
  checkterm$term<-fct_inorder(checkterm$term)
  tb<-data.frame()
  for (i in checkgenes2) {
    tb1<-subset(checkterm[grepl(i,checkterm$genes),],
                select=c('genes','key','Enrichindex','term','pvalue','FoldEnrichment'))
    tb1$genes<-i
    tb<-rbind(tb,tb1)
  } 
  
  checkgenes_tb3<-checkgenes_tb2 %>% subset(select=c('SYMBOL','logFC','Freq'))
  tb<-merge(tb,checkgenes_tb3,by.x='genes',by.y='SYMBOL')

  tb$Enrichindex<-factor(tb$Enrichindex,level=c( "high", "medium","low"))

  tb2<-tb
  tb2$term <- factor(tb2$term, levels = checkterm$term)
  tb2$logFC<-round(tb2$logFC,2)
  df2 <- tb2 %>%
    make_long(logFC,Freq,genes, term, Enrichindex, key)
  df2.1 <- tb2 %>%
    make_long(logFC,Freq,genes,term)
  df2.1$node <- factor(df2.1$node, levels = c(unique(tb2$logFC),unique(tb2$Freq),unique(tb2$genes), levels(tb2$term)))
  
  
  # title<-paste(pro,GOtype,check,sep='_')
  names(df2)[1]<-'Var'
  s1<-ggplot(df2, aes(x =Var, 
                      next_x = next_x, 
                      node = node,
                      next_node = next_node,
                      fill = factor(node),
                      label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1) +
    geom_sankey_label(size = 4, color = 1, fill = "white",
                      hjust=0.5) + 
    scale_fill_viridis_d(option = "A", alpha = 0.95) +
    theme_sankey(base_size = 16) +
    theme(legend.position = "none")+ 
    labs(x = NULL)
  
  # for sankeydotplot
  names(df2.1)[1]<-'Var'
  s1.1<-ggplot(df2.1, aes(x =Var, 
                          next_x = next_x, 
                          node = node,
                          next_node = next_node,
                          fill = factor(node),
                          label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1) +
    geom_sankey_label(size = 4, color = 1, fill = "white",
                      hjust=1) + 
    scale_fill_viridis_d(option = "A", alpha = 0.95) +
    theme_sankey(base_size = 16) +
    theme(legend.position = "none")+ 
    labs(x = NULL)
  
  tb3 <- tb2 %>%
    group_by(term) %>%
    summarise(count = n(), .groups = 'drop') %>%
    left_join(checkterm %>% dplyr::select(term, pvalue,FoldEnrichment), by = "term")
  tb3$logP <- -log10(tb3$pvalue)
  
  tb3<- tb3 %>%
    mutate(ymax = cumsum(count)) %>% 
    mutate(ymin = ymax -count) %>%
    mutate(label =  cumsum(count) - count / 2)
  
  s2<-ggplot() +
    geom_point(data = tb3,
               aes(x = logP,
                   y= label,
                   size= FoldEnrichment,
                   color= FoldEnrichment)) +
    scale_size_continuous(range=c(2,8)) +
    scale_y_continuous(expand = c(0,0),limits = c(0,max(tb3$label+1))) +
    scale_colour_distiller(palette = "Reds", direction = 1) +
    labs(x = "-log10(Pvalue)",
         y= "") +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 11))
  
  
  s1.1<- s1.1 + theme(plot.margin = unit(c(0,10,0,0),units="cm"))
  s3<- ggdraw() + draw_plot(s1.1, 0, 0, 0.85, 1) + draw_plot(s2, scale = 1.1,0.55, 0.05, 0.4, 0.88)
  
  ggsave(filename=file.path(path,pro,paste0(KEGGtype,"_",check,'_sankeyplot.png')),plot=s1,width = 12, height = 8, scale = 1)
  ggsave(filename=file.path(path,pro,paste0(KEGGtype,"_",check,'_sankeyplot.pdf')),plot=s1,width = 12, height = 8, units = "in", dpi = 300)
  
  ggsave(filename=file.path(path,pro,paste0(KEGGtype,"_",check,'_sankeydotplot.png')),plot=s3,width = 16, height = 8, scale = 1)
  ggsave(filename=file.path(path,pro,paste0(KEGGtype,"_",check,'_sankeydotplot.pdf')),plot=s3,width = 16, height = 8, units = "in", dpi = 300)
  print("KEGG-freq plot is finished!")
}