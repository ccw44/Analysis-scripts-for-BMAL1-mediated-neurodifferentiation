run_DESeq2DEG <- function(dat,grouplist,
                          target_gene,
                          pro, path,
                          case, control,
                          logFC=1,color1){
  library(DESeq2) 
  library(clusterProfiler)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(ggplotify)
  library(patchwork)
  dir.create(file.path(path,pro))
  
  colData <- data.frame(row.names=grouplist$sample, 
                        grouplist1=grouplist$group,
                        batch=grouplist$batch
  )
  colData<-subset(colData,grouplist1 %in% c(control,case))
  grouplist<-subset(grouplist,group %in% c(control,case))
  grouplist$group<-factor(grouplist$group,level=c(control,case))
  grouplist1 <- grouplist$group

  dat<-dat[match(rownames(colData), names(dat))]
  if (length(unique(colData$batch))>1){
    dds<- DESeqDataSetFromMatrix(countData = dat,
                                 colData = colData,
                                 design = ~ grouplist1+batch)
  } else {
  dds<- DESeqDataSetFromMatrix(countData = dat,
                               colData = colData,
                               design = ~ grouplist1)
  }
  dds$grouplist1<-relevel(factor(grouplist1),ref=control)
  dds2 <- DESeq(dds)
  table(grouplist1)
  
  tmp <- results(dds2,contrast=c("grouplist1",case,control))
  # tmp <- results(dds2)
  DEG_DESeq2<-as.data.frame(tmp[order(tmp$padj),])
  DEG_DESeq2<-na.omit(DEG_DESeq2)
  deg<-DEG_DESeq2
  colnames(deg)<-c("AveExpr", "logFC"  ,'lfcSE'   , 'stat',"P.Value"  , "adj.P.Val")
  deg<-deg[order(deg$logFC),] 
  head(deg)
  ### volcano 
  if(T){
    nrDEG<-deg
    head(nrDEG)
    attach(nrDEG)
    df<-nrDEG
    # cut_logFC <- with(df,mean(abs(df$logFC)) + 2*sd(abs(df$logFC)) )
    # if(cut_logFC > 1){logFC_t = 1}
    # if(cut_logFC < 1){logFC_t = cut_logFC}
    logFC_t <- logFC
    pvalue_t <- 0.05
    df$v<- -log10(adj.P.Val) 
    df$g<-ifelse(df$adj.P.Val > pvalue_t,'stable',
                 ifelse( df$logFC > logFC_t,'up', 
                         ifelse( df$logFC < -logFC_t,'down','stable') )
    )
    df$g <- factor(df$g, levels = c("up","down","stable"))
    df$name=rownames(df)
    bitrdf1 <- bitr(df$name,fromType="ENSEMBL",toType="SYMBOL",OrgDb ='org.Hs.eg.db')
    names(bitrdf1)[which(names(bitrdf1) =='ENSEMBL')]<-'name'
    df<-merge(df,bitrdf1,by='name',all.x=T)
    head(df)
    
    this_tile <- paste0('Threshold of logFC is ',round(logFC_t,3),
                        '\nThe number of up gene is ',nrow(df[df$g == 'up',]) ,
                        '\nThe number of down gene is ',nrow(df[df$g == 'down',])
    )
    for_label <-  df %>% 
      dplyr::filter(df$SYMBOL %in% target_gene)
    # p5
    p5 <- ggplot(data = df, 
                 aes(x = logFC, 
                     y = -log10(adj.P.Val))) +
      geom_point(alpha=0.6, size=1.5, 
                 aes(color=g)) +
      geom_point(size = 3, shape = 1, data = for_label) +
      ggrepel::geom_label_repel(aes(label = name),data = for_label,
                                max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                                color="black") +
      ylab("-log10(Padj)")+
      scale_color_manual(values=c("#ff6633","#34bfb5", "#828586"))+
      # geom_vline(xintercept= 0,lty=4,col="grey",lwd=0.8) +
      geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(pvalue_t),lty=4,col="black",lwd=0.8) +
      xlim(-5, 5)+
      theme_bw()+
      ggtitle(this_tile )+
      theme(panel.grid = element_blank(),
            plot.title = element_text(size=8,hjust = 0.5),
            legend.title = element_blank(),
            legend.text = element_text(size=8))
    # volcano p7
    for_label_up<-subset(for_label,g=='up')
    for_label_down<-subset(for_label,g=='down')
    for_label_none<-subset(for_label,g=='stable')
    
    mytheme <- theme_classic() +
      theme(axis.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            legend.position = 'none',
            plot.margin = margin(10,5,5,5) 
      )
    mycol <-c("#EB4232","#a4e1da","#d8d8d8")
    p7 <- ggplot(data = df, 
                 aes(x = logFC, 
                     y = -log10(adj.P.Val),
                     color = g)) + 
      geom_point(size = 1.5) + 
      scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + 
      scale_x_continuous(limits = c(-16, 16),
                         breaks = seq(-15, 15, by = 5)) + 
      scale_y_continuous(expand = expansion(add = c(0.5, 0.5)),
                         limits = c(0, 170), 
                         breaks = seq(0, 170, by = 20)) +
      geom_hline(yintercept = c(-log10(pvalue_t)),linewidth = 0.7,color = "black",lty = "dashed") + 
      geom_vline(xintercept = c(-logFC_t, logFC_t),linewidth = 0.7,color = "black",lty = "dashed") + 
      ggtitle(this_tile)+
      mytheme
    p7<-p7+
      geom_point(data = for_label_down,
                 aes(x = logFC, 
                     y = -log10(adj.P.Val)),
                 color = '#2263f5',size = 2, alpha = 1) +
      geom_text_repel(data = for_label_down,
                      aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL),
                      seed = 233,
                      size = 8,
                      fontface='bold',
                      color = '#2263f5',
                      min.segment.length = 0,
                      force = 8, 
                      force_pull = 2, 
                      box.padding = 0.1, 
                      max.overlaps = Inf, 
                      segment.linetype = 3, 
                      segment.color = '#2263f5', 
                      segment.alpha = 0.8,
                      nudge_x = -4 - for_label_down$logFC,
                      direction = "y", 
                      hjust = 1 
      )+
      geom_point(data = for_label_none,
                 aes(x = logFC, 
                     y = -log10(adj.P.Val)),
                 color = 'black',size = 1, alpha = 1) +
      geom_text_repel(data = for_label_none,
                      aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL),
                      seed = 233,
                      size = 5,
                      fontface='italic',
                      color = 'black',
                      min.segment.length = 0,
                      force = 8,
                      force_pull = 2,
                      box.padding = 0.1,
                      max.overlaps = Inf,
                      segment.linetype = 3,
                      segment.color = 'black',
                      segment.alpha = 0.8,
                      nudge_x = 0 - for_label_none$logFC,
                      direction = "y",
                      hjust = 0 
      )+
      geom_point(data = for_label_up,
                 aes(x = logFC, 
                     y = -log10(adj.P.Val)),
                 color = '#925092',size = 2, alpha = 1) +
      geom_text_repel(data = for_label_up,
                      aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL),
                      seed = 233,
                      size = 8,
                      fontface='bold',
                      color = '#925092',
                      min.segment.length = 0,
                      force = 2,
                      force_pull = 2,
                      box.padding = 0.1,
                      max.overlaps = Inf,
                      segment.linetype = 3,
                      segment.color = '#925092',
                      segment.alpha = 0.8,
                      nudge_x = 5 - for_label_up$logFC,
                      direction = "y",
                      hjust = 0 
      )
  }
  
  if(T){ 
    x=deg$logFC 
    names(x)=rownames(deg) 
    cg=c(names(head(sort(x),100)),
         names(tail(sort(x),100)))
    library(pheatmap)
    n=t(scale(t(dat[cg,])))
    n[n>2]=2
    n[n< -2]= -2
    ac=data.frame(group=grouplist$group)
    rownames(ac)=colnames(n)
    palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]
    names(palette) <- names(table(grouplist$group))
    p6 <- pheatmap(n,show_colnames =F,
                   show_rownames = F,
                   cluster_cols = T,
                   clustering_method='average',
                   main = pro,
                   annotation_colors = list(group = palette),
                   annotation_col=ac) 

    ha <- HeatmapAnnotation(

      foo = anno_block(
        gp = gpar(fill = color1),# group
        labels = c(control,case), 
        labels_gp = gpar(col = "white", fontsize = 10)
      )
    )
    grouplist2<-factor(grouplist1,levels = c(control,case))
    p8<-Heatmap(n,
                column_title =paste0(control,'vs',case,'_DEG_top100'),
                color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                show_row_names = F,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 6),
                cluster_columns =T,
                clustering_method_columns = "average",
                cluster_rows = T,
                column_split = grouplist2,
                top_annotation = ha,
                column_names_rot=45,
                column_names_gp = gpar(fontsize = 8),
                row_title = NULL,
                show_row_dend = F,
                heatmap_legend_param = list(
                  at = c(-2, 0, 2),
                  labels = c("low", "zero", "high"),
                  title = 'Expression level',
                  legend_height = unit(4, "cm"),
                  title_position = "leftcenter-rot"
                )
    )
  }
  p_DEG_1 <- p5+as.ggplot(p6)+plot_layout(widths = c(2.5,4))
  p_DEG_2<- p7+as.ggplot(p8)+plot_layout(widths = c(2.5,4))
  ggsave(p_DEG_1,path=file.path(path,pro),
         filename=paste0('DEG_plot_1.png'),
         width=14,height=8)
  ggsave(p_DEG_2,path=file.path(path,pro),
         filename=paste0('DEG_plot_2.png'),
         width=14,height=8)
  write.csv(df, file = file.path(path,pro, paste0('DEG_FC',logFC,'.csv')))
  save(df,file =file.path(path,pro, paste0('DEG_FC',logFC,'.Rdata')))
  return(df)
}