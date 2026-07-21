run_pcahclust<-function(
    path,pro,
    dat, type='count',group, out=NULL,order,
    method1='average',batch=F,
    color,color2){
  library(scatterplot3d)
  expreset <- read.csv(dat,check.names = F)
  # expreset <- read.csv(dat,check.names = F)
  
  grouplist<-read.csv(group,check.names = F)
  grouplist<-grouplist %>% filter(group %in% order)
  grouplist$group<-factor(grouplist$group,levels = order)
  grouplist <- grouplist %>% arrange(group)
  grouplist<-grouplist[!grouplist$sample %in% out,]
  
  # Group
  grouplist1 <- grouplist$group
  names(grouplist1)<-grouplist$sample
  grouplist2<-factor(grouplist1,levels = order)
  group_list <- as.data.frame(grouplist2)
  names(group_list)[1]<-'group'
  
  expreset_tmp<-expreset[match(grouplist$sample, names(expreset))]
  rownames(expreset_tmp)<-expreset[,1]
  
  expreset_tmp2<-t(expreset_tmp)
  constant_columns <- apply(expreset_tmp2, 2, var) == 0
  expreset_tmp2 <- expreset_tmp2[,!constant_columns]
  pca.res <- prcomp(expreset_tmp2, scale. = T, # 标准化
                    center = T # 中心化
  )
  screeplot(pca.res, type = "lines") 
  pcap<- as.data.frame(pca.res$x)
  pcap$group<-grouplist$group
  
  cv <- c()
  pv <- c()
  groups <- unique(grouplist$group)
  
  for (i in 1:length(groups)) {
    gsize <- table(grouplist$group)[groups[i]]
    cv <- c(cv, rep(color[i], each = gsize))
    pv <- c(pv, rep(14+i, each = gsize))
  }
  
  png(file.path(path,paste0(pro,'_pca3d_',type,'(no ',out,').png')),width=8,height = 6,units = 'in',res=300)
  
  scatterplot3d(
    pcap[, 1:3], 
    color = cv,
    pch = pv,    
    lty.hide = 2,
    cex.symbols = 2
  )
  
  legend("topleft",
         order,
         fill=color,
         box.col=NA)
  
  dev.off()
  
  library(vegan)
  library(ape)
  library(phangorn)
  library(dendextend)
  
  expreset2<-expreset_tmp
  
  par(oma = c(0, 0,0, 0))
  par(mar = c(3, 1,2, 5))
  
  dist <- vegdist(t(expreset2),method = 'bray')
  
  
  
  
  hc1 <- hclust(dist,method=method1)
  
  dend <- as.dendrogram(hc1)
  dend <- dend %>%
    set("labels_colors", cv[order.dendrogram(dend)]) %>%
    set("labels_cex", 0.7)
  dend <- dend %>%
    set("branches_k_color", cv[order.dendrogram(dend)])
  labels_cex(dend) <- 0.8
  
  png(file.path(path,paste0(pro,'_hclust_',method1,'_',type,'(no ',out,').png')),
      width=8,height = 6,units = 'in',res=300) 
  # bottom left top right
  plot(dend, type = "rectangle", horiz = TRUE)
  legend("topleft", legend = order,
         fill = color, box.col = NA)
  # average complete mcquitty ward.D
  dev.off()
  
  if (isTRUE(batch)){
    library(limma)
    batch<-grouplist$batch
    design<-model.matrix(~grouplist2)
    expreset_batch<-removeBatchEffect(expreset_tmp, 
                                      batch=batch,design = design)
    expreset_batch2<-t(expreset_batch)
    constant_columns <- apply(expreset_batch2, 2, var) == 0
    expreset_batch2 <- expreset_batch2[,!constant_columns]
    pca.res <- prcomp(expreset_batch2, scale. = T, 
                      center = T
    )
    pcap<- as.data.frame(pca.res$x)
    pcap$group<-grouplist$group
    
    png(file.path(path,paste0(pro,'_batch_pca3d_',type,'(no ',out,').png')),width=12,height = 12,units = 'in',res=300)
    
    scatterplot3d(
      pcap[, 1:3],  
      color = cv,
      pch = pv,  
      lty.hide = 2,
      cex.symbols =3,
      cex.axis=1.5, 
      cex.lab=2	
    )
    
    legend("topleft",
           order,
           fill=color,
           box.col=NA,
           cex = 2
           )
    
    dev.off()
    
    pdf(file.path(path,paste0(pro,'_batch_pca3d_',type,'(no ',out,').pdf')),width=8,height = 8)
    
    scatterplot3d(
      pcap[, 1:3],  
      color = cv,  
      pch = pv,    
      lty.hide = 2,
      cex.symbols =3,
      cex.axis=1.5, 
      cex.lab=2
      
    )
    
    legend("topleft",
           order,
           fill=color,
           box.col=NA,
           cex = 2)
    
    dev.off()
    
    dist <- vegdist(t(expreset_batch),method = 'bray')
    hc1 <- hclust(dist,method=method1)
    
    dend <- as.dendrogram(hc1)
    dend <- dend %>%
      set("labels_colors", cv[order.dendrogram(dend)]) %>%
      set("labels_cex", 0.7)
    dend <- dend %>%
      set("branches_k_color", cv[order.dendrogram(dend)])
    labels_cex(dend) <- 0.8
    
    png(file.path(path,paste0(pro,'_batch_hclust_',method1,'_',type,'(no ',out,').png')),
        width=8,height = 6,units = 'in',res=300) 
    # bottom left top right
    plot(dend, type = "rectangle", horiz = TRUE)
    legend("topleft", legend = order,
           fill = color, box.col = NA)
    # average complete mcquitty ward.D
    dev.off()
  }
}


run_pcahclust2<-function(
    path,pro,
    dat, type='count',group, gtype,out=NULL,order,
    method1='average',batch=F,
    color,color2){
  library(scatterplot3d)
  expreset <- read.csv(dat,check.names = F)

  grouplist<-read.csv(group,check.names = F)
  grouplist<-grouplist %>% filter(group %in% order)
  grouplist$group<-factor(grouplist$group,levels = order)
  grouplist <- grouplist %>% arrange(group)
  grouplist<-grouplist[!grouplist$sample %in% out,]
 
  grouplist1 <- grouplist$group
  
  names(grouplist1)<-grouplist$sample
  grouplist2<-factor(grouplist1,levels = order)
  group_list <- as.data.frame(grouplist2)
  names(group_list)[1]<-'group'
  
  expreset_tmp<-expreset[match(grouplist$sample, names(expreset))]
  rownames(expreset_tmp)<-expreset[,1]
  
  expreset_tmp2<-t(expreset_tmp)
  constant_columns <- apply(expreset_tmp2, 2, var) == 0
  expreset_tmp2 <- expreset_tmp2[,!constant_columns]
  pca.res <- prcomp(expreset_tmp2, scale. = T,
                    center = T
  )
  screeplot(pca.res, type = "lines") 
  pcap<- as.data.frame(pca.res$x)
  pcap$group<-grouplist$group
  
  cv <- c()
  pv <- c()
  groups <- unique(grouplist$group)
  
  for (i in 1:length(groups)) {
    gsize <- table(grouplist$group)[groups[i]]
    cv <- c(cv, rep(color[i], each = gsize))
    pv <- c(pv, rep(14+i, each = gsize))
  }
  
  png(file.path(path,paste0(pro,'_pca3d_',type,'(no ',paste(out,collapse = '&'),').png')),width=8,height = 6,units = 'in',res=300)
  
  scatterplot3d(
    pcap[, 1:3],            
    color = cv,    
    pch = pv,      
    lty.hide = 2,
    cex.symbols = 2
  )
  
  legend("topleft",
         order,
         fill=color,
         box.col=NA)
  
  dev.off()
  
  library(vegan)
  library(ape)
  library(phangorn)
  library(dendextend)
  
  expreset2<-expreset_tmp
  
  par(oma = c(0, 0,0, 0)) 
  par(mar = c(3, 1,2, 5)) 
  
  dist <- vegdist(t(expreset2),method = 'bray')
  
  
  
  
  hc1 <- hclust(dist,method=method1)
  
  dend <- as.dendrogram(hc1)
  dend <- dend %>%
    set("labels_colors", cv[order.dendrogram(dend)]) %>%
    set("labels_cex", 0.7)
  dend <- dend %>%
    set("branches_k_color", cv[order.dendrogram(dend)])
  labels_cex(dend) <- 0.8 # 放大字体
  
  png(file.path(path,paste0(pro,'_hclust_',method1,'_',type,'(no ',paste(out,collapse = '&'),').png')),
      width=8,height = 6,units = 'in',res=300) 
  # bottom left top right
  par(mar = c(5, 4, 4, 8))
  plot(dend, type = "rectangle", horiz = TRUE)
  legend("topleft", legend = order,
         fill = color, box.col = NA)
  # average complete mcquitty ward.D
  dev.off()
  
  if (isTRUE(batch)){
    library(limma)
    batch<-grouplist$batch
    design<-model.matrix(~grouplist2)
    expreset_batch<-removeBatchEffect(expreset_tmp, 
                                      batch=batch,design = design)
    expreset_batch2<-t(expreset_batch)
    constant_columns <- apply(expreset_batch2, 2, var) == 0
    expreset_batch2 <- expreset_batch2[,!constant_columns]
    pca.res <- prcomp(expreset_batch2, scale. = T,
                      center = T
    )
    pcap<- as.data.frame(pca.res$x)
    pcap$group<-grouplist$group
    
    png(file.path(path,paste0(pro,'_batch_pca3d_',type,'(no ',paste(out,collapse = '&'),').png')),width=8,height = 6,units = 'in',res=300)
    
    scatterplot3d(
      pcap[, 1:3],     
      color = cv,   
      pch = pv,   
      lty.hide = 2,
      cex.symbols =3,
      cex.axis=1.5, 
      cex.lab=2	
    )
    
    legend("topleft",
           order,
           fill=color,
           box.col=NA,
           cex = 2
    )
    
    dev.off()
    
    pdf(file.path(path,paste0(pro,'_batch_pca3d_',type,'(no ',paste(out,collapse = '&'),').pdf')),width=8,height = 8)
    
    scatterplot3d(
      pcap[, 1:3],   
      color = cv,    
      pch = pv,         
      lty.hide = 2,
      cex.symbols =3,
      cex.axis=1.5, 
      cex.lab=2
      
    )
    
    legend("topleft",
           order,
           fill=color,
           box.col=NA,
           cex = 2)
    
    dev.off()
    
    dist <- vegdist(t(expreset_batch),method = 'bray')
    hc1 <- hclust(dist,method=method1)
    
    dend <- as.dendrogram(hc1)
    dend <- dend %>%
      set("labels_colors", cv[order.dendrogram(dend)]) %>%
      set("labels_cex", 0.7)
    dend <- dend %>%
      set("branches_k_color", cv[order.dendrogram(dend)])
    labels_cex(dend) <- 0.8 
    
    png(file.path(path,paste0(pro,'_batch_hclust_',method1,'_',type,'(no ',paste(out,collapse = '&'),').png')),
        width=8,height = 6,units = 'in',res=300) 
    # bottom left top right
    par(mar = c(5, 4, 4, 8))
    plot(dend, type = "rectangle", horiz = TRUE)
    legend("topleft", legend = order,
           fill = color, box.col = NA)
    # average complete mcquitty ward.D
    dev.off()
  }
}