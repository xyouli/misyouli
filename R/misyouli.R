#' @title miscellaneous functions I create and use
#'
#' @description
#'
#' @param
#'
#' @return NULL
#'
#' @examples
#'
#' @export
#'
#'

# test assocation between a single signature and a single gene
fisherTest <- function(module_score,CN,one_side = T) {

  topQ <- quantile(module_score,0.75)
  table <- matrix(rep(0,4),nrow = 2,ncol = 2)
  module <- ifelse(module_score>=topQ,"high","low")
  temp <- data.frame(module=module,CN=CN)

  table[1,1] <- length(which(temp$module=="high" & temp$CN=="mut"))
  table[1,2] <- length(which(temp$module=="high" & temp$CN=="wt"))
  table[2,1] <- length(which(temp$module=="low" & temp$CN=="mut"))
  table[2,2] <- length(which(temp$module=="low" & temp$CN=="wt"))

  if(one_side) {
    fish <- fisher.test(table,alternative = 'greater')
  } else {
    fish <- fisher.test(table)
  }
  return(c(fish$p.value,fish$estimate))
}

# test associations between a signature and all genes, both fisher and spearman correlation/lm
moduleCNTest <- function(module_score,CN_score,CN_gain,CN_loss,index,lm = F) {

  score <- unlist(module_score[index,])

  if( lm ) {
    spearman_pos <- c()
    spearman_neg <- c()
    spearman_cor <- c()

    lm_pos <- c()
    lm_neg <- c()
    beta_coeff <- c()
    r.squared <- c()

    gain <- c()
    loss <- c()
    OR <- c()

    for(j in 1:nrow(CN_score)){
      CN <- unname(unlist(CN_score[j,]))

      # spearman rank correlation
      pos <- cor.test(score,CN,method = "spearman",alternative = "greater")
      neg <- cor.test(score,CN,method = "spearman",alternative = "less")
      spearman_pos <- c(spearman_pos,pos$p.value)
      spearman_neg <- c(spearman_neg,neg$p.value)
      spearman_cor <- c(spearman_cor,cor(score,CN,method = "spearman"))

      # linear model
      fit <- lm(score ~ CN + subtype_variable$basal + subtype_variable$her2 + subtype_variable$lumA + subtype_variable$lumB)
      sum <- summary(fit)
      beta <- sum$coefficients[,1][2]
      p <- sum$coefficients[,4][2]
      if(beta>0) {
        lm_pos <- c(lm_pos,p/2)
        lm_neg <- c(lm_neg,1-p/2)
      }
      else if(beta<0) {
        lm_pos <- c(lm_pos,1-p/2)
        lm_neg <- c(lm_neg,p/2)
      }
      beta_coeff <- c(beta_coeff,beta)
      r.squared <- c(r.squared,sum$r.squared)
    }

    # fisher's exact test
    gain <- apply(CN_gain,1,function(x){return(fisherTest(score,unlist(x)))})
    loss <- apply(CN_loss,1,function(x){return(fisherTest(score,unlist(x)))})

    p <- data.frame(spearman_pos = spearman_pos,spearman_neg = spearman_neg, CN_gain = gain[1,], CN_loss = loss[1,],lm_pos = lm_pos,lm_neg = lm_neg,spearman_cor = spearman_cor,lm_beta = beta_coeff,lm_r_squared = r.squared, OR_gain = gain[2,],OR_loss = loss[2,])
    write.table(p,paste("p_value_for_",rownames(module_score)[index],".txt",sep = ""),sep = '\t',col.names = NA)
    return(p)
  }
  else if ( !lm ) {
    spearman_pos <- c()
    spearman_neg <- c()
    spearman_cor <- c()

    # spearman rank correlation
    for(j in 1:nrow(CN_score)){
      CN <- unname(unlist(CN_score[j,]))
      pos <- cor.test(score,CN,method = "spearman",alternative = "greater")
      neg <- cor.test(score,CN,method = "spearman",alternative = "less")
      spearman_pos <- c(spearman_pos,pos$p.value)
      spearman_neg <- c(spearman_neg,neg$p.value)
      spearman_cor <- c(spearman_cor,cor(score,CN,method = "spearman"))
    }

    # fisher's exact test
    gain <- apply(CN_gain,1,function(x){return(fisherTest(score,unlist(x)))})
    loss <- apply(CN_loss,1,function(x){return(fisherTest(score,unlist(x)))})

    p <- data.frame(spearman_pos = spearman_pos,spearman_neg = spearman_neg, CN_gain = gain[1,], CN_loss = loss[1,], spearman_cor = spearman_cor,OR_gain = gain[2,],OR_loss = loss[2,])
    write.table(p,paste("p_value_for_",rownames(module_score)[index],".txt",sep = ""),sep = '\t',col.names = NA)
    return(p)
  }
}

# read in p values for all signatures under the path storing single txt files for each signature
read.p <- function(path,lm = F) {
  p <- list()
  i <- 1
  module <- c()
  if(lm) {
    for (file in list.files(path) ) {
      print(i)
      p[[i]] <- read.table(file,sep = '\t',header = T,row.names = 1,check.names = F)
      p[[i]] <- p[[i]][,1:6] # first 6 columns contain p value, last several contain coefficients
      module_name <- strsplit(file,split = "p_value_for_",fixed = T)[[1]][2]
      module_name <- strsplit(module_name,split = ".txt")[[1]][1]
      module <- c(module,module_name)
      print(module_name)
      i <- i+1
    }
    names(p) <- module
  }
  else if(!lm) {
    for (file in list.files(path) ) {
      print(i)
      p[[i]] <- read.table(file,sep = '\t',header = T,row.names = 1,check.names = F)
      p[[i]] <- p[[i]][,1:4] # first 4 columns contain p value
      module_name <- strsplit(file,split = "p_value_for_",fixed = T)[[1]][2]
      module_name <- strsplit(module_name,split = ".txt")[[1]][1]
      module <- c(module,module_name)
      print(module_name)
      i <- i+1
    }
    names(p) <- module
  }
  return(p)
}

# -log 10 transformation for all p values
logP <- function(p) {
  return(lapply(p,function(x){
    return(-log(x,10))}))
}

# bonferroni correction for all p values
bonferroniByGene <- function(p) {
  return(lapply(p,function(x){
    return(apply(x,2,function(x){return(p.adjust(x,method = "bonferroni"))}))
  }))
}

# bonferroni correction for all p values specifying number of segments to correct
bonferroniBySeg <- function(p,num) {
  return(lapply(p,function(x){
    return(apply(x,2,function(x){return(ifelse(x*num>1,1,x*num))}))
  }))
}

# BH correction for all p values
BH <- function(p) {
  return(lapply(p,function(x){
    return(apply(x,2,function(x){return(p.adjust(x,method = "BH"))}))
  }))
}

# plot p values for all genes for a specific signature
PPlot <- function(p,main,vertical_lines) {
  temp <- get(main,p)
  pos <- temp[,1]
  neg <- temp[,2]
  gain <- temp[,3]
  loss <- temp[,4]
  total_gene <- nrow(temp)
  m1 <- max(pos,gain)
  m2 <- max(neg,loss)
  vertical <- vertical_lines[-c(1,24)]
  text_pos <- vertical[1]/2
  for(i in 2:length(vertical)){
    thispos <- vertical[i] - (vertical[i] - vertical[i-1])/2
    text_pos <- rbind(text_pos,thispos)
  }
  text_pos <- rbind(text_pos,total_gene-(total_gene-vertical[22])/2)

  png(filename = paste(main,'.png',sep = ''),width = 2000,height = 1000)
  par(mfrow = c(2,1),lwd=2,cex.main = 3,cex.lab = 2,cex.axis = 2)

  par(mar=c(0,5,10,3))
  plot(c(1:total_gene),pos,type = "l", xaxt = "n",yaxt = "n",
       ylim = c(0,m1),col = "red",ylab = expression("-log"[10]*'q'),cex = 0.5)
  title(main=main,line = 5)
  lines(c(1:total_gene),gain,xaxt = "n",yaxt = "n",col = "orange")
  axis(2,pos = 0)
  abline(v = vertical,lwd=1)


  par(mar = c(10,5,0,3),lwd=2)
  plot(c(1:total_gene),neg,type = "l", xaxt = "n",yaxt = "n",
       ylim = c(m2,0),col = "darkblue",ylab = expression("-log"[10]*'q'),cex = 0.5,xlab = NA)
  lines(c(1:total_gene),loss,xaxt = "n",yaxt = "n",col = "midnightblue")
  axis(2,pos = 0)
  abline(v = vertical,lwd=1)
  text(text_pos,y = 0.9*m,labels = c(1:22,'x'),cex = 2)
  dev.off()
}

# get significant gene heatmap matrix specifying a significant threshold (alpha)
sigGene <- function(p,alpha,index=c(1:length(p))) {
  sub_p <- p[index]
  hp <- matrix(nrow = length(sub_p),ncol = nrow(p[[1]]))
  threshold <- -log(alpha,10)
  # 1 and 3 column:pos   2 and 4 column:neg
  for ( i in 1:length(sub_p) ) {
    for (j in 1:nrow(p[[1]])) {
      if (sub_p[[i]][j,1]>threshold & sub_p[[i]][j,3]>threshold) {hp[i,j] <- 1}
      else if(sub_p[[i]][j,2]>threshold & sub_p[[i]][j,4]>threshold) {hp[i,j] <- -1}
      else {hp[i,j] <- 0}
    }
  }
  rownames(hp) <- names(sub_p)
  colnames(hp) <- rownames(p[[1]])
  return(hp)
}

heatmap.sigGene <- function(hp,main,labrow = T,q) {
  colBreak <- c(-1,-0.9,0.9,1)
  colramp <- c("midnightblue","white","red4")

  text_pos <- vertical_lines[2]/2
  for(i in 2:(length(vertical_lines)-1)){
    thispos <- vertical_lines[i+1] - (vertical_lines[i+1] - vertical_lines[i])/2
    text_pos <- c(text_pos,thispos)
  }

  labcol <- rep("",vertical_lines[24])
  labcol[round(text_pos)] <- c(1:22,"X")

  hp.wrap <- function(hp) {
    heatmap.2(hp,Rowv = F,Colv = F,dendrogram = "none",scale = "none",trace = "none",
              labRow = labrow,labCol = labcol,density.info = "none",
              breaks = colBreak,col = colramp,key = F,colsep = vertical_lines,sepcolor = "black",
              margins = c(20,60),cexRow = 2,cexCol = 3,srtCol = 360,
              lwid = c(1,4.5),lhei = c(1,8))
  }

  legend.wrap <- function(q) {
    pos <- paste("positive correlation/gain(q<",q,")",sep = "")
    neg <- paste("negative correlation/loss(q<",q,")",sep = "")
    legend("bottom",legend = c(pos,neg),fill = c("red4","midnightblue"),ncol = 2,cex = 3)
  }

  if(labrow) {
    labrow <- rownames(hp)
  }

  height <- 35*nrow(hp)
  png(paste(main,".png",sep = ""),width = 4000,height)
  hp.wrap(hp)
  legend.wrap(q)
  dev.off()
}

immuneNameConversion <- function(x) {
  rownames(x) <- immune_module$short_name[match(rownames(x),immune_module$name)]
  return(x)
}

ampNameConversion <- function(x) {
  rownames(x) <- amplicon_module$short_name[match(rownames(x),amplicon_module$name)]
  return(x)
}

permP <- function(p,perm_p) {
  # perm_p input: a list [[1]]:spearman_pos [[2]]:spearman_neg[[3]]gain[[4]]loss([[5]]lm_pos[[6]]lm_neg)
  p.adj <- function(p,perm_p) {
    # given a vector p and perm_p, adjust p according to perm_p
    return(unlist(lapply(p,function(x){return((sum(perm_p<x)+1)/(length(perm_p)+1))})))
  }
  for (i in 1:length(p)) {
    for(j in 1:length(perm_p)) {
      # ith module jth statistic
      print(paste("adjusting for module:",names(p)[i],"statistic:",names(perm_p)[j]))
      p[[i]][,j] <- p.adj(p[[i]][,j],perm_p[[j]][,i])
    }
  }
  return(p)
}

hp.clust <- function(hp,main,labrow=T,q) {

  colBreak <- c(-1,-0.9,0.9,1)
  colramp <- c("midnightblue","white","red4")
  par(mar = c(2,2,2,2))

  text_pos <- vertical_lines[2]/2
  for(i in 2:(length(vertical_lines)-1)){
    thispos <- vertical_lines[i+1] - (vertical_lines[i+1] - vertical_lines[i])/2
    text_pos <- c(text_pos,thispos)
  }

  labcol <- rep("",vertical_lines[24])
  labcol[round(text_pos)] <- c(1:22,"X")

  hp.wrap <- function(hp) {
    # euclidean distance between module
    dist <- dist(hp,method = "euclidean")
    hc <- hclust(dist,method = 'complete')
    heatmap.2(hp,Rowv = as.dendrogram(hc),Colv = F,dendrogram = "row",scale = "none",trace = "none",
              labRow = labrow,labCol = labcol,density.info = "none",
              breaks = colBreak,col = colramp,key = F,colsep = vertical_lines,sepcolor = "black",
              margins = c(20,60),cexRow = 2.8,cexCol = 3,srtCol = 360,
              lhei = c(1,8))
  }

  legend.wrap <- function(q) {
    pos <- paste("positive correlation/gain(q<",q,")",sep = "")
    neg <- paste("negative correlation/loss(q<",q,")",sep = "")
    legend("bottom",legend = c(pos,neg),fill = c("red4","midnightblue"),ncol = 2,cex = 3)
  }

  if(labrow) {
    labrow <- rownames(hp)
  }

  height <- 35*nrow(hp)
  png(paste(main,".png",sep = ""),width = 4000,height)
  hp.wrap(hp)
  legend.wrap(q)
  dev.off()

}

hp.clust.all <- function(hp,method="complete",main,labrow=F,n=10) {
  # euclidean distance between modules
  dist <- dist(hp,method = "euclidean")
  hc <- hclust(dist,method)
  g <- cutree(hc,k=n)
  write.table(g,paste(main,".txt",sep = ""),sep = '\t')


  colBreak <- c(-1,-0.9,0.9,1)
  colramp <- c("midnightblue","white","red4")
  par(mar = c(2,2,2,2))

  text_pos <- vertical_lines[2]/2
  for(i in 2:(length(vertical_lines)-1)){
    thispos <- vertical_lines[i+1] - (vertical_lines[i+1] - vertical_lines[i])/2
    text_pos <- c(text_pos,thispos)
  }

  labcol <- rep("",vertical_lines[24])
  labcol[round(text_pos)] <- c(1:22,"X")

  if (labrow ) {
    labrow = rownames(hp)
    if(sum(grep("0.01",main))) {
      png(paste(main,".png",sep = ""),width = 4000,height = 3000)
      heatmap.2(hp,Rowv = as.dendrogram(hc),Colv = F,dendrogram = "row",scale = "none",trace = "none",
                labRow = labrow,labCol = labcol,density.info = "none",
                breaks = colBreak,col = colramp,key = F,colsep = vertical_lines,sepcolor = "black",
                margins = c(20,70),cexRow = 3,cexCol = 3,srtCol = 360)
      legend("bottom",legend = c("positive correlation/gain(q<0.01)","negative correlation/loss(q<0.01)"),
             fill = c("red4","midnightblue"),ncol = 2,cex = 3)
      dev.off()
    }
    else if (sum(grep("0.05",main))) {
      png(paste(main,".png",sep = ""),width = 4000,height = 3000)
      heatmap.2(hp,Rowv = as.dendrogram(hc),Colv = F,dendrogram = "row",scale = "none",trace = "none",
                labRow = labrow,labCol = labcol,density.info = "none",
                breaks = colBreak,col = colramp,key = F,colsep = vertical_lines,sepcolor = "black",
                margins = c(20,70),cexRow = 3,cexCol = 3,srtCol = 360)
      legend("bottom",legend = c("positive correlation/gain(q<0.05)","negative correlation/loss(q<0.05)"),
             fill = c("red4","midnightblue"),ncol = 2,cex = 3)
      dev.off()
    }
  }

  else {
    if(sum(grep("0.01",main))) {
      png(paste(main,".png",sep = ""),width = 4000,height = 3000)
      heatmap.2(hp,Rowv = as.dendrogram(hc),Colv = F,dendrogram = "row",scale = "none",trace = "none",
                labRow = F,labCol = labcol,density.info = "none",
                breaks = colBreak,col = colramp,key = F,colsep = vertical_lines,sepcolor = "black",
                margins = c(20,70),cexRow = 3,cexCol = 3,srtCol = 360)
      legend("bottom",legend = c("positive correlation/gain(q<0.01)","negative correlation/loss(q<0.01)"),
             fill = c("red4","midnightblue"),ncol = 2,cex = 3)
      dev.off()
    }
    else if (sum(grep("0.05",main))) {
      png(paste(main,".png",sep = ""),width = 4000,height = 3000)
      heatmap.2(hp,Rowv = as.dendrogram(hc),Colv = F,dendrogram = "row",scale = "none",trace = "none",
                labRow = F,labCol = labcol,density.info = "none",
                breaks = colBreak,col = colramp,key = F,colsep = vertical_lines,sepcolor = "black",
                margins = c(20,70),cexRow = 3,cexCol = 3,srtCol = 360)
      legend("bottom",legend = c("positive correlation/gain(q<0.05)","negative correlation/loss(q<0.05)"),
             fill = c("red4","midnightblue"),ncol = 2,cex = 3)
      dev.off()
    }
  }

}

hp.common <- function(hp1,hp2) {
  gene <- intersect(colnames(hp1),colnames(hp2))
  module <- intersect(rownames(hp1),rownames(hp2))
  hp1 <- hp1[module,gene]
  hp2 <- hp2[module,gene]
  hp <- matrix(0,nrow = length(module),ncol = length(gene))
  for(i in 1:length(module)) {
    for(j in 1:length(gene)) {
      if(hp1[i,j]==1 & hp2[i,j]==1) {
        hp[i,j] <- 1
      } else if (hp1[i,j]==-1 & hp2[i,j]==-1) {
        hp[i,j] <- -1
      }
    }
  }
  colnames(hp) <- gene
  rownames(hp) <- module
  return(hp)
}

hp.keep <- function(hp) {
  rmrow <- c()
  rmcol <- c()
  for (i in 1:nrow(hp)) {
    if (sum(abs(hp[i,]))==0) {
      rmrow <- c(rmrow,i)
    }
  }
  for (j in 1:ncol(hp)) {
    if (sum(abs(hp[,j]))==0) {
      rmcol <- c(rmcol,j)
    }
  }
  return(hp[-rmrow,-rmcol])
}

###########################################################################################
# permutation functions
spearmanPerm <- function(module_score,CN_score,index,nPerm) {
  score <- unlist(module_score[index,])
  for(i in 1:nPerm) {
    temp <- sample(score,length(score),replace = F)
    # for each permutation, record p values of all genes
    spearman_pos <- c()
    spearman_neg <- c()
    for(j in 1:nrow(CN_score)){
      CN <- unname(unlist(CN_score[j,]))
      pos <- cor.test(temp,CN,method = "spearman",alternative = "greater")
      neg <- cor.test(temp,CN,method = "spearman",alternative = "less")
      spearman_pos <- c(spearman_pos,pos$p.value)
      spearman_neg <- c(spearman_neg,neg$p.value)
    }
    # keep the smallest p values
    pos_p <- min(spearman_pos)
    neg_p <- min(spearman_neg)
    write.table(t(c(pos_p,neg_p)),file = paste("spearman_perm_p_value_for_",rownames(module_score)[index],".txt",sep = ""),sep = '\t',col.names = F,row.names = F,append = T)
  }
  return(1)
}

lmPerm <- function(module_score,CN_score,index,nPerm) {
  score <- unlist(module_score[index,])
  for(i in 1:nPerm) {
    temp <- sample(score,length(score),replace = F)
    lm_pos <- c()
    lm_neg <- c()
    for(j in 1:nrow(CN_score)) {
      CN <- unname(unlist(CN_score[j,]))
      fit <- lm(temp ~ CN + subtype_variable$basal + subtype_variable$her2 + subtype_variable$lumA + subtype_variable$lumB)
      sum <- summary(fit)
      beta <- sum$coefficients[,1][2]
      p <- sum$coefficients[,4][2]
      if(beta>0) {
        lm_pos <- c(lm_pos,p/2)
        lm_neg <- c(lm_neg,1-p/2)
      }
      else if(beta<0) {
        lm_pos <- c(lm_pos,1-p/2)
        lm_neg <- c(lm_neg,p/2)
      }
    }
    pos_p <- min(lm_pos)
    neg_p <- min(lm_neg)
    write.table(t(c(pos_p,neg_p)),file = paste("lm_perm_p_value_for_",rownames(module_score)[index],".txt",sep = ""),sep = '\t',col.names = F,row.names = F,append = T)
  }
  return(1)
}

fisherGainPerm <- function(module_score,CN_gain,index,nPerm) {
  score <- unlist(module_score[index,])
  for(i in 1:nPerm) {
    temp <- sample(score,length(score),replace = F)
    gain <- apply(CN_gain,1,function(x){return(fisherTest(temp,unlist(x)))})
    write.table(min(gain[1,]),file = paste("fisher_gain_perm_p_value_for_",rownames(module_score)[index],".txt",sep = ""),sep = '\t',col.names = F,row.names = F,append = T)
  }
  return(1)
}

fisherLossPerm <- function(module_score,CN_loss,index,nPerm) {
  score <- unlist(module_score[index,])
    temp <- sample(score,length(score),replace = F)
    loss <- apply(CN_loss,1,function(x){return(fisherTest(temp,unlist(x)))})
    write.table(min(loss[1,]),file = paste("fisher_loss_perm_p_value_for_",rownames(module_score)[index],".txt",sep = ""),sep = '\t',col.names = F,row.names = F,append = T)
  }
  return(1)
}

perm.wrap <- function(module_score,CN_score,CN_gain,CN_loss,index,nPerm,home) {
  score <- unlist(module_score[index,])
  module_name <- rownames(module_score)[index]

  for(i in 1:nPerm) {
    print(paste("permutation:",i))
    temp <- sample(score,length(score),replace = F)

    #spearman
    spearman_pos <- c()
    spearman_neg <- c()
    for(j in 1:nrow(CN_score)){
      CN <- unname(unlist(CN_score[j,]))
      pos <- cor.test(temp,CN,method = "spearman",alternative = "greater")
      neg <- cor.test(temp,CN,method = "spearman",alternative = "less")
      spearman_pos <- c(spearman_pos,pos$p.value)
      spearman_neg <- c(spearman_neg,neg$p.value)
    }
    # keep the smallest p values
    pos_p <- min(spearman_pos)
    neg_p <- min(spearman_neg)
    setwd(paste(home,'spearman',sep = '/'))
    write.table(t(c(pos_p,neg_p)),file = paste("spearman_perm_p_value_for_",module_name,".txt",sep = ""),sep = '\t',col.names = F,row.names = F,append = T)

    # lm
    lm_pos <- c()
    lm_neg <- c()
    for(j in 1:nrow(CN_score)) {
      CN <- unname(unlist(CN_score[j,]))
      fit <- lm(temp ~ CN + subtype_variable$basal + subtype_variable$her2 + subtype_variable$lumA + subtype_variable$lumB)
      sum <- summary(fit)
      beta <- sum$coefficients[,1][2]
      p <- sum$coefficients[,4][2]
      if(beta>0) {
        lm_pos <- c(lm_pos,p/2)
        lm_neg <- c(lm_neg,1-p/2)
      }
      else if(beta<0) {
        lm_pos <- c(lm_pos,1-p/2)
        lm_neg <- c(lm_neg,p/2)
      }
    }
    pos_p <- min(lm_pos)
    neg_p <- min(lm_neg)
    setwd(paste(home,'lm',sep = '/'))
    write.table(t(c(pos_p,neg_p)),file = paste("lm_perm_p_value_for_",module_name,".txt",sep = ""),sep = '\t',col.names = F,row.names = F,append = T)

    # fisher
    gain <- apply(CN_gain,1,function(x){return(fisherTest(temp,unlist(x)))})
    loss <- apply(CN_loss,1,function(x){return(fisherTest(temp,unlist(x)))})
    setwd(paste(home,'fisher',sep = '/'))
    write.table(t(c(min(gain[1,]),min(loss[1,]))),file = paste("fisher_perm_p_value_for_",module_name,".txt",sep = ""),sep = '\t',col.names = F,row.names = F,append = T)

  }
}


# convert gene symbol to gene id for matrix x
symbol_2_id <- function(x,symbol_id) {
  gene <- rownames(x)
  # gene <- unlist(lapply(gene_name,function(x){strsplit(x,split = '_',fixed = T)[[1]][1]}))
  index <- match(gene,symbol_id$hgnc_symbol)
  sub_gene <- gene[!(is.na(index))]
  sub_x <- x[!is.na(index),]
  new_index <- match(sub_gene,symbol_id$hgnc_symbol)
  gene_id <- symbol_id$entrezgene[new_index]
  rownames(sub_x) <- gene_id
  return(sub_x)
}

# calculate segment score from gene-level CN score
calc_segments<-function(x, gmtFile, method="mean", scale=F, gsaObj=NA){

  geneset.obj<- GSA.read.gmt(gmtFile)
  genenames<-row.names(x)
  np=length(geneset.obj$genesets)

  if(scale){xs=t(scale(t(x),center=T,scale=T))}else{xs<-x}

  if(method!="gsa"){
    val=matrix(0,nrow=np,ncol=ncol(x))
    for(i in 1:np){
      gene.set=which(genenames %in% geneset.obj$genesets[[i]])
      gene.set=gene.set[!is.na(gene.set)]
      if(length(gene.set)>1){
        if(method=="mean"){
          val[i,]=colSums(xs[gene.set,,drop=F])/length(gene.set)
        }
        if(method=="median"){
          val[i,]=t(apply(xs[gene.set,,drop=F],2,median))
        }
        if(method=="pca"){
          y<-prcomp(as.matrix(xs[gene.set,,drop=F]))
          val[i,]=y$rotation[,1]
        }
      } else if (length(gene.set) == 1) {
        val[i,] <- unlist(xs[gene.set,])
      }
    }
  }else{
    val<-GSA.make.features(gsaObj,xs,geneset.obj$genesets,genenames)
  }

  dimnames(val)<-list(geneset.obj$geneset.names,dimnames(x)[[2]])
  return(val)
}

# calculate signature
calc_modules<-function(x, geneset.obj, method="mean", scale=F, gsaObj=NA){

  # geneset.obj<- GSA.read.gmt(gmtFile)
  genenames<-row.names(x)
  np=length(geneset.obj$genesets)

  if(scale){xs=t(scale(t(x),center=T,scale=T))}else{xs<-x}

  if(method!="gsa"){
    val=matrix(NA,nrow=np,ncol=ncol(x))
    for(i in 1:np){
      gene.set=which(genenames %in% geneset.obj$genesets[[i]])
      gene.set=gene.set[!is.na(gene.set)]
      if(length(gene.set)>1){
        if(method=="mean"){
          val[i,]=colSums(xs[gene.set,,drop=F])/length(gene.set)
        }
        if(method=="median"){
          val[i,]=t(apply(xs[gene.set,,drop=F],2,median))
        }
        if(method=="pca"){
          y<-prcomp(as.matrix(xs[gene.set,,drop=F]))
          val[i,]=y$rotation[,1]
        }
      } else if (length(gene.set) == 1) {
        val[i,] <- unlist(xs[gene.set,])
      }
    }
  }else{
    val<-GSA.make.features(gsaObj,xs,geneset.obj$genesets,genenames)
  }
  dimnames(val)<-list(geneset.obj$geneset.names,dimnames(x)[[2]])
  return(val)
}

medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}

standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}

overlapSets<-function(x,y){
  x<-x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
  y<-y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]

  x<-x[sort.list(row.names(x)),]
  y<-y[sort.list(row.names(y)),]

  return(list(x=x,y=y))
}

assignDiffScore.dwd<-function(x,y){
  both<-overlapSets(x,y) 	# get the overlap of genes
  both$x<- apply(both$x,2,function(x){sign(x)*sqrt(x^2/sum(x^2))}) 	# set the distance to the origin to 1
  both$y<- apply(both$y,2,function(x){sign(x)*sqrt(x^2/sum(x^2))})	# set the distance to the origin to 1
  msproj<- apply(both$y,2,function(x,y){x%*%y},both$x[,1])	# project the samples on the MS-pL axis
  mlproj<- apply(both$y,2,function(x,y){x%*%y},both$x[,2])	# project the samples on the pL-mL axis
  diffScore<- mlproj - msproj
  return( diffScore )	# return the point on the differentiation axis
}

# calculate GHI_RS signature score
GHI_RS <- function(edata) {
  getGene <- function(edata,gene){
    return(unlist(edata[rownames(edata)==gene,]))
  }
  gene <- c(2597,2990,60,7037,6175,2886,2064,596,5241,57758,2099,6790,4605,891,332,4288,4320,1515,968,2944,573)
  for(i in 1:length(gene)){
    assign(paste('X',gene[i],sep = ''),getGene(edata,gene[i]))
  }
  # normalized according to reference gene
  reference.Avg <- (X2597+X2990+X60+X7037+X6175)/5
  refNorm <- function(x,ref) {
    x <- x-ref
    x <- x-min(x)
    x <- x*15/max(x)
    return(x)
  }
  for(i in 1:length(gene)){
    assign(paste('X',gene[i],sep = ''),refNorm(get(paste('X',gene[i],sep = '')),reference.Avg))
  }

  GRB7_Group <- 0.9*X2886 + 0.1*X2064
  for(i in 1:ncol(edata))
  {
    if(GRB7_Group[i]<8)
    {GRB7_Group[i]<-8}
  }

  ER_Group<- (X596+1.2*X5241+X57758+0.8*X2099)/4

  Prolif_Group<- (X6790+X4605+X891+X332+X4288)/5

  for(i in 1:ncol(edata))
  {
    if(Prolif_Group[i]<6.5)
    {Prolif_Group[i]<-6.5}
  }

  Invasion_Group <- (X4320+X1515)/2

  CD68 <- X968
  GSTM1 <- X2944
  BAG1 <- X573

  RSU = 0.47*GRB7_Group - 0.34*ER_Group + 1.04*Prolif_Group + 0.10*Invasion_Group + 0.05*CD68 - 0.08*GSTM1 - 0.07*BAG1
  RS = 20*(RSU - 6.7)
  # GHI_RS <- 20*(RSU - 6.7)
  # GHI_RS_3Group <- 20*(RSU - 6.7)
  #
  # for(i in 1:ncol(edata))
  # {
  #   if(GHI_RS[i]<0){GHI_RS[i]<-0}
  #   if(GHI_RS[i]>100){GHI_RS[i]<-100}
  #
  #   if(GHI_RS[i]<18){GHI_RS_3Group[i]<-"Low"}
  #   if(GHI_RS[i]>=18 && GHI_RS[i]< 31){GHI_RS_3Group[i]<- "Intermediate"}
  #   if(GHI_RS[i]>= 31){GHI_RS_3Group[i]<- "High"}
  # }

  return(RS)
}

# wrap for process PanCan gene expression data
exp_wrap <- function(edata) {
  keep <- c('29126','1493','5133')
  edata70 <- edata[rowSums(edata<=2)<(0.3*ncol(edata)),]
  for(i in 1:3) {
    if(!(keep[i] %in% rownames(edata70))) {
      edata70 <- rbind(edata70,edata[keep[i],])
    }
  }
  edata70[edata70<=2] <- 0  # missing data marked with 0
  edata70log2 <- log(edata70,base = 2)
  edata70log2[edata70log2=="-Inf"] <- 0

  exp <- medianCtr(edata70log2)
  exp <- standardize(exp)
  return(exp)
}


# wrap for Elastic Net model, for multiple signatures and collect beta/AUC to result directory
caret_wrap <- function(trainX,trainY,testX,testY,bi,type,working_dir,result_dir) {
  if(!bi) {
    setwd(working_dir)

    # set cross validation resampling method
    train_control <- trainControl(method = 'LGOCV',number = 200,classProbs = F)

    alpha <- seq(0.1,0.9,by=0.1)
    lambda <- list()
    for(i in 1:9) {
      init <- glmnet(trainX,trainY,alpha = alpha[i])
      lambda[[i]] <- init$lambda
    }
    lambda_min <- min(unlist(lapply(lambda,min)))
    lambda_max <- max(unlist(lapply(lambda,max)))

    tune_grid = expand.grid(alpha = seq(0.1,0.9,by=0.1), lambda = seq(lambda_min,lambda_max,length.out = 100))

    # train model
    glmnet_obj <- train(trainX, trainY, method = "glmnet", metric = "RMSE",
                        trControl = train_control,tuneGrid = tune_grid)

    save(glmnet_obj,file = paste(type,'glmnet_obj.rda',sep = '_'))

    # model performance
    pred_train <- predict(glmnet_obj,newdata = trainX)
    pred_test <- predict(glmnet_obj,newdata = testX)
    cor_train <- cor(pred_train,trainY)
    cor_test <- cor(pred_test,testY)

    # collect results for different modules
    setwd(result_dir)

    # collect beta coefficients
    beta <- as.matrix(coef(glmnet_obj$finalModel,glmnet_obj$bestTune$lambda))
    beta <- t(beta)
    rownames(beta) <- this_module
    write.table(beta,paste(type,'beta_coefficients.txt',sep = '_'),sep = '\t',col.names = F,append = T)

    stat <- t(c(cor_train,cor_test))
    write.table(stat,paste(type,'pred_cor_report.txt',sep = '_'),sep = '\t',col.names = F,row.names = this_module,append = T)

  }  else {
    setwd(working_dir)

    trainY2 <- ifelse(trainY==1,'high','low')

    # set cross validation resampling method
    train_control <- trainControl(method = 'LGOCV',number = 200,classProbs = T)

    alpha <- seq(0.1,0.9,by=0.1)
    lambda <- list()
    for(i in 1:9) {
      init <- glmnet(trainX,trainY,alpha = alpha[i],family = 'binomial')
      lambda[[i]] <- init$lambda
    }
    lambda_min <- min(unlist(lapply(lambda,min)))
    lambda_max <- max(unlist(lapply(lambda,max)))

    tune_grid = expand.grid(alpha = seq(0.1,0.9,by=0.1), lambda = seq(lambda_min,lambda_max,length.out = 100))

    # train model
    glmnet_obj <- train(trainX, trainY2, method = "glmnet", metric = "Accuracy",
                        trControl = train_control,tuneGrid = tune_grid)

    save(glmnet_obj,file = paste(type,'glmnet_obj.rda',sep = '_'))

    # roc curve
    pred_train <- predict(glmnet_obj,newdata = trainX,type = 'prob')
    pred_test <- predict(glmnet_obj,newdata = testX,type = 'prob')

    roc_train<-roc.area(trainY, pred_train[,1])
    roc_test<-roc.area(testY, pred_test[,1])

    # collect results for different modules
    setwd(result_dir)
    # collect beta coefficients
    beta <- as.matrix(coef(glmnet_obj$finalModel,glmnet_obj$bestTune$lambda))
    beta <- t(beta)
    rownames(beta) <- this_module
    write.table(beta,paste(type,'beta_coefficients.txt',sep = '_'),sep = '\t',col.names = F,append = T)

    stat <- t(c(roc_train$A,roc_test$A))
    write.table(stat,paste(type,'AUC_report.txt',sep = '_'),sep = '\t',col.names = F,row.names = this_module,append = T)
  }
}

# wrap for Elastic Net model, classification, suitable for single signature
caret_wrap_ss <- function(trainX,trainY,testX,testY) {
  trainY2 <- ifelse(trainY==1,'high','low')

  # set cross validation resampling method
  train_control <- trainControl(method = 'LGOCV',number = 200,classProbs = T)

  alpha <- seq(0.1,0.9,by=0.1)
  lambda <- list()
  for(i in 1:9) {
    init <- glmnet(trainX,trainY,alpha = alpha[i],family = 'binomial')
    lambda[[i]] <- init$lambda
  }
  lambda_min <- min(unlist(lapply(lambda,min)))
  lambda_max <- max(unlist(lapply(lambda,max)))

  tune_grid = expand.grid(alpha = seq(0.1,0.9,by=0.1), lambda = seq(lambda_min,lambda_max,length.out = 100))

  # train model
  glmnet_obj <- train(trainX, trainY2, method = "glmnet", metric = "Accuracy",
                      trControl = train_control,tuneGrid = tune_grid)

  save(glmnet_obj,file = 'glmnet_obj.rda')

  # collect beta coefficients
  beta <- as.matrix(coef(glmnet_obj$finalModel,glmnet_obj$bestTune$lambda))
  write.table(beta,'beta_coefficients.txt',sep = '\t',col.names = NA)

  # roc curve
  pred_train <- predict(glmnet_obj,newdata = trainX,type = 'prob')
  pred_test <- predict(glmnet_obj,newdata = testX,type = 'prob')
  pred_train <- prediction(pred_train$high,labels = trainY)
  pred_test <- prediction(pred_test$high,labels = testY)
  perf_train <- performance(pred_train,measure = 'tpr',x.measure = 'fpr')
  perf_test <- performance(pred_test,measure = 'tpr',x.measure = 'fpr')
  auc_train <- signif(performance(pred_train,measure = 'auc')@y.values[[1]][1],2)
  auc_test <- signif(performance(pred_test,measure = 'auc')@y.values[[1]][1],2)
  plot_ROC(perf_train,perf_test,auc_train,auc_test,'ROC')

  stat <- t(c(auc_train,auc_test))
  write.table(stat,'AUC_report.txt',sep = '\t',col.names = F,row.names = F)
}

# plot ROC with two ROCs
plot_ROC <- function(perf1,perf2,a1,a2,main) {
  tiff(paste0(main,'.tiff'),width = 1.5,height = 1.5,units = 'in',res = 300)
  par(mai = c(0.2,0.2,0.05,0.05),cex.axis = 0.3)
  plot(perf1,col = 'red',lwd = 1.2,xlab = "",ylab = "",box.lwd=0.8,
       xaxis.xaxt = 'n',yaxis.yaxt = 'n')
  plot(perf2,col = 'darkblue',lwd = 1.2,add = T)
  abline(a=0,b=1,lwd = 0.8)
  axis(1,tck=(-0.02),lwd = 0.8)
  axis(2,tck=(-0.02),lwd = 0.8)
  mtext(side = 1,text = seq(0,1,by=0.2),at = seq(0,1,by=0.2),cex = 0.5,line = (-0.2))
  mtext(side = 2,text = seq(0,1,by=0.2),at = seq(0,1,by=0.2),cex = 0.5,line = 0.1)

  legend('bottomright',legend = c(paste0('AUC = ',a1),paste0('AUC = ',a2)), lty = c(1,1),lwd = c(1,1) ,
         col = c('red','darkblue'),cex = 0.6,bty = 'n')
  dev.off()
}

x_y <- function(beta,segment_anno) {
  x <- c()
  y <- c()
  for(i in 1:nrow(segment_anno)) {
    this_seg <- rownames(segment_anno)[i]
    this_start <- segment_anno[i,2]
    this_end <- segment_anno[i,3]
    this_x <- seq(this_start,this_end,by = 1)
    this_y <- rep(beta[this_seg],length(this_x))
    x <- c(x,this_x)
    y <- c(y,this_y)
  }
  return(data.frame(x=x,y=y))
}

# plot model feature for single signature
plot_seg_ss <- function(beta_seg,main) {
  total_gene <- 24776
  vertical <- vertical_lines[-c(1,24)]
  text_pos <- vertical[1]/2
  for(i in 2:length(vertical)){
    thispos <- vertical[i] - (vertical[i] - vertical[i-1])/2
    text_pos <- rbind(text_pos,thispos)
  }
  text_pos <- rbind(text_pos,total_gene-(total_gene-vertical[22])/2)

  index <- match(main,colnames(beta_seg))
  beta <- beta_seg[,index]
  names(beta) <- rownames(beta_seg)

  min_y <- min(beta)
  max_y <- max(beta)
  # beta whole arm
  index <- grep('wholearm',names(beta))
  beta_wholearm <- beta[index]
  beta <- beta[-index]

  pos_beta <- beta[beta>0]
  neg_beta <- beta[beta<0]

  pos_seg <- names(pos_beta)
  neg_seg <- names(neg_beta)

  pos_anno <- segment_anno[pos_seg,]
  neg_anno <- segment_anno[neg_seg,]

  # pos regions excluding whole arms
  pos_coor <- x_y(beta,pos_anno)
  neg_coor <- x_y(beta,neg_anno)

  x <- c(pos_coor$x,neg_coor$x)
  y <- c(pos_coor$y,neg_coor$y)
  color <- c(rep('red',nrow(pos_coor)),rep('darkblue',nrow(neg_coor)))

  # whole arm beta
  beta_pos_wholearm <- beta_wholearm[beta_wholearm>0]
  beta_neg_wholearm <- beta_wholearm[beta_wholearm<0]

  if( length(beta_pos_wholearm) == 0) {
    pos_wholearm_coor <- data.frame(x=0,y=0)
  } else {
    pos_wholearm_anno <- segment_anno[names(beta_pos_wholearm),]
    pos_wholearm_coor <- x_y(beta_wholearm,pos_wholearm_anno)
  }

  if(length(beta_neg_wholearm)==0) {
    neg_wholearm_coor <- data.frame(x=0,y=0)
  } else {
    neg_wholearm_anno <- segment_anno[names(beta_neg_wholearm),]
    neg_wholearm_coor <- x_y(beta_wholearm,neg_wholearm_anno)
  }
  x_wholearm <- c(pos_wholearm_coor$x,neg_wholearm_coor$x)
  y_wholearm <- c(pos_wholearm_coor$y,neg_wholearm_coor$y)
  color_wholearm <- c(rep('pink',nrow(pos_wholearm_coor)),rep('lightblue',nrow(neg_wholearm_coor)))
  par(cex.axis = 2)
  png(filename = paste(main,'.png',sep = ''),width = 28,height = 5,res = 72,units = 'in')
  par(cex.axis = 2,cex.lab = 2.5,mai=c(0.6,1.5,0.6,0.5))
  plot(x_wholearm,y_wholearm,type = 'h',col = color_wholearm,xlim=c(1,24776),ylim=c(min_y,max_y),lwd = 1,xaxs="i",xaxt = 'n',ylab = '',xlab = "")
  abline(v=vertical_lines,lwd = 1)
  abline(h = 0,lwd = 1)
  lines(x,y,xaxt = "n",yaxt = "n",col = color,type = 'h')
  mtext(c(1:20,"",22,'x'),side = 1,at = text_pos,line = 1.5,cex = 2.5)
  mtext('weight',side = 2,at = 0,line = 4,cex = 2.5)
  dev.off()
}

# deal with NAs in CNGene data matrix
# extreme method, a gene is assigned the greatest amplification or the least deletion

CN_Seg2Gene <- function(data = 'seg.txt', anno = 'hg18.txt', index) {
  data <- read.table(data,sep = '\t',header = T,stringsAsFactors = F)
  anno <- read.table(anno,sep = '\t',header = T,stringsAsFactors = F)
  sample_list <- unique(data$sample)

  this_sample <- sample_list[,1][index]
  print(paste(index,this_sample))

  seg <- filter(data,Sample==this_sample)
  seg[,3] <- as.integer(seg[,3])
  seg[,4] <- as.integer(seg[,4])
  CN_score <- rep(0,nrow(anno))

  # process data segment by segment, should be faster than gene by gene
  for( i in 1:nrow(seg)) {
    segment_chro <- seg[i,2]
    segment_start <- seg[i,3]
    segment_end <- seg[i,4]
    segment_score <- seg[i,6]

    # find genes completely fall in this segment
    this_anno <- filter(anno,Chromosome==segment_chro) %>% filter(between(Gene_Start,segment_start,segment_end) & between(Gene_End,segment_start,segment_end))
    if(nrow(this_anno) > 0) {
      gene_index <- match(this_anno[,5],anno[,5])
      # assign scores to these genes
      CN_score[gene_index] <- segment_score
    }
  }

  # deal with genes fall in multiple segments
  # process gene by gene
  NA_index <- which(CN_score==0)
  for ( i in NA_index) {
    gene_name <- anno[i,5]
    gene_chro <- anno[i,2]
    gene_start <- anno[i,3]
    gene_end <- anno[i,4]

    # find segments this gene falls in
    this_seg <- filter(seg,Chromosome==gene_chro) %>% filter(between(Segment_Start,gene_start,gene_end) | between(Segment_End,gene_start,gene_end))
    # find the greatest amplification or the least deletion
    if (nrow(this_seg) > 0) {
      abs_score <- abs(this_seg[,6])
      segment_index <- which.max(abs_score)
      if (length(segment_index)==1) {
        CN_score[i] <- this_seg[segment_index,6]
      }
    }
  }

  CN_score <- as.data.frame(t(CN_score))
  rownames(CN_score) <- this_sample
  write.table(CN_score,paste(this_sample,'.txt',sep = ''),sep = '\t',col.names = F,row.names = T)
}

opt_cut <- function(perf) {
  df <- data.frame(cut = perf@alpha.values[[1]],spec = perf@x.values[[1]],sens = perf@y.values[[1]])
  J = df$sens + df$spec
  index <- which.max(J)
  return(df$cut[index])
}

uq_norm<- function(x,y=NA){
  uqs<- apply(x,2,function(x){ quantile(x[x>0 & !is.na(x)],0.75)})
  if(is.na(y)){
    y<- median(uqs)
  }
  x.norm <- t(apply(x,1,function(x,y){x*y},y/uqs))
  dimnames(x.norm)<- dimnames(x)
  return(x.norm)
}

calc_signature_wrap <- function(edata,geneset.obj,diff_centroid) {
  module_score <- calc_modules(edata,geneset.obj,method = "median")

  NAmodule <- rownames(module_score[is.na(module_score[,1]),])
  if(!is.null(NAmodule)) {
    print(NAmodule)
    index <- match(NAmodule,rownames(module_score))
    module_score <- module_score[-index,]
  }

  # add PD1,PDL1,CTLA4
  PDL1 <- unlist(edata[rownames(edata)==29126,])
  module_score <- rbind(module_score,PDL1)
  CTLA4 <- unlist(edata[rownames(edata)==1493,])
  module_score <- rbind(module_score,CTLA4)
  PD1 <- unlist(edata[rownames(edata)==5133,])
  module_score <- rbind(module_score,PD1)

  # CD103_Ratio
  CD103_pos <- module_score[rownames(module_score)=="CD103_Positive_Median_Cancer.Cell.2014_PMID.25446897"]
  CD103_neg <- module_score[rownames(module_score)=="CD103_Negative_Median_Cancer.Cell.2014_PMID.25446897"]
  CD103_ratio <- CD103_pos - CD103_neg # log2 scale division
  module_score <- rbind(module_score,CD103_ratio)
  rownames(module_score)[nrow(module_score)] <- "CD103_Ratio_Cancer.Cell.2014_PMID.25446897"

  # Differentiation score
  diff_score <- assignDiffScore.dwd(diff_centroid,edata)
  module_score <- rbind(module_score,diff_score)
  rownames(module_score)[nrow(module_score)] <- 'UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035'

  # GHI_RS_Model
  GHI_RS_Model_NJEM.2004_PMID.15591335 <- GHI_RS(edata)
  module_score <- rbind(module_score,GHI_RS_Model_NJEM.2004_PMID.15591335)


  # as.numeric
  dim_name <- dimnames(module_score)
  module_score <- apply(module_score,2,as.numeric)
  dimnames(module_score) <- dim_name
  return(module_score)
}

plot_surv <- function(t,delta,group,col,xlab = 't',ylab = 'survival',main) {
  df <- length(table(group))-1

  legend.group <- function(diff) {
    N <- df+1
    name <- names(diff$n)
    name <- unlist(lapply(name,function(x){return(strsplit(x,split = '=')[[1]][2])}))
    text <- paste(name[1],paste(diff$obs[1],diff$n[1],sep = '/'))
    for(i in 2:N) {
      text <- c(text,paste(name[i],paste(diff$obs[i],diff$n[i],sep = '/')))
    }
    return(text)
  }

  fit <- survfit(Surv(t,delta) ~ group,conf.type='none')
  diff <- survdiff(Surv(t,delta) ~ group)
  p <- signif(1-pchisq(diff$chisq,df=df),digits=3)
  plot(fit,xlab = xlab,ylab = ylab,mark.time = T,col = col, main = main)
  legend('bottomright',paste('log rank p =',p),text.col = ifelse(p<0.05,'red','black'))
  legend('bottomleft',legend = legend.group(diff),lty = c(1,1), col = col)
}

get_surv_stat <- function(t,delta,module_score) {
  sig <- colnames(module_score)
  module_score_bi <- apply(module_score,2,function(x){ifelse(x>=quantile(x,0.67),1,0)})
  colnames(module_score_bi) <- paste0('X',1:ncol(module_score))
  data <- data.frame(t=t,delta=delta)
  data <- cbind(data,module_score_bi)
  covariates <- colnames(module_score_bi)
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(t, delta)~', x)))

  univ_models <- lapply( univ_formulas, function(x){coxph(x,data)})
  # Extract data
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           p_value<-signif(x$wald["pvalue"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR_lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR_upper <- signif(x$conf.int[,"upper .95"],2)
                           res<-c(beta, HR, HR_lower, HR_upper,p_value)
                           names(res)<-c("beta", "HR", "HR_lower", "HR_upper", "p_value")
                           return(res)
                         })
  res <-  as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
  res <- res %>% mutate(Signature = sig)
  return(res)
}

# elastic net prediction to survival
get_surv_stat_EN <- function(t,delta,module_score,CN_score) {
  sig <- colnames(module_score)
  res <- lapply(sig,
                function(s){
                  load(paste0("E:/longleaf/elastic_net_glmnet_obj/whole_genome/all_module/glmnet_all_module/TCGA_GISTIC2_segment/",s,'/ca_co_glmnet_obj.rda'))
                  pred <- predict(glmnet_obj,newdata = CN_score,type = 'prob')
                  prob <- pred$high
                  prob_bi <- ifelse(prob >= quantile(prob,0.67),1,0)
                  cox <- coxph(Surv(t,delta) ~ prob_bi)
                  x <- summary(cox)
                  p_value<-signif(x$wald["pvalue"], digits=2)
                  beta<-signif(x$coef[1], digits=2);#coeficient beta
                  HR <-signif(x$coef[2], digits=2);#exp(beta)
                  HR_lower <- signif(x$conf.int[,"lower .95"], 2)
                  HR_upper <- signif(x$conf.int[,"upper .95"],2)
                  res<-c(beta, HR, HR_lower, HR_upper,p_value)
                  names(res)<-c("beta", "HR", "HR_lower", "HR_upper", "p_value")
                  return(res)
                })
  res <- as.data.frame(t(as.data.frame(res)))
  res <- res %>% mutate(Signature = sig)
  return(res)
}

get_pCR_stat <- function(pCR,module_score) {
  sig <- colnames(module_score)
  module_score_bi <- apply(module_score,2,function(x){ifelse(x>=quantile(x,0.67),1,0)})

  get_table <- function(pCR,score_bi) {
    table <- matrix(0,nrow = 2,ncol = 2)
    table[1,1] <- length(which(score_bi == 1 & pCR == 1))
    table[1,2] <- length(which(score_bi == 1 & pCR == 0))
    table[2,1] <- length(which(score_bi == 0 & pCR == 1))
    table[2,2] <- length(which(score_bi == 0 & pCR == 0))
    return(table)
  }

  univ_table <- lapply(sig, function(s){get_table(pCR,module_score_bi[,s])})
  univ_models <- lapply(univ_table, function(x){fisher.test(x)})
  # Extract data
  univ_results <- lapply(univ_models,
                         function(x){
                           p_value <- x$p.value
                           OR <- x$estimate
                           res<-c(OR,p_value)
                           names(res)<-c("OR","p_value")
                           return(res)
                         })
  res <-  as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
  res <- res %>% mutate(Signature = sig)
  return(res)
}

get_pCR_stat_EN <- function(pCR,module_score, CN_score) {
  sig <- colnames(module_score)

  get_table <- function(pCR,score_bi) {
    table <- matrix(0,nrow = 2,ncol = 2)
    table[1,1] <- length(which(score_bi == 1 & pCR == 1))
    table[1,2] <- length(which(score_bi == 1 & pCR == 0))
    table[2,1] <- length(which(score_bi == 0 & pCR == 1))
    table[2,2] <- length(which(score_bi == 0 & pCR == 0))
    return(table)
  }

  res <- lapply(sig,
                function(s){
                  load(paste0("E:/longleaf/elastic_net_glmnet_obj/whole_genome/all_module/glmnet_all_module/TCGA_GISTIC2_segment/",s,'/ca_co_glmnet_obj.rda'))
                  pred <- predict(glmnet_obj,newdata = CN_score,type = 'prob')
                  prob <- pred$high
                  prob_bi <- ifelse(prob >= quantile(prob,0.67),1,0)
                  table <- get_table(pCR,prob_bi)
                  x <- fisher.test(table)
                  p_value <- x$p.value
                  OR <- x$estimate
                  res<-c(OR,p_value)
                  names(res)<-c("OR","p_value")
                  return(res)
                })

  res <-  as.data.frame(t(as.data.frame(res)))
  res <- res %>% mutate(Signature = sig)
  return(res)
}

get_pCR_AUC_EN <- function(pCR,module_score, CN_score) {
  sig <- colnames(module_score)

  res <- lapply(sig,
                function(s){
                  load(paste0("E:/longleaf/elastic_net_glmnet_obj/whole_genome/all_module/glmnet_all_module/TCGA_GISTIC2_segment/",s,'/ca_co_glmnet_obj.rda'))
                  pred <- predict(glmnet_obj,newdata = CN_score,type = 'prob')
                  pred <- prediction(pred$high,labels = pCR)
                  perf <- performance(pred,measure = 'tpr',x.measure = 'fpr')
                  auc <- signif(performance(pred,measure = 'auc')@y.values[[1]][1],2)
                  return(auc)
                })

  res <-  as.data.frame(t(as.data.frame(res)))
  colnames(res) <- 'AUC'
  res <- res %>% mutate(Signature = sig)
  return(res)
}

Maf2Matrix <- function(x) {
  sample <- unique(x$Tumor_Sample_Barcode)
  gene <- unique(x$Hugo_Symbol)
  result <- matrix(0,nrow = length(gene),ncol = length(sample))
  for ( i in 1:length(sample) ) {
    print(paste('Sample:',i))
    this_sample <- sample[i]
    this_MAF <- filter(x,Tumor_Sample_Barcode == this_sample)
    this_gene <- unique(this_MAF$Hugo_Symbol)
    index <- match(this_gene,gene)
    result[index,i] <- 1
  }
  rownames(result) <- gene
  colnames(result) <- sample
  return(result)
}
