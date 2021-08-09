#Routine R functions used in RNA-seq analyses

inputStringTieRPKMfiles <- function(sampleName, RPKMDataPath){
  
  suppressMessages(library(dplyr))
  
  rpkm <- data.frame()
  for (i in 1:length(sampleName)){
    inFile <- paste0(RPKMDataPath, sampleName[i])
    
    rpkm.tmp <- read.table(inFile, sep = "\t", header = T, stringsAsFactors = FALSE)
    
    rpkm.tmp <- rpkm.tmp[,c("Gene.ID", "Gene.Name", "FPKM")]
    #for duplicated IDs, just sum up.
    rpkm.tmp <- aggregate(FPKM~Gene.ID+Gene.Name, data = rpkm.tmp, FUN=sum)
    
    if(nrow(rpkm) == 0){
      rpkm <- rpkm.tmp[,c("Gene.ID", "Gene.Name", "FPKM")] 

    } else {
      rpkm <- left_join(rpkm, rpkm.tmp[,c("Gene.ID", "FPKM")], by = c("Gene.ID"))
    }
  }
  colnames(rpkm) <- c("id", "name",
                      as.character(sampleName))
  return(rpkm)
}

countsToDEseq2FDR <- function(counts, CGroup = 0, TGroup = 0, 
                              min_read = 10){
  
  suppressMessages(library(DESeq2))
  
  counts <- counts[apply(counts, 1, function(x){sum(x)}) > min_read, ]
  groups <- factor(c(rep("CGroup",CGroup),rep("TGroup",TGroup)))
  sampleInfo <- data.frame(groups, row.names = colnames(counts))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleInfo, design = ~ groups)
  dds$groups = relevel(dds$groups,ref="CGroup")
  dds <- DESeq(dds)
  res <- results(dds,independentFiltering=F)
  
  #normalized counts obtained based on DESeq2 "plotCounts" source code
  normCounts <- counts(dds, normalized = TRUE)
  results <- as.data.frame(cbind(normCounts, res))
  results$id <- row.names(results)  
  return (results)
}

classifyDEG <- function(df, ctr.rpkm, trt.rpkm, FDR.col, log2FC.col,  
                        log2FC = 1, RPKM = 1, FDR = 0.05){
  
  df$ctr.mean <- rowMeans(df[ctr.rpkm])
  df$trt.mean <- rowMeans(df[trt.rpkm])
  group <- rep(0, nrow(df)) 
  
  for (i in 1:nrow(df)){
    if(df$ctr.mean[i] > RPKM | df$trt.mean[i] > RPKM){
      if(df[log2FC.col][i,] > log2FC & df[FDR.col][i,] < FDR){
        group[i] <- "up-regulated"
      }
      else if (df[log2FC.col][i,] < -log2FC & df[FDR.col][i,] < FDR){
        group[i] <- "down-regulated"
      }
      else{
        group[i] <- "similar_level"
      }
    } 
    else{
      group[i] <- "low_expression_level"
    }
  }
  return(group)
}

ggScatterplot <- function(df, x, y, group, gene,
                          my.color=c("blue", "grey50", "grey50", "red3"),
                          label.up, label.down, xlab, ylab, title,
                          genes4Label = NULL,
                          FC.line = 2){
  
  suppressMessages(library(ggplot2))
  suppressMessages(library(cowplot))
  suppressMessages(library(ggrepel))
  suppressMessages(library(gridExtra))
  suppressMessages(library(grid))
  
  g <- ggplot(df, aes_string(x=x, y=y, color=group,label=gene)) + 
    geom_point() + scale_color_manual(values=my.color) + 
    xlab(xlab) + ylab(ylab) + ggtitle(title) + theme_cowplot(13) + 
    theme(legend.position = "none") + xlim(0, 20) + ylim(0, 20) + 
    geom_abline(intercept = log2(FC.line), slope = 1, linetype = 2) + 
    geom_abline(intercept = -log2(FC.line), slope = 1, linetype = 2) + 
    annotate("text", x = -Inf, y = Inf, label = label.up, 
             col = "red3", size = 3.8, hjust = -0.2, vjust = 1.5) + 
    annotate("text", x= Inf, y =-Inf, label = label.down, 
             col = "blue", size = 3.8, hjust = 1.2, vjust = -1) + 
    geom_text_repel(
      data = subset(df, df[, gene] %in% genes4Label),
      size = 3.8, segment.size = 0.3, segment.color = "black",
      direction = "x", nudge_y = 4.5, nudge_x =-3.5,
      point.padding = 0.25, box.padding = 0.25) +
    geom_point(
      data = subset(df, df[, gene] %in% genes4Label), 
      col = "black", size = rel(1.5))
  return(g)  
}