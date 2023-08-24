
##DESeq2 analysis

x<-read.table("mydata.txt",sep="\t",header=T, row.names=1)
x<-as.matrix(x)
y<-read.table("coldata.txt",sep="\t",header=T, row.names=1)

dds<-DESeqDataSetFromMatrix(countData=x,colData=y,design=~infection_status)
dds<-DESeq(dds)

keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

contrast<-c("infection_status","infected","uninfected")
res<-results(dds,contrast) 

DE_results<-results(dds,contrast)

##Volcano plots

# Custom colors for volcano plot
keyvals <- ifelse(
  res$log2FoldChange < -1.2 & res$pvalue < 1e-2, '#0099FF',
  ifelse(res$log2FoldChange > 1.2 & res$pvalue < 1e-2, 'red',
         'grey'))
         keyvals[is.na(keyvals)] <- 'grey'
           names(keyvals)[keyvals == 'red'] <- 'high'
           names(keyvals)[keyvals == 'grey'] <- 'mid'
           names(keyvals)[keyvals == '#0099FF'] <- 'low'

# Code to render volcano plot
plot_final <- EnhancedVolcano(DE_results,
                              lab = rownames(DE_results),
                              x = 'log2FoldChange',
                              y = 'pvalue',
                              title = 'CD8 T cells',
                              pCutoff = 1e-2,
                              FCcutoff = 1.2,
                              pointSize = 3.0,
                              labSize = 4.0,
                              selectLab = rownames(res)[which(names(keyvals) %in% c('high', 'low'))], 
                              colCustom = keyvals, 
                              colAlpha = 0.9, 
                              cutoffLineType = 'longdash', 
                              cutoffLineCol = 'black', 
                              gridlines.major = FALSE,
                              gridlines.minor = FALSE, 
                              legendPosition = 'none', 
                              legendLabSize = 14,
                              legendIconSize = 4.0,
                              axisLabSize = 12) 

summary(res$log2FoldChange > 1.2 & res$pvalue < 1e-2)
summary(res$log2FoldChange < -1.2 & res$pvalue < 1e-2)

## Heatmaps

# Filtering for top differentially expressed genes
sigs <- na.omit(DE_results)
sigs <- sigs[sigs$padj < 0.1,] # p-value cut-off
sigs <- sigs[sigs$log2FoldChange > 3,] # log fold change cut-off
sigs <- sigs[sigs$baseMean > 500,] # min mean gene count cut-off

summary(sigs)
df <- as.data.frame(sigs)

# Normalizing gene expression data

rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object
matrix <- assay(rlog_out)[rownames(df), rownames(y)] #selects top genes 
matrix.scaled <- t(apply(matrix, 1, scale)) # creates z-score for each column
colnames(matrix.scaled) <- colnames(matrix)

# Custom colors for heatmap
ann <- data.frame(y$infection_status)
colnames(ann) <- c("infection_status")
colours <- list("infection_status" =c("infected"="red", "uninfected"="#0099FF"))
colAnn <- HeatmapAnnotation(df=ann, which="col",col=colours, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))

# Code to render heatmap
h <- Heatmap(matrix.scaled, 
             cluster_rows = T,
             cluster_columns = T,
             name = "Z-score",
             show_column_names = FALSE, 
             column_km = 2,
             top_annotation = colAnn)

## PCA plots
rld <- rlog(dds)

PCA <- plotPCA(rld,
               intgroup = c('cell_type', 'infection_status'),
               returnData = FALSE)

# PCA loadings
rv = rowVars(assay(rld)) 
select = order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pc = prcomp(t(assay(rld)[select,]))

loadings = as.data.frame(pc$rotation)
write.table(as.data.frame(pc$rotation),file="PCA_loadings.txt",sep="\t")
