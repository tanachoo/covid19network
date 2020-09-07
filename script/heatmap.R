# Author: yoshi
# Date: 5/24/2019
# Update: 12/06/2019
# Project: Ecv
# Script: heatmap
# Usage: Rscript heatmap.R input_filename output_filename dist_method hclust_method
# Rscript heatmap.R ECv_STAD_th050.txt heatmap_ECv_th050_STAD_eucl_complete.png euclidean complete
# You need to set 4 arguments when running this script.
# Clustering distance method: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
# Clustering method: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"

library(gplots)
library(viridis)

# Set args
args <- commandArgs(trailingOnly=T)
input <- args[1]
output <- args[2]
dist_method <- args[3]
hclust_method <- args[4]

# Load data
message(paste('Load: ', input, sep=''))
data <- read.table(input, header=T, row.name=1, sep='\t', quote='', as.is=T, check.names=F)
data <- round(data, digits=6)
message('ECv data matlix:')
print(dim(data))

## Clustering
message(paste('distance method: ', dist_method, sep=''))
message(paste('clustering method: ', hclust_method, sep=''))

# clustering oriented gene rows
message("Cluster genes...")
dist <- dist(data, method=dist_method)
hclust_row <- hclust(dist, method=hclust_method)
#hclust_cutree_row <- cutree(hclust_row, k=2) # To stain color at rowside  
#hclust_row_color <- c("#66ff00", "#db36a4") # Select color (green, magenta)
plot(hclust_row, hang=-1)

# clustering oriented sample columns
message("Cluster samples...")
transpose_data <- t(data)
Dist <- dist(transpose_data, method=dist_method)
hclust_col <- hclust(Dist, method=hclust_method)
#hclust_cutree_col <- cutree(hclust_col, k=2) # To stain color at rowside  
#hclust_col_color <- c("#0080ff", "#ff8000") # Select color (blue, orange)
#plot(hclust_col, hang=-1)

## Depict heatmap
message("Depict heatmap...")
row_gene_dend <- as.dendrogram(hclust_row)
col_patient_dend <- as.dendrogram(hclust_col)
mat <- as.matrix(data)
png(output, height=9000, width=12000, res=1500) #height=9000, width=12000, res=800-1500
heatmap.2(mat, col=magma, scale='none', Rowv=row_gene_dend, Colv=col_patient_dend,
          dendrogram='both', density.info='none', trace='none', cexRow=0.1, cexCol=0.1, key=TRUE,
          keysize=1.0, key.xlab='ECv')
#heatmap.2(mat, col=magma, scale='none', Rowv=row_gene_dend, Colv=col_patient_dend,
#          dendrogram='both', density.info='none', trace='none', cexRow=0.1, cexCol=0.1, key=TRUE,
#          RowSideColors=hclust_row_color[hclust_cutree_row],
#          ColSideColors=hclust_col_color[hclust_cutree_col], keysize=1.0, key.xlab='ECv')
#col=colorpanel(256, '#5614B0', '#DBD65C')
##col=viridis, magma, plasma, inferno
#RowSideColors=hclust_row_color[hclust_cutree], ColSideColors=EMTgroup_color
#legend("topleft", legend=c("EMT high", "EMT middle", "EMT low"), fill=c("#f7797d", "#FBD786", "#C6FFDD"), border=FALSE, bty="n", y.intersp=0.7, cex=0.7)
#legend("topleft", legend=c("high", "low"), fill=c("yellow", "blue"), border=FALSE, bty="n", y.intersp=0.7, cex=0.7)
#legend("topleft", legend=c("group 1", "group 2", "ECv high", "ECv low"), fill=c("#0080ff", "#ff8000", "#66ff00", "#db36a4"), border=FALSE, bty="n", y.intersp=0.7, cex=0.7)
dev.off()
message(paste('Save file as: ', output, sep=''))
# message('== fin ==')

