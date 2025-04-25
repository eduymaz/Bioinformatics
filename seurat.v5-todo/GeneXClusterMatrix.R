'''
The pipeline argument names implemented in the Seurat pipeline web address were used. 
For easy and efficient reading of the results, two matrices clustered by gene names and 
cell barcodes are presented as a single matrix.
'''

library(dplyr)
library(Seurat)
library(patchwork)

# This function can be used when apply Seurat pipeline according to Pbmc data on Seurat Web Address.
SeuratObjectTransGenexClusterMatrix = function(){
  
  ClusterList = list()
  ClusterFiltered = list()
  MergedResult = data.frame()
  initialFlag = TRUE
  
  for (i in (0:8)){
    RowInfo = rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters == i,])
    ClusterFiltered = df_GeneToCellName[,RowInfo]
    row_sum <- as.data.frame(apply(ClusterFiltered, MARGIN = 1, sum))
    colnames(row_sum) = paste("Cluster",i)
    if(initialFlag){
      initialFlag = FALSE
      MergedResult = row_sum
      next
    }
    
    MergedResult = cbind(MergedResult, row_sum)
  }
  return(MergedResult)
}

#   Ex. Result
#    
#   Two dataframes are merged; cell_x_Cluster and GeneToCellName
#
#         Cluster0 Cluster1 Cluster2 Cluster3 Cluster4 ....
#
#  Gene1    34        65      56        76      90
#  Gene2    45        67      85        76      23
#  Gene3    32        55      56        76      93
#  Gene4
#  .
#  .      