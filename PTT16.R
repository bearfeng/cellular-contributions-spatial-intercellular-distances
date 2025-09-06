rm(list=ls())
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(tidyverse)
library(stringr)
library(harmony)
library(data.table)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library("beyondcell")
library(SeuratObject) 

setwd("F:/keti/奥芬达唑")
scRNA <- readRDS("renamesOS.rds")
meta <- scRNA@meta.data

scRNA<- JoinLayers(scRNA, what = "all")
meta <- scRNA@meta.data
table(scRNA$celltype)
Idents(scRNA) <- "celltype"

OC <- subset(scRNA,ident = c("Osteoclasts"))
table(OC$celltype)
OC<- JoinLayers(OC, what = "all")
sc <-OC
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 4500)
sc <- ScaleData(sc, verbose = FALSE)
sc <- RunPCA(sc, npcs = 30, verbose = FALSE)
pc.num = 1:30
sc <- sc %>%
  RunUMAP(dims = pc.num) %>%#  RunTSNE( dims = pc.num) %>%
  FindNeighbors(dims = pc.num)
save(sc,file = "OC_subcluster.Rdata")

load(file = "OC_subcluster.Rdata")
scRNA <- sc

# 初始化总结果数据框
results <- data.frame(resolution = numeric(),
                      cluster = numeric(),
                      up_ratio = numeric(),
                      down_ratio = numeric())

# 循环处理不同分辨率
for (res in seq(0.01, 0.02, 0.01)) {
  
  # 设置差异基因筛选条件
  logFCfilter = 0.25
  adjPvalFilter = 0.05
  
  # 根据当前分辨率运行FindClusters
  scRNA <- FindClusters(scRNA, resolution = res)
  
  # 更新seurat_clusters列
  scRNA@meta.data$seurat_clusters <- scRNA@meta.data[[paste0("RNA_snn_res.", res)]]
  scRNA@meta.data$seurat_clusters <- as.factor(scRNA@meta.data$seurat_clusters)
  Idents(scRNA) <- "seurat_clusters"
  
  # 找到所有差异基因
  pbmc.markers <- FindAllMarkers(object = scRNA, assay = "RNA", slot = "data", 
                                 only.pos = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)
  
  # 筛选符合条件的差异基因
  sig <- pbmc.markers[
    abs(as.numeric(as.vector(pbmc.markers$avg_log2FC))) >= 1 & 
      as.numeric(as.vector(pbmc.markers$p_val_adj)) < 0.05, 
  ]
  
  gs<-GetCollection(PSc)
  # 获取signature的上调和下调基因列表
  signature_up <- gs@genelist$`sig-2647`$up
  signature_down <- gs@genelist$`sig-2647`$down
  
  # 提取所有群体
  clusters <- unique(sig$cluster)
  
  # 循环按群提取上调和下调基因的交集及占比
  for (cluster in clusters) {
    
    # 提取该群的上调基因
    upregulated <- sig[sig$cluster == cluster & sig$avg_log2FC > 0, ]
    upregulated_genes_list <- upregulated$gene  # 获取上调基因的名称
    
    # 提取该群的下调基因
    downregulated <- sig[sig$cluster == cluster & sig$avg_log2FC < 0, ]
    downregulated_genes_list <- downregulated$gene  # 获取下调基因的名称
    
    # 计算上调基因与signature_up的交集
    up_intersect <- intersect(upregulated_genes_list, signature_up)
    up_ratio <- length(up_intersect) / length(signature_up)
    
    # 计算下调基因与signature_down的交集
    down_intersect <- intersect(downregulated_genes_list, signature_down)
    down_ratio <- length(down_intersect) / length(signature_down)
    
    # 每次循环创建一个新的result1矩阵
    result1 <- data.frame(resolution = res,
                          cluster = cluster,
                          up_ratio = up_ratio,
                          down_ratio = down_ratio)
    
    # 将result1与总的results数据框合并
    results <- rbind(results, result1)
  }
}

# 输出最终的results数据框
print(results)

# 你可以将最终结果导出为CSV文件
write.csv(results, "all_resolutions_results.csv", row.names = FALSE)


library("beyondcell")
library("Seurat")
library("ggplot2")
library(Seurat)
library(SeuratObject) 

set.seed(4622)
rm(list = ls())
setwd("F:/keti/奥芬达唑")
load(file = "OC_subcluster.Rdata")

scRNA <- sc
# 设置差异基因筛选条件
logFCfilter = 0.25
adjPvalFilter = 0.05
# 根据当前分辨率运行FindClusters
res <- 0.01#上面确定的分辨率
scRNA <- FindClusters(scRNA, resolution = res)

sc <- scobj
sc[["RNA"]] <- as(object =sc[["RNA"]], Class = "Assay")
DefaultAssay(sc) <- "RNA"

gs<-GetCollection(PSc)
bc <- bcScore(sc, gs, expr.thres = 0.01) 
bc <- bcUMAP(bc, k.neighbors = 4, res = 0.2)
bc <- bcUMAP(bc, pc = 10, k.neighbors = 4, res = 0.2)
bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.2")
bc <- bcRanks(bc, idents = "Octype")

celltype <- bc@ranks$Octype
write.csv(celltype,file = "beyoudcellcelltype.csv")
cluster <-  bc@ranks$bc_clusters_res.0.2
write.csv(cluster,file = "beyoudcellcluster.csv")


pdf("1.治疗分群.pdf", width = 7, height = 4) 
bcClusters(bc, UMAP = "beyondcell", idents = "Octype", pt.size = 2)
dev.off()

bcClusters(bc, UMAP = "Seurat", idents = "seurat_clusters")
pdf("2.4方图.pdf", width =10 , height = 10)
bc4Squares(bc, idents = "Octype", lvl="C2_Mature osteoclasts",top = 3)
dev.off()

library(Matrix)

pdf("3.top 1.pdf", width =4 , height = 4)
bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = c("sig-2649", "sig-10651")), pt.size = 1.5)[1]
dev.off()


pdf("4.top 1 Histogram.pdf", width =3.5 , height = 7)
bcHistogram(bc, signatures = "sig-2649", idents = "Octype")
dev.off()



