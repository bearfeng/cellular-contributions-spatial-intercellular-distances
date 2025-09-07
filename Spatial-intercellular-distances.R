library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggplot2)
set.seed(4622)
# Read seurat object
PT_5_OS_and_ohters <- readRDS("./CARD/PT_5_OS_and_ohters.rds")
setwd("./CARD/out")
CARD_obj <- readRDS("./CARD/PT09_SP_CARD_obj.rds") 
prop_matrix <- as.data.frame(CARD_obj@Proportion_CARD)
prop_matrix$ID <- rownames(prop_matrix)
location_matrix <- CARD_obj@spatial_location
location_matrix$ID <- rownames(location_matrix)
spatial_data <- merge(location_matrix,prop_matrix, by = "ID")

# Calculate the weighted centroid of each cell type
cell_types <- c("Cluster0", "Cluster1", "Cluster2","B_cells","Endothelial_cells", "Fibroblasts",
                "Low_quality_cells", "Myeloid_cells", "NK/T_cells", "Osteoclasts")

centroids <- sapply(cell_types, function(ctype) {
  total <- sum(spatial_data[[ctype]])
  w_x = sum(spatial_data$x * spatial_data[[ctype]]) / total
  w_y = sum(spatial_data$y * spatial_data[[ctype]]) / total
  c(w_x, w_y)
})

# Calculate the distance between centroids
dist_matrix <- as.matrix(dist(t(centroids)))
rownames(dist_matrix) <- colnames(dist_matrix) <- cell_types

# Result visualization
print(dist_matrix)
p<-pheatmap::pheatmap(dist_matrix, display_numbers = TRUE,
                      color=colorRampPalette(c("#6788C2","white" ,"#1aafd0"))(100),
                      treeheight_row = 20 ,
                      treeheight_col = 20 ,
                      fontsize_number =15, 
                      fontsize_row =15,
                      fontsize_col =15,
                      cellwidth = 40,cellheight = 40, 
                      #border="white" ,
                      angle_col = 45, 
                      
)

ggsave("dist_matrix.jpg",p,width = 10,height = 9)
ggsave("dist_matrix.pdf",p,width = 10,height = 9)
ggsave("dist_matrix.eps",p,width = 10,height = 9)