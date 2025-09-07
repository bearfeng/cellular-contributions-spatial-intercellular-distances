library(tidyverse)
library(ggpubr)
library(ggplot2)
set.seed(4622) 
#Read Seurat object
PT_5_OS_and_ohters <- readRDS("./CARD/PT_5_OS_and_ohters.rds")
setwd("./CARD/out")
sc_data <- as.data.frame(GetAssayData(object = PT_5_OS_and_ohters, 
                                      assay = "RNA",  
                                      slot = "counts"))
cell_metadata <- PT_5_OS_and_ohters@meta.data
cell_metadata$CellID <- rownames(cell_metadata)
cell_metadata$CellType <- cell_metadata$cellType
cell_types <- unique(cell_metadata$CellType)
CARD_obj <- readRDS("./CARD/PT09_SP_CARD_obj.rds") 
prop_matrix <- as.data.frame(CARD_obj@Proportion_CARD)

# Calculate the average expression level of each cell type
avg_exp <- as.data.frame(sapply(cell_types, function(ct){
  cells <- cell_metadata$CellID[cell_metadata$CellType == ct]
  if(length(cells) > 0) rowMeans(sc_data[, cells, drop=FALSE])
  else numeric(length(genes))
}))
# Only retain the target gene
target_genes <- c("IBSP", "TWIST1", "SPARC")
avg_exp_target <- avg_exp[target_genes, ]


# Z-score standardization
z_score <- function(x) (x - mean(x)) / sd(x)
z_scores <- apply(avg_exp_target, 1, function(gene_exp){
  z_score(gene_exp)
}) %>% as.data.frame()

# Calculate the comprehensive score for each cell type
z_scores$CompositeScore <- z_scores$IBSP*6/38 + z_scores$TWIST1*19/38 + z_scores$IBSP*13/38

print("细胞类型综合得分:")
print(z_scores)
prop_df <- as.data.frame(prop_matrix)
prop_df$SpotID <- rownames(prop_df)

# Calculate the contribution of each cell type to each spot
contribution_df <- prop_df %>%
  pivot_longer(-SpotID, names_to = "CellType", values_to = "Proportion") %>%
  left_join(z_scores %>% select(CompositeScore) %>% 
              rownames_to_column("CellType"), by = "CellType") %>%
  mutate(RawContribution = Proportion * CompositeScore) %>%
  group_by(SpotID) %>%
  mutate(RelativeContribution = RawContribution / sum(RawContribution)) %>%
  ungroup() %>%
  select(SpotID, CellType, Proportion, CompositeScore, RelativeContribution)

result <- contribution_df %>%
  group_by(CellType) %>%
  summarise(
    Avg_Proportion = mean(Proportion, na.rm = TRUE),
    Avg_CompositeScore = mean(CompositeScore, na.rm = TRUE),
    Avg_RelativeContribution = mean(RelativeContribution, na.rm = TRUE)
  ) %>% 
  ungroup() %>%
  arrange(desc(Avg_CompositeScore)) 
name <- result$CellType
result <- result[,-1]
rownames(result) <- name

# Output result
p <- result %>%
  rownames_to_column("CellType") %>%
  pivot_longer(-c(CellType,Avg_RelativeContribution), names_to = "Metrics", values_to = "Z_score") %>%
  ggplot(aes(x = CellType, y = Z_score, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_line(aes(x = CellType, y = Avg_RelativeContribution, group = 1, color = "Avg_RelativeContribution"), 
            size = 1) +
  scale_fill_manual(values = c(
    "Avg_CompositeScore" = "#6788C2",  
    "Avg_Proportion" = "#E0CCE0"   
  ))+
  scale_color_manual(name = "", values = c("Avg_RelativeContribution" = "#DC143C")) +
  labs(title = "PT09_axis_contributions",
       y = "Output values") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Output result
write_csv(contribution_df, "pt09_celltype_contributions.csv")
write_csv(result, "result.csv")

ggsave("Cells_contribution_PT09.png",p,width = 8,height = 6)
ggsave("Cells_contribution_PT09.pdf",p,width = 8,height = 6)
ggsave("Cells_contribution_PT09.jpg",p,width = 8,height = 6)

