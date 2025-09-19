# remove the environment
rm(list = ls())

# loading required libraries
library(pheatmap)
library(igraph)

# loading the data
PG_pro <- read.table("./normalised_data/Heat/heat_normalised_PG_protein.tsv", sep = "\t",
                     header = T, check.names = F, row.names = 1)

PG_lip <- read.table("./normalised_data/Heat/heat_normalised_PG_lipid_smal.tsv", sep = "\t",
                     header = T, check.names = F, row.names = 1)

# thy_pro <- read.table("./normalised_data/Drought/drought_normalised_thy_protein.tsv", sep = "\t",
#                       header = T, check.names = F, row.names = 1)
# 
# thy_lip <- read.table("./normalised_data/Drought/drought_normalised_thy_lipid_small.tsv", sep = "\t",
#                       header = T, check.names = F, row.names = 1)
# 

#setting output directory
outdir <- "/Users/mohitmahey/Library/CloudStorage/OneDrive-MichiganStateUniversity/SIDE_PROJECTS/peters_mix_omics/figures/pearson_correlation/correlation_matix/heat"

# making all of them as matirx
PG_pro_mat <- as.matrix(PG_pro)
PG_lip_mat <- as.matrix(PG_lip)

# thy_lip_mat <- as.matrix(thy_lip)
# thy_pro_mat <- as.matrix(thy_pro)

# Correlation between columns
cor_PG_pro <- cor(PG_pro_mat, method = "pearson") 

cor_PG_lip <- cor(PG_lip_mat, method = "pearson") 

# cor_thy_pro <- cor(thy_pro_mat, method = "pearson") 
# 
# cor_thy_lip <- cor(thy_lip_mat, method = "pearson") 

# Heatmap of correlations between columns
# PG_protein-protein
png(file.path(outdir, "PG_protein-protein.png"), width = 3000, height = 2000, res = 300)
pheatmap(cor_PG_pro, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "PG_protein-protein")
dev.off()

#PG_lipid-lipid
png(file.path(outdir, "PG_lipid-lipid.png"), width = 6000, height = 4000, res = 300)
pheatmap(cor_PG_lip, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "PG_lipid-lipid")
dev.off()

# #Thy_lipid-lipid
# png(file.path(outdir, "Thy_lipid-lipid.png"), width = 7000, height = 5000, res = 300)
# pheatmap(cor_thy_lip, 
#          cluster_rows = TRUE, 
#          cluster_cols = TRUE, 
#          color = colorRampPalette(c("blue", "white", "red"))(100),
#          main = "Thy_lipid-lipid")
# dev.off()
# 
# # Thy_protein-protein
# png(file.path(outdir, "Thy_protein-protein.png"), width = 7000, height = 5000, res = 300)
# pheatmap(cor_thy_pro, 
#          cluster_rows = TRUE, 
#          cluster_cols = TRUE, 
#          color = colorRampPalette(c("blue", "white", "red"))(100),
#          main = "Thy_protein-protein",
#          show_rownames = F,
#          show_colnames = F)
# dev.off()


# saving all the cor matrix
write.csv(cor_PG_pro,file.path(outdir, "PG_protein-protein.csv"),
          row.names = TRUE)

write.csv(cor_PG_lip,file.path(outdir, "PG_lipid-lipid.csv"),
          row.names = TRUE)

# write.csv(cor_thy_lip,file.path(outdir, "Thy_lipid-lipid.csv"),
#           row.names = TRUE)
# 
# write.csv(cor_thy_pro,file.path(outdir, "Thy_protein-protein.csv"),
#           row.names = TRUE)

# 1. Remove self-correlations
diag(cor_PG_pro) <- NA

# 2. Pick a correlation threshold
r_thresh <- 0.7

# 3. Get upper triangle indices (avoid duplicates)
ut <- which(upper.tri(cor_PG_pro), arr.ind = TRUE)

# 4. Keep only pairs with |r| >= threshold
keep <- !is.na(cor_PG_pro[ut]) & abs(cor_PG_pro[ut]) >= r_thresh
#keep <- ut

# 5. Build edge list
edges <- data.frame(
  from   = rownames(cor_PG_pro)[ut[keep, "row"]],
  to     = colnames(cor_PG_pro)[ut[keep, "col"]],
  weight = cor_PG_pro[ut][keep],
  sign   = ifelse(cor_PG_pro[ut][keep] >= 0, "pos", "neg"),
  stringsAsFactors = FALSE
)

# 6. Create graph
g <- graph_from_data_frame(edges, directed = FALSE)

# 7. Add a simple node attribute (degree)
V(g)$degree <- degree(g)

# 8. Save network
write_graph(g, "PG_protein-protein_cor_0.7.graphml", format = "graphml")
write.csv(edges, "PG_protein-protein_cor_0.7.csv", row.names = FALSE)


# 1. Remove self-correlations
diag(cor_PG_lip) <- NA

# 2. Pick a correlation threshold
l_thresh <- 0.7

# 3. Get upper triangle indices (avoid duplicates)
ut <- which(upper.tri(cor_PG_lip), arr.ind = TRUE)

# 4. Keep only pairs with |r| >= threshold
keep <- !is.na(cor_PG_lip[ut]) & abs(cor_PG_lip[ut]) >= l_thresh
#keep <- ut

# 5. Build edge list
edges <- data.frame(
  from   = rownames(cor_PG_lip)[ut[keep, "row"]],
  to     = colnames(cor_PG_lip)[ut[keep, "col"]],
  weight = cor_PG_lip[ut][keep],
  sign   = ifelse(cor_PG_lip[ut][keep] >= 0, "pos", "neg"),
  stringsAsFactors = FALSE
)

# 6. Create graph
g <- graph_from_data_frame(edges, directed = FALSE)

# 7. Add a simple node attribute (degree)
V(g)$degree <- degree(g)

# 8. Save network
write_graph(g, "PG_lipid-lipid_cor_0.7.graphml", format = "graphml")
write.csv(edges, "PG_lipid-lipid_cor_0.7.csv", row.names = FALSE)

#write_graph(g, "PG_lipid-lipid_cor_all.graphml", format = "graphml")
#write.csv(edges, "PG_protein-protein_cor_all.csv", row.names = FALSE)
