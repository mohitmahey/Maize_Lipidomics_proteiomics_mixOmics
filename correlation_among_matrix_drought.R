# remove the environment
rm(list = ls())

# loading required libraries
library(pheatmap)

# loading the data
PG_pro <- read.table("./normalised_data/Drought/drought_normalised_PG_protein.tsv", sep = "\t",
                     header = T, check.names = F, row.names = 1)

PG_lip <- read.table("./normalised_data/Drought/drought_normalised_PG_lipid_small.tsv", sep = "\t",
                     header = T, check.names = F, row.names = 1)

thy_pro <- read.table("./normalised_data/Drought/drought_normalised_thy_protein.tsv", sep = "\t",
                      header = T, check.names = F, row.names = 1)

thy_lip <- read.table("./normalised_data/Drought/drought_normalised_thy_lipid_small.tsv", sep = "\t",
                      header = T, check.names = F, row.names = 1)


#setting output directory
outdir <- "./figures/Drought/correlation_matix"

# making all of them as matirx
PG_pro_mat <- as.matrix(PG_pro)
PG_lip_mat <- as.matrix(PG_lip)

thy_lip_mat <- as.matrix(thy_lip)
thy_pro_mat <- as.matrix(thy_pro)

# Correlation between columns
cor_PG_pro <- cor(PG_pro_mat, method = "pearson") 

cor_PG_lip <- cor(PG_lip_mat, method = "pearson") 

cor_thy_pro <- cor(thy_pro_mat, method = "pearson") 

cor_thy_lip <- cor(thy_lip_mat, method = "pearson") 

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

#Thy_lipid-lipid
png(file.path(outdir, "Thy_lipid-lipid.png"), width = 7000, height = 5000, res = 300)
pheatmap(cor_thy_lip, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Thy_lipid-lipid")
dev.off()

# Thy_protein-protein
png(file.path(outdir, "Thy_protein-protein.png"), width = 7000, height = 5000, res = 300)
pheatmap(cor_thy_pro, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Thy_protein-protein",
         show_rownames = F,
         show_colnames = F)
dev.off()


# saving all the cor matrix
write.csv(cor_PG_pro,file.path(outdir, "PG_protein-protein.csv"),
          row.names = TRUE)

write.csv(cor_PG_lip,file.path(outdir, "PG_lipid-lipid.csv"),
          row.names = TRUE)

write.csv(cor_thy_lip,file.path(outdir, "Thy_lipid-lipid.csv"),
          row.names = TRUE)

write.csv(cor_thy_pro,file.path(outdir, "Thy_protein-protein.csv"),
          row.names = TRUE)

