## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
library(igraph)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', out.width = '90%')

rm(list=ls())
## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(5249) # for reproducibility, remove for normal use

## -------------------------------------------------------------------------------------------------------------------
# loading the data
PG_pro <- read.table("./normalised_data/Drought/drought_normalised_PG_protein.tsv", sep = "\t",
                      header = T, check.names = F, row.names = 1)

PG_lip <- read.table("./normalised_data/Drought/drought_normalised_PG_lipid_small.tsv", sep = "\t",
                      header = T, check.names = F, row.names = 1)

metadata <- read.csv("working_data/Drought/metadata_drought.csv")

#setting output directory
outdir <- "./figures/Drought/PG"

# assigning data
X <- PG_pro
Y <- PG_lip

dim(X) # check the dimensions of the X dataframe
dim(Y) # check the dimensions of the Y dataframe

## ---- eval = FALSE--------------------------------------------------------------------------------------------------
png(file.path(outdir, "Drought_initial_cor_heatmap.png"), width = 2000, height = 2000, res = 300)
imgCor(X, Y, sideColors = c("purple", "green")) # produce a heat map of the cross correlation matrix
dev.off()

## ---- out.width='80%', fig.cap= "FIGURE 2: Heatmap of lambda1 and lambda2 values coloured by the resulting cross-validation score"----
grid1 <- seq(0.001, 0.2, length = 10) # set grid search values for each regularisation parameter
grid2 <- seq(0.001, 0.2, length = 10)

cv.tune.rcc.drought <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, validation = "loo") # optimise the regularisation parameter values


## -------------------------------------------------------------------------------------------------------------------
sink(file = "cv.tune.rcc.drought_PG.txt")
cv.tune.rcc.drought # examine the results of CV tuning
sink()

opt.l1 <- cv.tune.rcc.drought$opt.lambda1 # extract the optimal lambda values
opt.l2 <- cv.tune.rcc.drought$opt.lambda2

CV.rcc.drought <- rcc(X, Y, method = "ridge", lambda1 = opt.l1, lambda2 = opt.l2) # formed optimised CV rCCA


## -------------------------------------------------------------------------------------------------------------------
# shrink.rcc.drought <- rcc(X,Y, method = 'shrinkage') # run the rCCA method using shrinkage
# shrink.rcc.drought$lambda # examine the optimal lambda values after shrinkage 

# plotting indiv graph
png(file.path(outdir, "plot_indiv_CV.png"), width = 2000, height = 2000, 
      res = 300)
plotIndiv(CV.rcc.drought, comp = 1:2, # plot the projection of samples for shrinkage rCCA data
          ind.names = metadata$Sample,
          group = metadata$Conditon, rep.space = "XY-variate", # used averaged variate subspace
          legend = TRUE, title = 'Drought PG, rCCA CV')
dev.off()

# plotting the network, cutoff = 0.5
network(CV.rcc.drought, comp = 1:2, interactive = FALSE,
lwd.edge = 0.2, size.node = 0.5,alpha.node = 0.5,row.names = F, col.names = F,
graph.scale = 0.3,cex.edge.label = 0.6,keysize.label = 0.5,
cutoff = 0.5,save = "png", name.save = "PG_Drought_network_CV_0.50")

# plotting the network, cutoff = 0.7
network(CV.rcc.drought, comp = 1:2, interactive = FALSE,
        lwd.edge = 0.2, size.node = 0.5,alpha.node = 0.5,row.names = F, col.names = F,
        graph.scale = 0.3,cex.edge.label = 0.6,keysize.label = 0.5,
        cutoff = 0.7,save = "png", name.save = "PG_Drought_network_CV_0.70")

# plotting the network, cutoff = 0.9
network(CV.rcc.drought, comp = 1:2, interactive = FALSE,
        lwd.edge = 0.2, size.node = 0.5,alpha.node = 0.5,row.names = F, col.names = F,
        graph.scale = 0.3,cex.edge.label = 0.6,keysize.label = 0.5,
        cutoff = 0.9,save = "png", name.save = "PG_Drought_network_CV_0.90")

# plotting the network no cutoff
network(CV.rcc.drought, comp = 1:2,
        lwd.edge = 0.2, size.node = 0.5,alpha.node = 0.5,row.names = F, col.names = F,
        graph.scale = 0.3,cex.edge.label = 0.6,keysize.label = 0.5,
        cutoff = 0, plot.graph = F,interactive = FALSE,
        save = "png", name.save = "PG_Drought_network_CV_0")

# gettinig files out for cytoscape
# cutoff = 0.5
pdf(NULL); on.exit(dev.off(), add = TRUE)  
net_0.5 <- network(CV.rcc.drought,
               comp = 1:2,
               cutoff = 0.5,
               interactive = FALSE,
               plot.graph = FALSE,
               name.save = NULL)
dev.off()

# cutoff = 0.7
pdf(NULL); on.exit(dev.off(), add = TRUE)  
net_0.7 <- network(CV.rcc.drought,
               comp = 1:2,
               cutoff = 0.7,
               interactive = FALSE,
               plot.graph = FALSE,
               name.save = NULL)

dev.off()

# cutoff = 0.9
pdf(NULL); on.exit(dev.off(), add = TRUE)  
net_0.9 <- network(CV.rcc.drought,
               comp = 1:2,
               cutoff = 0.9,
               interactive = FALSE,
               plot.graph = FALSE,
               name.save = NULL)

dev.off()

# cut 0 
pdf(NULL); on.exit(dev.off(), add = TRUE)  
net_0 <- network(CV.rcc.drought,
               comp = 1:2,
               cutoff = 0,
               interactive = FALSE,
               plot.graph = FALSE,
               name.save = NULL)

dev.off()

# saving the network files
g_0.5 <- net_0.5$gR  
# export
write_graph(g_0.5,
            file = file.path(outdir, "Drought_network_CV_0.5.graphml"),
            format = "graphml")

g_0.7 <- net_0.7$gR  
# export
write_graph(g_0.7,
            file = file.path(outdir, "Drought_network_CV_0.7.graphml"),
            format = "graphml")

g_0.9 <- net_0.9$gR  
# export
write_graph(g_0.9,
            file = file.path(outdir, "Drought_network_CV_0.9.graphml"),
            format = "graphml")


g_0 <- net_0$gR  
# export
write_graph(g_0,
            file = file.path(outdir, "Drought_network_CV_0.graphml"),
            format = "graphml")

# plotting heat map
png(file.path(outdir, "Drought_cim_CV_0.9.png"), width = 3400, height = 2400, res = 300)
cim(CV.rcc.drought, comp = 1:2, ylab = "proteins", xlab = "lipids",
              cutoff = 0.9, margins = c(12, 14))
dev.off()

# plotting full heat map
png(file.path(outdir, "Drought_cim_CV.png"), width = 3400, height = 2400, res = 300)
cim(CV.rcc.drought, comp = 1:2, ylab = "proteins", xlab = "lipids", margins = c(12, 14),
    cutoff = 0)
dev.off()


# getting the data from the heatmap
# for the full heat map, cutoff = 0
pdf(NULL); on.exit(dev.off(), add = TRUE)  
cim_no_cutoff <- cim(CV.rcc.drought, comp = 1:2, ylab = "proteins",
                     xlab = "lipids",
                     margins = c(12, 14),
                     cutoff = 0)

# Convert vector to matrix
mat <- matrix(cim_no_cutoff$mat,
              nrow = length(cim_no_cutoff$row.names),
              ncol = length(cim_no_cutoff$col.names),
              byrow = FALSE,
              dimnames = list(cim_no_cutoff$row.names, cim_no_cutoff$col.names))

# Save to CSV if needed
write.csv(mat, file = file.path(outdir, "cim_no_cutoff_matrix_drought.csv"))

# with cutoff of 0.9
cim_0.9_cutoff <- cim(CV.rcc.drought, comp = 1:2, ylab = "proteins", 
                      xlab = "lipids",
                      cutoff = 0.9,
                      margins = c(12, 14))

# Convert vector to matrix
mat_0.9_cut <- matrix(cim_0.9_cutoff$mat,
                      nrow = length(cim_0.9_cutoff$row.names),
                      ncol = length(cim_0.9_cutoff$col.names),
                      byrow = FALSE,
                      dimnames = list(cim_0.9_cutoff$row.names, cim_0.9_cutoff$col.names))

# Save to CSV if needed
write.csv(mat_0.9_cut, file = file.path(outdir, "cim_0.9_cutoff_matrix_drought.csv"))

# get the loadings
# Export protein (X) loadings
load_X <- as.data.frame(CV.rcc.drought$loadings$X)
write.csv(load_X, file = file.path(outdir, "Proteins_drought_all_loadings.csv"), 
          row.names = TRUE)

# Export lipid (Y) loadings
load_Y <- as.data.frame(CV.rcc.drought$loadings$Y)
write.csv(load_Y, file = file.path(outdir, "Lipids_all_loadings.csv"), 
          row.names = TRUE)


# plotting indvidual plots
png(file.path(outdir, "plot_indiv.png"), width = 3400, height = 2400, res = 300)
plotIndiv(CV.rcc.drought, comp = 1:2, 
          ind.names = metadata$Sample,
          group = metadata$Conditon, rep.space = "XY-variate", 
          legend = TRUE, title = 'Multivariant projections based on rCCA analysis')
dev.off()

# plotting arrow plots
png(file.path(outdir, "plot_arrow.png"), width = 3400, height = 2400, res = 300)
plotArrow(CV.rcc.drought, group = metadata$Conditon, 
          col.per.group = color.mixo(1:6),
          title = 'Multivariant projections based on rCCA analysis') 
dev.off()
  
