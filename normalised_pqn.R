# trying to normalise the data, 
# apparantly  Probabilistic Quotient Normalization should work

# cleaning the environment
rm(list=ls())

# loading the libraries
# List of packages required
packages <- c("tidyverse","preprocessCore", "pcaMethods", 
              "multtest","Rcpm")

# Rcpm has all the functions 
#devtools::install_github("ricoderks/Rcpm")

# if not installed, install and load it
for (package in packages) {
  if (!(package %in% rownames(installed.packages()))) {
    BiocManager::install(package)
  }
  library(package, character.only = TRUE)
}

# loading my data
data <- read.table("working_data/Heat/heat_thy_protein.txt", sep = "\t",
                          header = T, check.names = F, row.names = 1)

data

# transform
data <- t(data)

data[1:5,1:5]

# making it into matrix
data_mat <- as.matrix(data)

# checking the matrix
any(is.na(data_mat))  
any(data_mat == 0)

# Show original matrix size and summary
cat("🟡 BEFORE removing all-zero proteins:\n")
dim(data_mat)  

# Remove proteins (columns) where all values are zero
data_mat <- data_mat[, colSums(data_mat != 0, na.rm = TRUE) > 0]

# Show new matrix size and summary
cat("\n🟢 AFTER removing all-zero proteins:\n")
dim(data_mat)         

data_mat[1:5,1:5]


##### data preprsing looks good

# performing normalisation
# assuming total area normalisation has been done
normalised_data <- pqn(X = data_mat, n = "mean")

normalised_data[1:5,1:5]
#View(normalised_data)

#dir.create("normalised_data")
write.csv(normalised_data,file ="normalised_data/heat_normalised_thy_protein.csv")


### test code
# NEW: also remove zero-variance columns (constant features)
# for deciding variance
# sds <- apply(data_mat, 2, sd, na.rm = TRUE)
# hist(sds, breaks = 50, main = "Feature standard deviations", xlab = "SD")
# summary(sds)
#var_threshold <- 0.1
#var_threshold <- 0
#data_mat <- data_mat[, apply(data_mat, 2, sd, na.rm = TRUE) > var_threshold]
# Step 2: Remove proteins where 80% or more of samples are zero
#data_mat <- data_mat[, (colSums(data_mat == 0, na.rm = TRUE) / nrow(data_mat)) < 0.99]

