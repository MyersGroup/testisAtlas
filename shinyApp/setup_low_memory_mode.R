load("cached_objects.rds", envir = .GlobalEnv)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
data <- as.matrix(data)
data <- as.big.matrix(data, backingfile = "data.big")
data <- attach.big.matrix("data.big.desc")

# keep these as easter eggs
#cell_data$group <- NULL
#cell_data$msci_ratio <- NULL
#cell_data$PseudoTime <- NULL
cell_data$log_library_size <- log10(cell_data$library_size)

cell_data$Tsne1 <- NULL
cell_data$Tsne1_QC2 <- NULL
cell_data$Tsne2 <- NULL
cell_data$Tsne2_QC2 <- NULL
cell_data$experiment <- NULL
cell_data$library_size <- NULL
cell_data$somatic <- NULL
cell_data$somatic2 <- NULL
cell_data$somatic3 <- NULL
# keep somatic4 as somatic
cell_data$somatic <- cell_data$somatic4
cell_data$somatic4 <- NULL
cell_data$sex <- NULL
cell_data$sex_predictions <- NULL
cell_data$autosomal <- NULL
cell_data$autosomal_predictions <- NULL
cell_data$X_predictions <- NULL
cell_data$Y_predictions <- NULL
cell_data$PseudoTime1 <- NULL
cell_data$PseudoTime2 <- NULL
cell_data$hclust_group <- NULL
cell_data$`mt-Rnr2` <- NULL
cell_data$`mt-Rnr2_pred` <- NULL

str(cell_data)

principal_curves <- principal_curves["df_9"]

str(principal_curves[["df_9"]])

library(data.table)
library(testisAtlas)
load_component_orderings()

save(chromosome.lengths, cell_data, principal_curves, SDAresults, rna_locations, component_order_dt, file = "cached_objects_small.rds")