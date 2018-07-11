load("SDA_objects.rds", envir = .GlobalEnv)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
data <- as.matrix(data)
data <- as.big.matrix(data, backingfile = "data.big")
data <- attach.big.matrix("data.big.desc")

datat$Tsne1 <- NULL
datat$Tsne1_QC2 <- NULL
datat$Tsne2 <- NULL
datat$Tsne2_QC2 <- NULL
datat$experiment <- NULL
datat$library_size <- NULL
datat$somatic <- NULL
datat$somatic2 <- NULL
datat$somatic3 <- NULL
datat$somatic4 <- NULL
datat$sex <- NULL
datat$sex_predictions <- NULL
datat$autosomal <- NULL
datat$autosomal_predictions <- NULL
datat$X_predictions <- NULL
datat$Y_predictions <- NULL
datat$PseudoTime1 <- NULL
datat$PseudoTime2 <- NULL
datat$hclust_group <- NULL
datat$`mt-Rnr2` <- NULL
datat$`mt-Rnr2_pred` <- NULL


# keep these as easter eggs
#datat$group <- NULL
#datat$msci_ratio <- NULL
#datat$PseudoTime <- NULL

str(datat)

principal_curves <- principal_curves["df_9"]

str(principal_curves[["df_9"]])

library(data.table)
library(testisAtlas)
load_component_orderings()
#source("functions_ordering.R", local = TRUE)

save(chromosome.lengths, datat, principal_curves, results, rna_locations, component_order_dt, file = "SDA_objects_small.rds")