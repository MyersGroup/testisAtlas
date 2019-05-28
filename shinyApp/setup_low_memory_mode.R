library(data.table)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
library(testisAtlas)

load2("cache")
load_component_orderings()

data <- as.matrix(data)
data <- as.big.matrix(data, backingfile = "data.big")
data <- attach.big.matrix("data.big.desc")

cell_data$log_library_size <- log10(cell_data$library_size)
cell_data$somatic <- cell_data$somatic4
cell_data <- cell_data[,c("cell","Tsne1_QC1","Tsne2_QC1","PseudoTime","msci_ratio","group","Umap1","Umap2","log_library_size","somatic",paste0("V",1:50)), with=F]

principal_curves <- principal_curves["df_9"]
str(principal_curves[["df_9"]])

save(cell_data, principal_curves, SDAresults, gene_annotations, component_order_dt, file = "cached_objects_small.rds")
