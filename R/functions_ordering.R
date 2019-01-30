
#' Order Components by PseudoTime
#'
#' @param component string; Component name as in datat
#'
#' @details 
#' The order is determined by using a weighted mean of the pseudotime values, where the weights are the cell scores of 
#' the component. Only cells with an absolute cell score of greater than 2 contribute to the mean.
#' 
#' @return weighted mean of pseudotime for a given component
#' @export

ptorder <- function(component, cell_metadata=datat, threshold=2){
  tmp <- cell_metadata[!is.na(PseudoTime)][abs(get(component)) > threshold][,.(PseudoTime, abs(get(component)))]
  return(weighted.mean(tmp$PseudoTime, tmp$V2))
}


#' Bootstrap Principal Curve Analysis
#'
#' @param input numeric matrix; a matrix of points in arbitrary dimension, e.g. column for PC1/Tsne1 and column for PC2/Tsne2
#' @param df_start integer; degrees of freedom to start at
#' @param df_end integer; degrees of freedom to finish at
#' @param plot logical; if TRUE (default) curve is plotted at end of each df iteration
#'
#' @details 
#' A principal curve is fitted first with df = df_start, then again for df_start+1 iteratively until reaching df_end
#' 
#' @return list object containing results of each sucessive princurve
#' @export
#' @import princurve

princurve_bootstrap <- function(input, df_start=4, df_end=10, plot=TRUE){
  
  principal_curves <- list()
  
  df_range <- df_start:df_end
  
  for(i in seq_along(df_range)){
    cat(paste("Fitting princurve with df", df_range[i],"\n"))
    
    if(df_range[i]==df_start){
      principal_curves[[i]] <- principal.curve(input, df = df_range[i])
    }else{
      principal_curves[[i]] <- principal.curve(input, df = df_range[i], start = principal_curves[[i-1]])
    }
    
    if(plot){
      plot(principal_curves[[i]])
    }
    
  }
  
  names(principal_curves) <- paste0("df_", df_range)
  
  return(principal_curves)
  
}

#' Load component orderings in global environ
#' 
#' @return !! This function modifies (populates) the global environment with objects relating to component names and orders !!
#' 
#' @export
#' @import data.table

load_component_orderings <- function(){

component_labels <<- c(
  "Lymphocytes"=3,
  "Macrophages"=11,
  "Fibroblast like"=32,
  "Leydig"=40,
  "Leydig (Kallikreins)"=19,
  "Leydig (Gstm3 high)"=24,
  "Leydig (Fabp3 low)"=26,
  "Sertoli & Leydig"=49,
  "Sertoli (Trim7 High)"=37,
  "Sertoli (Aard High)"=45,
  "Sertoli (rare)"=16,
  "Sertoli (rare)"=21,
  "Sertoli/Muscle?"=10,
  "Gfra1 stem cells"=50,
  "Undifferentiated Spermatogonia"=31,
  "Intermediate Spermatogonia"=7,
  "Differentiated Spermatogonia"=33,
  "Pre-leptotene"=2,
  "Leptotene"=5,
  "Leptotene-Zygotene"=44,
  "X activation (Hormad1)"=38,
  "Pre-Pachytene (& Hormad1)"=23,
  "Early General"=27,
  "Early Pachytene 1"=13,
  "Early Pachytene 2"=47,
  "Pachytene"=42,
  "Odd Mid Pachytene Clump"=48,
  "Late Pachytene"=39,
  "Meiotic Divisions"=20,
  "Cul4a"=25,
  "Acrosomal"=30,
  "Acrosomal V2"=28,
  "Late Acrosomal"=35,
  "Round Spermatid"=15,
  "Spermiogenesis V2"=36,
  "Spermiogenesis"=17,
  "Late Spermiogenesis 1"=18,
  "Late Spermiogenesis 2"=34,
  "Single Cell (Macrophage)"=14,
  "Single Cell (Macrophage)"=8,
  "Single Cell"=4,
  "Single Cell (Macrophage)"=46,
  "Single Cell (Macrophage)"=1,
  "Ribosomal"=43,
  "Respiration"=9,
  "Batch (Late)"=41,
  "Hormad1 Batch"=12,
  "Hormad1"=29,
  "Chemical Dissociation (WT)"=6,
  "CNP / Batch"=22)

QC_fail_components <<- c(22,6,25,29,12,28,41,1,46,4,8,14,9,43)
somatic_components <<- c(3,11,32,40,19,45,10,24,26,37,49,16,21, 14,8,4,46,1)

component_order_dt <<- data.table(component_number = component_labels, name = names(component_labels))

component_order_dt[, QC_fail := FALSE]
component_order_dt[, Somatic := FALSE]

component_order_dt[component_number %in% QC_fail_components, QC_fail := TRUE]
component_order_dt[component_number %in% somatic_components, Somatic := TRUE]

component_order_dt[, pseudotime_average := sapply(paste0("V",component_number), ptorder)]

component_order_dt[Somatic==TRUE | component_number %in% c(22,43), pseudotime_average := NA]

component_order <<- component_order_dt[!grep("Single",name)][!component_number %in% c(22,43)][order(-pseudotime_average, na.last = F)]$component_number
component_order_all <<- component_order_dt[order(-pseudotime_average, na.last = F)]$component_number

meiotic_component_order <<- component_order_dt[!is.na(pseudotime_average)][order(-pseudotime_average)]$component_number

PN_ratio <- apply(results$scores, 2, function(x) log(abs(min(x)/max(x))))

half_exclusion_comps <<- c(paste0(names(which(PN_ratio < -log(5))),"N"),
                     paste0(names(which(PN_ratio > log(5))),"P"))

ordering <<- c(rbind(paste0("V",component_order,"P"),paste0("V",component_order,"N")))  #component_labels
half_exclusions <<- which(ordering %in% half_exclusion_comps)

tmp <- data.table(component_number=as.numeric(gsub("V|P|N","",ordering)))
setkey(component_order_dt, component_number)

component_names <<- merge(tmp, component_order_dt[,.(component_number, name)],by="component_number", sort=F)$name

}