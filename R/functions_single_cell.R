
#' Rotate a tSNE embedding about the origin
#'
#' @param tsne object; output of Rtsne::Rtsne()
#' @param angle numeric; degrees to rotate the tSNE embedding by
#'
#' @return The Y matrix from Rtsne output, rotated by angle
#'
#' @export

rotate_tsne <- function(tsne, angle){
  angle <- (-angle * pi) / (-180)
  rotm <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol=2)
  tsne$Y %*% rotm
}

#' Generate rotation matrix / rotate a matrix
#'
#' @param object matrix; matrix with two columns
#' @param angle numeric; degrees to rotate the tSNE/UMAP/PCA embedding by
#' @param rotation logical; if TRUE returns rotation matrix instead of the rotated matrix (default: FALSE)
#'
#' @return Rotated matrix
#'
#' @export

rotate_matrix <- function(object, angle, rotation=FALSE){
  angle <- (-angle * pi) / (-180)
  rotm <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol=2)
  object %*% rotm
}


#' Print tsne based on raw expression rather than SDA scores
#'
#' @param gene string; gene symbol, of which expression values will be used to colour points (cells)
#' @param cell_metadata data.table with columns Tsne1, Tsne2, and PCs
#' @param expression_matrix cell by gene expressionn matrix
#'
#' @return ggplot2 object
#'
#' @export
#' @import ggplot2
print_raw_tsne <- function(gene, expression_matrix=data, cell_metadata=cell_data){
  ggplot(merge(cell_metadata, expression_dt(gene, expression_matrix)), aes(Tsne1, Tsne2, color=get(gene))) +
    geom_point(size=0.2) +
    scale_color_viridis(direction=-1) +
    theme(legend.position = "bottom") +
    ggtitle(paste0(gene," - t-SNE")) + simplify
}




#' Print PCA
#'
#' @param cell_metadata data.table with columns Tsne1, Tsne2, and PCs
#' @param pc string; name of principal compopnent (column of cell_data), of which cell score values will be used to colour points (cells)
#'
#' @return ggplot2 object
#'
#' @export
#' @import ggplot2
print_pca <- function(cell_metadata=cell_data, pc){
  ggplot(cell_metadata, aes(Tsne1, Tsne2, color=get(pc))) +
    geom_point(size=0.2) +
    scale_colour_distiller(palette="YlOrRd", direction=1) +
    theme(legend.position = "bottom") +
    ggtitle(paste0(pc," - t-SNE")) + simplify
}




#' Clustered heatmap (partially deprecated)
#'
#' @param cell_metadata data.table; cell_data, with subset of gene columns
#' @param name string; title
#' @param annotation.col data.table; subset of cell_data, passed to annCol of aheatmap
#' @param colv_order passed to Colv of aheatmap
#' @param col_lab passed to labCol of aheatmap
#' @param row_text_size passed to cexRow of aheatmap
#'
#' @details requires PCA object to be bound to cell_metadata
#'
#' @return aheatmap plot
#'
#' @export
#' @import NMF
clustered_heatmap <- function(cell_metadata=cell_data, name, annotation.col=NULL, colv_order=NULL, col_lab=NULL, row_text_size=0.5){
  
  cols3 <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
  #[, !c("cell", "cell_id"), with = FALSE]
  aheatmap(t(as.matrix(cell_metadata)),
           treeheight = 100,
           hclustfun = "ward.D",
           breaks=-1,
           color = "-RdBu:50",
           cexAnn = 1.5,
           cexRow=row_text_size,
           annCol = annotation.col, 
           annColors = list(Prm2=cols3, Uchl1=cols3, Sycp2=cols3, Nek2=cols3, Hemgn=cols3, Spaca1=cols3, Gapdhs=cols3, Nasp=cols3, Acrv1=cols3, Insl3=cols3, 'sqrt(library_size)'=cols3, library_size=cols3, cluster.SDA=c(brewer.pal(9, "Set1"),"#000000")),
           sub="Cell",
           main=paste(name,"Individual scores"),
           Rowv = NA,
           Colv = colv_order,
           labCol = col_lab
  )
}


#' Plot scores
#'
#' @param component string ; component name as in cell_metadata object e.g. "V12"
#' @param cell_metadata data.table with columns cell, Tsne1_QC1, Tsne2_QC2, and components V1, V2 etc.
#' @param point_size numeric
#'
#' @return ggplot object
#' 
#' @details
#' 
#' @export
#' 
#' @import ggplot2 data.table RColorBrewer

plot_cell_scores <- function(component="V25", point_size=0.6, cell_metadata=cell_data){
  tmp <- cell_metadata[,.(group,get(component))]
  tmp[group=="mj",group:="WT"]
  setnames(tmp, c("group",component))

return(ggplot(tmp, aes(get(component), group, colour=group)) +
  geom_jitter(size=point_size, alpha=0.5, stroke=0) +
  scale_color_brewer(palette = "Set1") + 
  theme_minimal() +
  theme(legend.position = "none", axis.title.y=element_blank() ) +
  xlab(paste0("Component ",component,"\n Cell Score")))
  }

#' Load Curve Data
#'
#' @param princurves list of outputs of princurve()
#' @param principal_curve string; name of list entry in principal_curves object, to use for pseudotime line
#'
#' @return data.table with columns for position of curve, annotated with Stage
#' 
#' @export
#' 
#' @import data.table
load_curve_data <- function(princurves=principal_curves, principal_curve="df_9"){
  
  curve_data <- data.table(princurves[["df_9"]]$s[princurves[["df_9"]]$tag, 1],
                           princurves[["df_9"]]$s[princurves[["df_9"]]$tag, 2])
  curve_data[,PseudoTime := 1:nrow(curve_data)]
  
  m <- nrow(curve_data)
  
  curve_data[m:(m*0.975),Stage := "Spermatogonia"]
  curve_data[(m*0.975):(m*0.96),Stage := "Leptotene"]
  curve_data[(m*0.96):(m*0.93),Stage := "Zygotene"]
  curve_data[(m*0.93):(m*0.62),Stage := "Pachytene"]
  curve_data[(m*0.62):(m*0.5),Stage := "Division II"]
  curve_data[(m*0.5):(m*0.2),Stage := "Round Spermatid"]
  curve_data[(m*0.2):0,Stage := "Elongating \nSpermatid"]
  
  curve_data[,Stage := factor(Stage, levels=c("Spermatogonia","Leptotene","Zygotene","Pachytene","Division II","Round Spermatid","Elongating \nSpermatid"))]
  
  return(curve_data)
}


#' Plot tSNE
#'
#' @param i ; either a numeric value corresponding to the component to plot, a gene name,
#' or the name of a variable stored in cell_data such as "group", "PseudoTime", or "library_size"
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param cell_metadata data.table with columns cell, Tsne1_QC1, Tsne2_QC2, and components V1, V2 etc.
#' @param expression_matrix cell by gene matrix of expression values
#' @param princurves list of outputs of princurve()
#' @param flip logical; should the cell scores be multiplied by -1 (score/gene loading sign combination are arbitrary like PCA scores).
#' @param predict logical; Should gene expression be the predicted (imputed) values or raw normalised
#' @param curve logical; Should principal curve (pseudotime) by plotted on top of the tSNE plot
#' @param stages logical; Should the stages be annotated
#' @param point_size numeric; size of point
#' @param principal_curve string; name of list entry in principal_curves object, to use for pseudotime line
#' @param dim1 string; name of variable in cell_metadata of the x axis to be used for the plot
#' @param dim2 string; name of variable in cell_metadata of the y axis to be used for the plot
#' @param log logical; should the varible mapped to the colour of the points be log transformed, default: F
#' @param curve_width numeric; width of the pseudotime curve
#'
#' @return The Y matrix from Rtsne output, rotated by angle
#' 
#' @details will add titles if you run load_component_orderings()
#' 
#' @export
#' 
#' @import ggplot2 data.table viridis ggnewscale

print_tsne <- function(i, factorisation=SDAresults, cell_metadata=cell_data, jitter=0, colourscale="diverging", expression_matrix=data, princurves=principal_curves, dim1="Tsne1_QC1", dim2="Tsne2_QC1", flip=FALSE, predict=FALSE, curve=FALSE, stages=FALSE, point_size=1, log=FALSE, principal_curve="df_9", curve_width=0.5){
  
  
  if(i %in% colnames(cell_metadata)){ # plot feature in cell_data
    
    tmp <- cell_metadata[,c(i,dim1, dim2),with=FALSE]
    names(tmp)[1] <- "feature"
    
    if(log){
      tmp$feature <- log(tmp$feature)
    }
    
    p <- ggplot(tmp[order(feature)], aes(get(dim1), get(dim2))) +
      geom_jitter(size=point_size, shape=21, stroke=0, aes(fill=feature), width=jitter, height=jitter) +
      ggtitle(paste(i))
    
    if(is.numeric(tmp$feature)){
      p <- p + scale_fill_viridis(guide = guide_colourbar(i), direction = -1)
    }else if(is.logical(tmp$feature)){
      p <- p + scale_fill_brewer(palette = "Set1")
    }else{
      p <- p + scale_fill_brewer(palette = "Paired") +
        guides(fill = guide_legend(override.aes = list(size=3, alpha=1), title = i))
    }
    
  }else if(mode(i)=="numeric"){# plot component
    
    tmp <- cell_metadata[,c(paste0("V",i), dim1, dim2),with=FALSE]
    names(tmp)[1] <- "score"
    
    if(flip){
      tmp[, score := score * (-1)]
    }
    
    p <- ggplot(tmp[order(abs(score))], aes(get(dim1), get(dim2))) +
      geom_point(size=point_size, stroke=0, aes(colour=score), position = position_jitter(width = jitter, height = jitter, seed=42L)) +
      ggtitle(paste(i,ifelse(exists("component_order_dt"),component_order_dt[component_number==i]$name,"")))
    
    
    if(colourscale=="diverging"){
      
      p <- p + scale_colour_gradientn(colours = log_colour_scale(range(tmp$score), scale=1, midpoint = "grey60", interpolate = "linear", asymetric = F),
                                      limits=c(-max(abs(tmp$score)),max(abs(tmp$score))),
                                      guide = guide_colourbar(paste0("Cell score\n(Component ",i,")"),
                                                              title.position = if(stages){"top"}else{"left"},
                                                              nbin=100))
    }else{
      
      if(tmp[,score][tmp[,which.max(abs(score))]] < 0 ){
        invert = 1
      }else{
        invert = -1
      }
      p <- p + scale_colour_viridis(direction = invert,
                                    guide = guide_colourbar(paste0("Cell score\n(Component ",i,")"),
                                                            title.position = if(stages){"top"}else{"left"}))
    }
    
    
  }else{ # plot gene not component
    gene = i

    if(predict){
      if(!gene %in% names(factorisation$loadings[[1]][1,])){
        return("Gene not found")
      }
      tmp <- merge(cell_metadata, sda_predict(gene, factorisation))[order(get(gene))]
    
    }else{
      if(!gene %in% colnames(expression_matrix)){
        return("Gene not found")
      }
      tmp <- merge(cell_metadata, expression_dt(gene, expression_matrix))[order(get(gene))]  
    }

    
    p <- ggplot(tmp, aes(get(dim1), get(dim2))) +
      geom_point(size=point_size, stroke=0, aes(colour=get(gene)), position = position_jitter(width = jitter, height = jitter, seed=42L)) +
      scale_color_gradient2(mid="lightgrey", high="midnightblue", low="lightgrey",
                            guide = guide_colourbar(paste0(gene," Expression"), title.position = if(stages){"top"}else{"left"})) +
      ggtitle(gene)
  }
  
  curve_data <- load_curve_data(princurves, principal_curve)
  colnames(curve_data)[1:2] <- c(dim1,dim2)
  
  if(curve){
    p <- p + new_scale_color() +
              geom_path(data = curve_data,
                       size = curve_width,
                       colour = 'black',
                       #aes(colour=Stage),
                       alpha=1,
                       arrow = arrow(angle = 12.5, ends = "first", type = "closed")) +
      scale_colour_brewer(palette = "Set1")
  }
  
  if(stages){
    
    p <- p + new_scale_color() +
            geom_path(data = curve_data[-c(1:290)],
              size = 4,
              aes(get(dim1)*1.7+3, get(dim2)*1.35-5, colour=Stage),
              alpha=0.5) +
      scale_colour_brewer(palette = "Set1")+
      guides(colour = guide_legend(ncol = 4, title.position = "top"))
  }
  
  p <- p + theme_minimal() +
    theme(legend.position = "bottom")
  
  if(dim1=="Tsne1_QC1"){
    p <- p + labs(x="t-SNE 1", y="t-SNE 2")
  }else{
    p <- p + labs(x=dim1, y=dim2)
  }
  
  return(p)

}




#' Plot tricolor tSNE
#'
#' @param df; string vector of two or three gene names
#' @param variables; string vector of two or three variable
#'  (gene) names present in df to create colours from
#' @param shift; logical, should the values be shifted so most 
#' negative becomes 0, for truely bi-signed variables (default: F)
#'
#' @return vector of colour rbg hex codes
#' 
#' @export
#' 
#' @import data.table
#' 
ternaryColours <- function(df, variables, shift=FALSE){
  df <- copy(df)
  for(i in variables){
    if(shift){
      df[, (i):= get(i) + abs(min(get(i)))] # shift so min value is 0  
    }
    df[get(i)<0, (i):= 0]
    df[, (i) := get(i) / max(unlist(list(1,get(i))))] # rescale to max of 1
  }
  return(
    rgb(red   = df[,variables[1], with=F][[1]],
        blue  = df[,variables[2], with=F][[1]],
        green = if(length(variables)==3){df[,variables[3], with=F][[1]]}else{0}
    )
  )
}


#' Plot tricolor tSNE
#'
#' @param genes; string vector of two or three gene names
#' @param returndf; logical, should the data be returned, or a plot (default)
#' @param predict; which factorisation to use for predicted expression
#' @param ptsize; numeric, size of the points in the plot
#' @param jitter; numeric, how jittered should be points be (to avoid overplotting)
#'
#' @return The Y matrix from Rtsne output, rotated by angle
#' 
#' @details print tsne with ternary colour scheme for expression
#' 
#' @export
#' 
#' @import ggplot2 data.table

tricolour_tsne <- function(genes, returndf=F, predict="SDA", cell_metadata=cell_data, ptsize=1.5, jitter=0.25){
  
  if(predict=="SDA"){
    cell_metadata <- merge(cell_metadata, sda_predict(genes))
  }else if(predict=="NNMF"){
    cell_metadata <- merge(cell_metadata, nmf_predict(genes))
  }else if(predict=="NNMF2"){
    cell_metadata <- merge(cell_metadata, nmf_predict(genes, nnmf_decomp5))
  }else{
    cell_metadata <- merge(cell_metadata, expression_dt(genes))
  }
  cell_metadata$rgb <- ternaryColours(cell_metadata, genes)
  
  if(returndf==T){
    return(cell_metadata)
  }else{
    
    set.seed(42) # for jitter consistency
    
    p <- ggplot(cell_metadata, aes(Tsne1_QC1, Tsne2_QC1)) +
      geom_jitter(aes(colour=rgb), stroke=0, size=ptsize, height=jitter, width=jitter) +
      scale_colour_identity() +
      annotate("text", -Inf, Inf, hjust = -0.5, vjust = 2, label = genes[1], colour="red", fontface="bold") +
      annotate("text", -Inf, Inf, hjust = -0.5, vjust = 4, label = genes[2], colour="blue", fontface="bold") +
      theme_dark()
    
    if(length(genes)==3){
      return(p + annotate("text", -Inf, Inf, hjust = -0.5, vjust = 6, label = genes[3], colour="green", fontface="bold"))
    }else{
      return(p)
    }
    
  }
}


#' Compute Imputed Expression values
#'
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param genes character vector; Gene Symbols of the genes to compute imputed expression for
#' @param name_extension string; will be appended to any gene names given (so resulting table could be joined with original values)
#'
#' @return A data.table with variables cell (the cell barcode) and a column for each gene containing the imputed expression values
#' 
#' @export
#' 
#' @import data.table

sda_predict <- function(genes, factorisation=SDAresults, name_extension=""){
  # use SDA parameters to create posterior prediction kind of
  stopifnot(!is.null(genes))
  
  predictions <- factorisation$scores %*% factorisation$loadings[[1]][, genes,drop=FALSE]
  predictions <- data.table(predictions, keep.rownames = T)
  setnames(predictions, c("cell",paste0(names(predictions)[-1],name_extension)))
  setkey(predictions, cell)
  return(predictions)
}

ica_predict <- function(genes, factorisation=nnmf_decomp, name_extension=""){
    # use ICA parameters to create posterior prediction
    stopifnot(!is.null(genes))
    
    predictions <- factorisation$S %*% factorisation$A[, genes,drop=FALSE]
    predictions <- data.table(predictions, keep.rownames = T)
    setnames(predictions, c("cell",paste0(names(predictions)[-1],name_extension)))
    setkey(predictions, cell)
    return(predictions)
  }

nmf_predict <- function(genes, factorisation=nnmf_decomp, name_extension=""){
  # use NNMF parameters to create posterior prediction
  stopifnot(!is.null(genes))
  
  predictions <- factorisation$W %*% factorisation$H[, genes,drop=FALSE]
  predictions <- data.table(predictions, keep.rownames = T)
  setnames(predictions, c("cell",paste0(names(predictions)[-1],name_extension)))
  setkey(predictions, cell)
  return(predictions)
}

pca_predict <- function(genes, factorisation=pcaresult, name_extension=""){
  # use PCA parameters to create posterior prediction
  stopifnot(!is.null(genes))
  
  predictions <- factorisation$projection %*% t(factorisation$loadings[genes,,drop=FALSE])
  predictions <- data.table(predictions, keep.rownames = T)
  setnames(predictions, c("cell",paste0(names(predictions)[-1],name_extension)))
  setkey(predictions, cell)
  return(predictions)
}

magic_predict <- function(gene, MAGIC_data, name_extension=""){
  predictions <- MAGIC_data$result[,gene]
  names(predictions) <- rownames(MAGIC_data$result)
  predictions <- as.data.table(predictions, keep.rownames = T)
  setnames(predictions,"rn","cell")
  setnames(predictions,"predictions",paste0(gene,name_extension))
}

#' Create data table of gene expression for a subset of genes
#'
#' @param expression_matrix cell by gene matrix of expression values
#' @param genes character vector; Gene Symbols of the genes to compute imputed expression for
#'
#' @details 
#' This functions requires the data matrix to be loaded
#' 
#' @return A data.table with variables cell (the cell barcode) and a column for each gene containing the raw normalised expression values
#' 
#' @export
#' 
#' @import data.table

expression_dt <- function(genes, expression_matrix=data){
  if(sum(!genes %in% colnames(expression_matrix))>0){return("Gene(s) not found")}
  expression <- data.table(as.matrix(expression_matrix[,genes,drop=FALSE]), keep.rownames = T)
  setnames(expression, "rn", "cell")
  setkey(expression, cell)
  return(expression)
}


# for plotting mutliple tsne plots, remove duplicated extra axes
simplify <- ggplot2::theme(legend.position = "none",
                           axis.title.x=element_blank(),
                           axis.title.y=element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks = element_blank())


#' Create list of grobs from a function and input
#'
#' @param fn; function which returns ggplot object
#' @param input; vector to be fed to the function
#'
#' @details 
#' Maybe this is redundant with lapply?
#' 
#' @export

create_grob_list <- function(fn=print_marker2, input=hist_subset){
  tmp <- list()
  
  for(i in seq_along(input)){
    tmp[[i]] <- fn(input[i]) + simplify
  }
  
  return(tmp)
}





#' Create long data table of gene expression for a subset of genes with pseudotimes
#'
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param cell_metadata data.table with columns cell, Tsne1_QC1, Tsne2_QC2, and components V1, V2 etc.
#' @param genes character vector; Gene Symbols of the genes to compute imputed expression for
#' 
#' @return A data.table with variables cell (the cell barcode), PseudoTime, Tsne coordinates,
#' value (containing imputed expression values), and Gene (a string indicator for which gene each value is for)
#' 
#' @export
#' 
#' @import data.table

gene_expression_pseudotime <- function(genes, factorisation=SDAresults, cell_metadata=cell_data){
  tmp <- merge(cell_metadata[somatic4==FALSE, c("cell","PseudoTime","Tsne1_QC1", "Tsne2_QC1","group"), with=FALSE],
               sda_predict(genes, factorisation, name_extension = ""))
  tmp <- melt(tmp, id.vars = c("cell","PseudoTime","Tsne1_QC1", "Tsne2_QC1","group"), variable.name = "Gene")
  tmp$Gene = factor(tmp$Gene, levels=genes)
  #tmp$variable = factor(tmp$variable, levels = sort(as.character(unique(tmp$variable))))
  return(tmp)
}





#' Plot pseudotime expression panel
#'
#' @param genes character vector; vector of gene symbols, of which expression values will be plotted over pseudotime
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param cell_metadata data.table with columns cell, Tsne1_QC1, Tsne2_QC2, and components V1, V2 etc.
#' @param ncol numeric; number of columns to use for faceting
#' @param title string; Title to be used in plot
#' @param gam_k numeric; passed to gam formula, controls smoothing
#' @param point_size numeric; size of each cell (point) in plot
#' @param highlight_reigon numeric vector; vector of size two with min and max pseudotime values between which to highlight in the plot
#'
#' @return ggplot2 object
#'
#' @export
#' @import ggplot2
plot_pseudotime_expression_panel <- function(genes, factorisation=SDAresults, cell_metadata=cell_data, ncol=7, title="Histone Genes", gam_k=5, point_size=0.2, highlight_reigon=NULL){
  
  if(mode(genes)=="character"){
    data <- gene_expression_pseudotime(genes, factorisation, cell_metadata)
  }
  
  p <-  ggplot(data, aes(-PseudoTime, value, colour=Gene)) +
      geom_point(alpha=0.2, size=point_size) +
      geom_smooth(method = "gam", formula = y ~ s(x, k = 50), se=FALSE, colour="black") +
      ylab("Predicted Gene Expression") + xlab("Pseudotime") +
      ggtitle(title) +
      facet_wrap(~Gene, scales="free_y", ncol = ncol) +
      theme(legend.position = "none")
  
  if(!is.null(highlight_reigon)){
    p <- p + geom_rect(data=data.frame(xmin=highlight_reigon[1], xmax=highlight_reigon[2], ymin=-Inf, ymax=Inf),
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                fill="grey20",
                colour=NA,
                alpha=0.2,
                inherit.aes = FALSE)
  }

  return(p)
  
}





#' Create long table combining gene expression and cell metadata
#'
#' @param genes character vector; vector of gene symbols, of which expression values will be plotted over pseudotime
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param cell_metadata data.table with columns cell, Tsne1_QC1, Tsne2_QC2, and components V1, V2 etc.
#' @param predict logical; if TRUE (default FALSE) imputed/predicted gene expression is used
#' @param expression_matrix matrix; rows=cells, columns=genes, unimputed but normalised (scaled) expression values
#'
#' @return data.table with columns 'cell', 'PsuedoTime', 'Tsne1_QC1', 'Tsne2_QC1', 'Gene', and 'Expression'
#'
#' @export
#' @import data.table
#' 
melt_genes <- function(genes, cell_metadata=cell_data, expression_matrix=data, predict=FALSE, factorisation=SDAresults){
  # generate melted version of cell_data on the fly for subset of genes
  # rather than keeping two 7Gb data.tables in memory
  #saveRDS(merge_sda3_melt2, "PseudoTime_Shiny/data_V3.rds")
  
  genes <- genes[genes %in% colnames(expression_matrix)]
  
  if(!all(genes %in% colnames(expression_matrix))){
    return("Gene not found")
  }
  
  if(predict){
    tmp <- merge(cell_metadata, sda_predict(genes, factorisation))
  }else{
    tmp <- merge(cell_metadata, expression_dt(genes, expression_matrix))
  }

  id_colums <- c("cell", "PseudoTime", "Tsne1_QC1", "Tsne2_QC1")
  
  merge_sda3_melt2 <- melt(tmp[,c(id_colums, genes), with=FALSE], id.vars = id_colums, value.name="Expression", variable.name="Gene")
  merge_sda3_melt2[,cell:=as.factor(cell)]
  return(merge_sda3_melt2)
}





#' Find the top x genes for a component
#'
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param component integer or string; number or name of component in 'SDAresults' object
#' @param n integer; number of gene symbols to retrieve
#' @param values logical; if TRUE (default FALSE) named gene loadings are returned
#'
#'
#' @return charachter vector of gene symbols, or named gene loadings are returned if values=TRUE
#'
#' @export
#'
get_top_genes <- function(component, factorisation=SDAresults, n=20, values=FALSE){
  
  # first decide which side we want
  highest_cell <- names(sort(-abs(factorisation$scores[,component]))[1])
  direction <- sign(factorisation$scores[highest_cell,component])
  
  top <- names(head(sort(-direction * factorisation$loadings[[1]][component, ]), n))
  
  if(values){
    top <- factorisation$loadings[[1]][component, top]
  }
  return(top)
}





#' Plot gene loadings with cell scores underneath
#'
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param gene_locations ouput of SDAtools::load_gene_locations()
#' @param i integer or string; number or name of component in 'SDAresults' object
#' @param max.items integer; number of genes to label in loadings plot
#'
#' @return ggplot2 object
#'
#' @export
#' 
#' @import ggplot2 gridExtra
print_loadings_scores <- function(i, factorisation=SDAresults, gene_locations=gene_annotations, max.items=30){
  grid.arrange(grobs=list(genome_loadings(factorisation$loadings[[1]][i,], label_both = FALSE, max.items = max.items, gene_locations = gene_locations) + ggtitle(i),
                          ggplot(data.table(cell_index=1:nrow(factorisation$scores), score=factorisation$scores[,paste0("V",i)],
                                            experiment=gsub("_.*","",gsub("[A-Z]+\\.","",rownames(factorisation$scores)))),
                                 aes(cell_index, score, colour=experiment)) +
                            geom_point(size=0.5, stroke=0) + xlab("Cell Index (by Library Size)") +
                            ylab("Score") +
                            scale_color_brewer(palette = "Paired") +
                            guides(color=guide_legend(ncol=2, override.aes = list(size=2)))),
               nrow=2, heights = c(5,2))
}





#' Print list of genes in component and their annotations
#'
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param i integer or string; number or name of component in 'SDAresults' object
#' @param n integer; number of genes to include in the list, default=100
#' @param annotations data.table with columns gene_symbol, fold_enrichment_fantom, fold_enrichment_gtex,
#'  Human_Orthologue, Inferitility_Gene, chromosome_name and start_position, created in 24_Gene-annotations.Rmd
#' @param order; order of ouput list by loading, "absolute", "decreasing" or "increasing"
#' 
#' @details See also get_top_genes() for a simpler function
#'
#' @return print of data.table
#'
#' @export
#' 
#' @import data.table
gene_list <- function(i, n=100, factorisation=SDAresults, annotations=gene_annotations, order="absolute") {
  tmp <- data.table(as.matrix(factorisation$loadings[[1]][i,]), keep.rownames = TRUE)
  setnames(tmp, c("gene_symbol","Loading"))
  setkey(tmp, gene_symbol)
  
  tmp <- merge(tmp, annotations, all.x = TRUE)
  tmp$Testis_Enriched <- tmp$fold_enrichment_gtex > 2 | tmp$fold_enrichment_fantom > 2
  
  setcolorder(tmp, c("gene_symbol","Human_Orthologue", "Loading", "Infertility_Gene", "Testis_Enriched", "fold_enrichment_gtex",  "fold_enrichment_fantom"))
  
  if(order=="absolute"){
    return(tmp[order(-abs(Loading))][1:n])
  }else if(order=="decreasing"){
    return(tmp[order(-Loading)][1:n])
  }else{
    return(tmp[order(Loading)][1:n])
  }
}




#' Volcano plot of gene ontology enrichment resultd
#'
#' @param data R object; GO results
#' @param component string; name of component e.g. "V5N" (N meaning negative size loadings)
#' @param top_n integer; how many GO terms to label
#' @param OR_threshold numeric; threshold on odds ratio of GO terms to label (applied after top_n)
#' @param label_size numeric; size of GO term labels
#'
#' @return ggplot2 object
#'
#' @export
#' 
#' @import ggplot2 data.table
go_volcano_plot <- function(data=GO_data, component="V5N", top_n=30, label_size=3, OR_threshold=1){
  
  if(class(data)[1]=="data.table"){
    tmp <- data[Component==component]
  }else{
    tmp <- data.table(data[[component]])
  }
  
  return(
    ggplot(tmp, aes(GeneOdds/BgOdds, -log(pvalue,10), size=Count)) +
      geom_point(aes(colour=p.adjust<0.05)) +
      scale_size_area() +
      geom_label_repel(data = tmp[order(p.adjust)][1:top_n][p.adjust<0.7][GeneOdds/BgOdds > OR_threshold],
                       aes(label = Description),
                       size = label_size,
                       force=5) + 
      ggtitle(paste0(component," Volcano plot for Gene Ontology (Biological Processes) enrichment analysis")) +
      xlab("Odds Ratio") +
      scale_x_log10(limits=c(1,NA), breaks=c(1,2,5,10,15,20))
  )
}





#' Plot gene expression through pseudotime (Raw & Imputed)
#'
#' @param genes character vector; Gene Symbols of the genes to compute imputed expression for
#' @param cell_metadata data.table with columns cell, Tsne1_QC1, Tsne2_QC2, and components V1, V2 etc.
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param expression_matrix matrix; rows=cells, columns=genes, unimputed but normalised (scaled) expression values
#' @param facet logical; if TRUE (default) each column is a faceted plot,
#' if FALSE just two plots are returned imputed & raw
#' 
#' @return A ggplot2 object
#' 
#' @export
#' 
#' @import ggplot2 data.table cowplot

imputed_vs_raw <- function(genes, cell_metadata=cell_data, factorisation=SDAresults, expression_matrix=data, facet=T){
  library(scales)
  
  tmp <- merge(cell_metadata[somatic4==FALSE], expression_dt(genes, expression_matrix))[,c(genes,"cell","PseudoTime","Tsne1_QC1", "Tsne2_QC1"), with=FALSE]
  tmp$Type = "Raw (Normalised)"
  tmp <- melt(tmp, id.vars = c("cell","PseudoTime","Tsne1_QC1", "Tsne2_QC1", "Type"))
  
  predicted_genes <- sda_predict(genes, factorisation)
  predicted_genes$Type = "Imputed"
  predicted_genes <- merge(cell_metadata[somatic4==FALSE,c("cell","PseudoTime","Tsne1_QC1","Tsne2_QC1"), with=FALSE], predicted_genes)
  predicted_genes <- melt(predicted_genes, id.vars = c("cell","PseudoTime","Tsne1_QC1", "Tsne2_QC1", "Type"))
  
  #scale_y_continuous(trans = 'asinh', breaks=c(0,2,5,10)) +
  #scale_y_continuous(trans = 'asinh', breaks=c(0,1,2,5,10)) +
  
  # Two columns - raw vs imputed
  if(facet==TRUE){
    return(
      plot_grid(
        ggplot(tmp, aes(-PseudoTime, value)) +
          scale_color_brewer(palette = "Set1") +
          geom_point(alpha=0.5, size=0.7, shape=20, stroke=0, colour=RColorBrewer::brewer.pal(3, name="Set1")[1]) +
          ylab("Normalised Gene Expression") + xlab("Pseudotime") +
          facet_wrap(~variable, scales = "free_y", ncol=1, strip.position = "left") +
          ggtitle("Unimputed") +
          theme_minimal() +
          theme(strip.background = element_blank(), strip.text.y = element_text(size=12), strip.placement = "outside")
        ,
        ggplot(predicted_genes, aes(-PseudoTime, value)) +
          scale_color_brewer(palette = "Set1") +
          geom_point(alpha=0.5, size=0.7, shape=20, stroke=0, colour=RColorBrewer::brewer.pal(3, name="Set1")[2]) +
          xlab("Pseudotime") + ylab("") +
          facet_wrap(~variable, scales = "free_y", ncol=1, strip.position = "right") +
          ggtitle("SDA Imputed") +
          theme_minimal()  + theme(strip.text.y = element_blank())
      , ncol=2, rel_widths = c(6,5))
    )
  }else{
    
    # print(
    # ggplot(rbind(tmp,predicted_genes), aes(-PseudoTime, asinh(value), colour=Type)) +
    #   scale_color_brewer(palette = "Set1") +
    #   geom_point(alpha=0.3, size=0.1) +
    #   ylab("Gene Expression") + xlab("Pseudotime") +
    #   facet_grid(variable~Type, scales = "free_y") +
    #   ggtitle("Raw vs Imputed Gene Expression") +
    #   theme_minimal() +
    #   theme(legend.position = "none")
    # )
    
    # two rows (genes overplotted) - raw vs imputed
    return(
      ggplot(rbind(tmp, predicted_genes), aes(-PseudoTime, asinh(value), colour=Type)) +
        geom_point(alpha=0.5, size=0.3) +
        ylab("Gene Expression") + xlab("Pseudotime") +
        facet_wrap(~Type, scales = "free_y", ncol = 1) +
        ggtitle("Raw vs Imputed Gene Expression") +
        scale_color_brewer(palette = "Set1") +
        theme_minimal()
    )
  }
}

asinh_trans <- function() {
  scales::trans_new("asinh",
            transform = asinh,
            inverse   = sinh)
}





#' Find pairs of cells that are close in tSNE space
#'
#' @param test_groups string; type of cell to pair cell to
#' @param cell_subset character vector; cell barcodes of cells to consider pairing, e.g. only cells with a certian component score >2
#' @param cell_metadata data.table with columns cell, Tsne1_QC1, Tsne2_QC2, and components V1, V2 etc.
#' @param swap_groups logical; control cells and test_group cells are swapped.
#' As this is a greedy algorithm, it can be better to swap the groups when there are few control cells
#'
#' @details 
#' For a given cell type, find the closest cells of a different type. e.g. find pairs of wt and Hormad1 mutatnt cells close together on tSNE.
#' This is useful for performing differentiall expression analysis whilst accounting for developmental stage of cells
#' The default control cells are genetically WT cells (WT,mj,SPG,SPD,SPCII,SPCI)
#' 
#' @return
#' 
#' @export
#' 
#' @import data.table

match_cells2 <- function(test_groups=NULL, cell_subset=NULL, cell_metadata=cell_data, swap_groups=F){
  
  if(!is.null(cell_subset)){
    cell_metadata <- cell_metadata[cell %in% cell_subset]
    print(nrow(cell_metadata))
  }
  
  control_groups <- c("WT","mj","SPG","SPD","SPCII","SPCI")
  
  if(is.null(test_groups)){
    test_groups <- c("Cul4a","Hormad1","CNP","Mlh3")
  }
  
  if(swap_groups){
    a <- control_groups
    b <- test_groups
    test_groups <- a
    control_groups <- b
  }
  
  test_cells <- cell_metadata[group %in% test_groups]$cell
  
  matched_cells <- matrix(nrow=length(test_cells), ncol=3)
  
  matched_cells[,1] <- test_cells
  
  colnames(matched_cells) <- c("Mutant","WT","distance")
  
  for(i in seq_along(test_cells)){
    
    Tsne1_cell_coord <- cell_metadata[cell==test_cells[i]]$Tsne1_QC1
    Tsne2_cell_coord <- cell_metadata[cell==test_cells[i]]$Tsne2_QC1
    
    tmp <- cell_metadata[Tsne1_QC1 < (Tsne1_cell_coord + 3) & Tsne1_QC1 > (Tsne1_cell_coord - 3)][Tsne2_QC1 < (Tsne2_cell_coord + 3) & Tsne2_QC1 > (Tsne2_cell_coord - 3)]
    
    # calculate distance to WT cells
    distances <- rowSums(cbind(Tsne1_cell_coord - tmp[group %in% control_groups]$Tsne1_QC1,
                               Tsne2_cell_coord - tmp[group %in% control_groups]$Tsne2_QC1)^2)
    
    names(distances) <- tmp[group %in% control_groups]$cell
    
    # remove WT cells already assigned neighbours
    distances <- distances[which(!tmp[group %in% control_groups]$cell %in% matched_cells[,2])]
    
    # assign closest cell is avaliable
    if (sum(!is.na(distances))==0){
      matched_cells[i,3] <- matched_cells[i,2] <- NA
    }else if(distances[which.min(distances)] > 3){
      matched_cells[i,3] <- matched_cells[i,2] <- NA
    }else{
      matched_cells[i,3] <- distances[which.min(distances)]
      matched_cells[i,2] <- names(which.min(distances))
    }
  }
  
  return(matched_cells)
}





#' Calculate enrichment for a set of genes in a single component
#'
#' @param component number or name of component
#' @param genes character vector of gene set to look for enrichment for
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param pos logical; if TRUE (default) the positive gene loadings are ranked highest
#' @param threshold numeric; how many genes should be included in the top genes list used to calculate enrichment
#' @param bg_genes which genes to use as the 'background' default uses all genes in the SDAresults object from SDA
#' @param test string, either binomial (from Exact::exact.test, slow) or fisher (R's fisher.test)
#'
#' @details 
#' Calculate p value of enrichemnt (fishers test) in a component, given a set of genes
#' see component_enrichment for a function to run on all components
#' 
#' @return
#' 
#' @export

single_component_enrichment <- function(gene_vector=SDAresults$loadings[[1]][1,], genes = Mybl1_genes,  pos=T, threshold=200, bg_genes=NULL, test="fisher"){
  
  if(is.null(bg_genes)){
    bg_genes <- names(gene_vector)
  }
  
  ranks <- sort(rank((if(pos){-1}else{1}) * gene_vector[bg_genes])[genes], na.last=T)
  
  #names(tmp) <- bg_genes
  
  yes <- sum(ranks <= threshold)
  
  contingency_table <- matrix(c(yes,
                                length(genes) - yes,
                                threshold - yes,
                                length(bg_genes) - length(genes) - (threshold - yes)),
                              nrow=2, ncol=2,
                              dimnames = list(ReachedTreshold = c(T,F),
                                              InList = c(T, F)))
  
  if(test=='binomial'){
    # this test doesn't assume all marginals are fixed, but is much slower
    test <- Exact::exact.test(contingency_table, alternative = "greater", model="binomial", cond.row = T, to.plot=F, npNumbers = 5)
  }else{
    test <- fisher.test(contingency_table, alternative="greater")
  }
  
  return(list("fisher_test"=test,
              "ranks"=ranks,
              "counts"=c(yes,
                length(genes),
                length(bg_genes),
                pos,
                threshold),
              "Positive"=pos,
              "contingency_table"=contingency_table
              ))
  
}






#' Calculate enrichment for a set of genes for each component
#'
#' @param genes character vector of gene set to look for enrichment for
#' @param threshold numeric; how many genes should be included in the top genes list used to calculate enrichment
#' @param factorisation SDA factorisation object, output of SDAtools::load_results()
#' @param bg_genes which genes to use as the 'background' default uses all genes in the SDAresults object from SDA
#' @param test string, either binomial (from Exact::exact.test, slow) or fisher (R's fisher.test)
#'
#' @details 
#' for each component (P&N seperately), calculate p value enrichemnt (fishers test) for a set of genes
#' This is a multi-component wrapper to AZ_pvalue
#' requires component_order_all object and component_order_dt to be loaded into the environment using load_component_orderings()
#' 
#' see single_component_enrichment for a function to run on a single gene vector/component
#' 
#' @return A data.table with a row for each component (+ve & -ve seperately) with columns giving the p.value of enrichment,
#' Odds Ratio, Component number and name, component type (Somatic or Meiotic), and component order
#' 
#' @export
#' 
#' @import data.table

component_enrichment <- function(genes, threshold=500, loadings_matrix=SDAresults$loadings[[1]], orderSDA=T, bg_genes=NULL, test="fisher"){
  
  stopifnot(nrow(loadings_matrix) < ncol(loadings_matrix))
  
  n <- nrow(loadings_matrix)
  
  pos_ad_e <- lapply(1:n, function(x) single_component_enrichment(loadings_matrix[x,], genes = genes, pos = T, threshold = threshold, bg_genes=bg_genes, test=test))
  neg_ad_e <- lapply(1:n, function(x) single_component_enrichment(loadings_matrix[x,], genes = genes, pos = F, threshold = threshold, bg_genes=bg_genes, test=test))
  
  combined <- c(pos_ad_e,neg_ad_e)
  names(combined) <- c(paste0(1:n,"P"),paste0(1:n,"N"))
  
  tmp <- data.table(
    component = names(combined),
    p.value = sapply(combined, function(x) x$fisher_test$p.value),
    OR = sapply(combined, function(x) x$fisher_test$estimate),
    hits = sapply(combined, function(x) x$counts[1]),
    hit_names = I(sapply(combined, function(x) names(x$ranks[x$ranks<x$counts[5]])))
  )
  
  #t(sapply(combined, function(x) c("p.value"=x$fisher_test$p.value, "OR"=x$fisher_test$estimate[[1]], "hits"=x$counts[1]))),
  
  if(orderSDA){
    
    tmp[, component := factor(component, levels=component[c(rbind(component_order_all,component_order_all+50))])]
    
    tmp$name <- rep(component_order_dt$name,2)
    
    tmp[order(component), rank := 1:100]
    tmp[rank >= 41, Type := "Meiotic"]
    tmp[rank < 41, Type := "Somatic"]
  }else{
    tmp$name <- tmp$component
    tmp$Type <- "Unk"
  }
  
  return(tmp)
}

#' Plot Manhatten Style Enrichment for each component
#'
#' @param enrichments result of the component_enrichment() function
#' @param topn int; how many components should be labeled
#' @param repel_force; degree of force to move points apart form each other
#' @param legend_position; numeric vector of length 2, giving positions
#' between 0 and 1 where the legend should be placed  within the plot
#'
#' @details 
#' 
#' @return A ggplot object
#' 
#' @export
#' 
#' @import ggplot2 SDAtools

manhatten_plot <- function(enrichments, topn=1, repel_force=1, legend_position=c(0.85,0.8)){
  ggplot(enrichments, aes(component, p.value, label=paste0(component," - ",name))) +
    geom_point(aes(colour=Type, size=OR)) +
    scale_size() +
    geom_label_repel(data=enrichments[order(p.value)][1:topn], force=repel_force, min.segment.length = 0) +
    geom_hline(yintercept = 0.05/100, col="red", size=0.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5),
          legend.position = legend_position,
          legend.box = "horizontal",
          legend.background = element_rect(colour="black")) +
    scale_color_brewer(palette = "Set1") +
    xlab("Components by Type and Pseudotime") +
    scale_y_continuous(trans=SDAtools::reverselog_trans())
}



#' Calculate log colour scale
#' 
#' Instead of taking log of the values and plotting that,
#' log the colour scale instead
#' 
#' @param values numeric vector, could be matrix or vector of values you will plot
#' used to calculate range for asymetrical 
#' @param scale use this to rescale the range, to change what part of the log scale you're on
#' @param midpoint colour of central value (0)
#' @param interpolate passed to colorRampPalette
#' @param asymetric logical; should 0 be at the center of offset to account for range of values
#' @param print logical; should colour scale be printed or should vector of hex values be returned
#' 

log_colour_scale <- function(values=vvv3, scale=0.05, midpoint="white", interpolate="spline", asymetric=T, print=F){
  newcols <- brewer.pal(11, "RdYlBu")
  newcols[6] <- midpoint
  #newcols <- colorRampPalette(c(newcols[1],newcols[3],newcols[6],newcols[9],newcols[11]), interpolate="spline", space="Lab")(2000)
  n=1000
  newcols <- colorRampPalette(newcols, interpolate=interpolate, space="Lab")(2*n)
  
  range = max(abs(values)) * scale
  logcolindex1 <- floor((log(seq(1, range, length.out = n))/log(range))*n)
  
  # index2 corrected for range being skewed to positive  values
  if(asymetric){
    index2length = floor(abs(min(values)/max(values))*n)
  }else{
    index2length = n
  }
  
  logcolindex2 <- floor((log(seq(1, range, length.out = index2length)) / log(range)) * n) #n should be index2length if not using full range of colours
  
  indexes <- c(abs(logcolindex1-n)[n:1],logcolindex2+n)
  newcolscale = newcols[indexes][(length(indexes)-1):1]
  
  if(print){
    # visualise scales
    cols <- function(a) image(1:length(a), 1, as.matrix(1:length(a)), col=a, axes=T , xlab="", ylab="")
    return(cols(newcolscale))
  }else{
    return(newcolscale)
  }
}


#' Calculate 2D density
#'
#' @param x numeric vector;
#' @param y numeric vector;
#' @param n integer; number of grid points in each direction
#' @param kern numeric; size of kernal to use for smoothing
#'
#' @details for colouring points in plot by the density of points, from http://slowkow.com/notes/ggplot2-color-by-density/
#' 
#' @return numeric vector of density at each point
#' 
#' @export
#' 
#' @importFrom MASS kde2d
get_density <- function(x, y, n = 1000, kern=0.02){
  dens <- MASS::kde2d(x = x, y = y, n = n, h=c(kern,kern))
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



#' Vectorised multionimial distribution
#'
#' @param counts numeric matrix; matrix of K columns of integers in 0:size. each sample is a row
#' @param prob numeric non-negative matrix with K columns, specifying the probability for the K classes; is internally normalized to sum 1. Infinite and missing values are not allowed. Each sample is a row.
#'
#' @details vectorised version of R base dmultinom
#' 
#' @return LOG probability density / likelihood of the data counts given the parameters in prob
#' 
#' @export
#' 
dmultinom_v <- function(counts, prob){
  lgamma(rowSums(counts)+1) + rowSums(counts * log(prob/rowSums(prob)) - lgamma(counts+1))
}


#' Reverse DGE Normalisation
#'
#' @param dge numeric matrix; rows are cells, columns are genes
#' @param sign logical; if TRUE (default) the sign of the values is retained (negatives are still negative)
#' @param cells character vector; strings corresponding to cell names to include
#' @param lib_correction numeric vector; size correction factor used for each cell
#' @param standard_dev numeric vector; standard deviations used to scale each gene
#'
#' @details Reverses the transformation of dropsim::normaliseDGE() to get back to a counts scale
#' 
#' @return Matrix of cells by genes, with values on the count scale
#' 
#' @export
#' 
reverse_normalisation <- function(dge, sign=TRUE, cells=cell_subset, lib_correction=lib_size_correction_factor, standard_dev=sds){
  if(sign){
    tmp <- t(((t(dge)*standard_dev)^2 * sign(t(dge))) / 10000) * lib_correction[cells]  
  }else{
    tmp <- t((t(dge)*standard_dev)^2 / 10000) * lib_correction[cells]
  }
  
  if(class(tmp)!="matrix"){
    names(tmp@x) <- NULL  
  }
  
  return(tmp)
}



#' Sort matrix to highest correlations are along the diagonal
#'
#' @param mat numeric matrix to be sorted
#' @param order1 order of rows
#' 
#' @details For a given row order, sorts the columns so that each row
#' has the highest correlation (in a greedy fashion)
#' 
#' @return numeric matrix
#' 
#' @export
#' 
sort_matrix <- function(mat=cross_cor, order1=max_cor_per_compoent){
  mix_comps <- rownames(mat)
  
  mix_order <- character()
  
  for (i in seq_along(order1)){
    tmp2 <- names(which.max(abs(mat[mix_comps, names(order1)[i] ])))
    if(is.null(tmp2)){
      tmp2 <- mix_comps
    }
    mix_order <- c(mix_order, tmp2)
    mix_comps <- mix_comps[!mix_comps %in% tmp2]
  }
  
  
  gene_loading_correlation = t(mat[mix_order, names(order1)][nrow(mat):1,]) #names(max_cor_per_compoent)
  
  return(gene_loading_correlation)
}


#' Plot correlation heatmap comparing two sets of loadings
#'
#' @param f1 gene loadings from factorisation method 1
#' @param f2 gene loadings from factorisation method 1
#' @param names string; names of f1 and f2 used in heatmap
#' @param method correlation method, note kendall not allowed
#' @param return_cor logical; should the ordered correlation matrix be returned,
#' if FALSE (default) heatmap is plotted using ComplexHeatmap package
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ComplexHeatmap
#' 
compare_factorisations <- function(f1=nnmf_decomp$H, f2=SDAresults$loadings[[1]], names=c("f1","f2"), method="spearman", return="hm", split=c(F,F), inject_matrix=NULL, randomise=FALSE){
  
  if(method=="kendall"){
    stop("kendall is too slow, not allowed it will crash R")
  }
  
  
  # check if whole SDAresults is passed, or just loadings
  if(is.list(f1)){
    f1loadings <- f1$loadings[[1]]
    f1scores <- f1$scores
  }else{
    f1loadings <- f1
  }
  
  if(is.list(f2)){
    f2loadings <- f2$loadings[[1]]
    f2scores <- f2$scores
  }else{
    f2loadings <- f2
  }
  
  # subset gene loadings to each has the same set of genes
  if(ncol(f1loadings)>ncol(f2loadings)){
    common_genes <- colnames(f1loadings)[colnames(f1loadings) %in% colnames(f2loadings)]
    f1loadings <- f1loadings[,common_genes]
    f2loadings <- f2loadings[,common_genes]
  }else{
    common_genes <- colnames(f2loadings)[colnames(f2loadings) %in% colnames(f1loadings)]
    f1loadings <- f1loadings[,common_genes]
    f2loadings <- f2loadings[,common_genes]
  }
  
  if(randomise){
    f2loadings <- f2loadings[,sample(1:ncol(f2loadings))]
  }

  # split positive and negative loadings
  if(split[1]){
    f1loadings <- t(split_loadings(f1loadings))
    if(is.list(f1)){
      f1scores <- split_loadings(t(f1scores))
    }
  }
  if(split[2]){
    f2loadings <- t(split_loadings(f2loadings))
    if(is.list(f2)){
      f2scores <- split_loadings(t(f2scores))
    }
  }
  
  
  cross_cor <- abs(cor(t(f1loadings), t(f2loadings),
                       method = method))
  
  #pt_order <- colnames(cross_cor)[component_order_all[50:1]]
  
  if(!is.null(inject_matrix)){
    cross_cor <- inject_matrix
  }
  
  max_cor_per_compoent <- sort(apply(cross_cor, 2, max), decreasing = T)
  
  gene_loading_correlation = sort_matrix(cross_cor, max_cor_per_compoent)
  
  max_cor_per_compoent_cols <- apply(gene_loading_correlation, 2, max)
  
  names(dimnames(gene_loading_correlation)) <- names[2:1]
  
  # create side annotation for heatmap
  make_annotation <- function(SDAscores, max_cor_per_compoent, side="row", return=return){
    
    tmp2 <- data.table(sum_abs_cell_score = colMeans(abs(SDAscores))[names(max_cor_per_compoent)],
                       cells_in_component = colSums(abs(SDAscores)>1)[names(max_cor_per_compoent)],
                       score_sum_2 = colMeans(SDAscores^2)[names(max_cor_per_compoent)],
                       score_sum_med_2 = apply(SDAscores^2, 2, function(x) quantile(x, probs=0.98))[names(max_cor_per_compoent)],
                       max_cor_per_compoent,
                       percent_nonWT_cells = 1 - (apply(abs(SDAscores), 2, function(x) sum(grepl("WT|mj|SPG|SPD|SPCII|SPCI",names(which(x>1))))) / 
                                                   apply(abs(SDAscores), 2, function(x) length(which(x>1))))[names(max_cor_per_compoent)],
                       mutant_contribution = 1 - (apply(abs(SDAscores), 2, function(x) sum(x[grepl("WT|mj|SPG|SPD|SPCII|SPCI",names(x))])) / 
                                                    apply(abs(SDAscores), 2, function(x) sum(x)))[names(max_cor_per_compoent)],
                       max_cell_score = apply(abs(SDAscores), 2, max)[names(max_cor_per_compoent)],
                       component=names(max_cor_per_compoent))
    
    if(return=="hm"){
      ha = HeatmapAnnotation(df = tmp2[,.(sum_abs_cell_score, max_cell_score)],
                             col = list(sum_abs_cell_score=circlize::colorRamp2(c(min(tmp2$sum_abs_cell_score), max(tmp2$sum_abs_cell_score)), c("white", "red")),
                                        max_cell_score = circlize::colorRamp2(c(0, max(max(tmp2$max_cell_score),0.1)), c("white", "blue"))),
                                        #mutant_contribution = circlize::colorRamp2(c(min(tmp2$mutant_contribution), max(max(tmp2$mutant_contribution),0.1)), c("white", "black"))),
                             annotation_legend_param=list(legend_direction = "horizontal",
                                                          legend_width = unit(5, "cm"), title_position = "lefttop"),
                             which = side,
                             show_legend = if(side=="row"){TRUE}else{FALSE}
      )
      
      return(ha)
    }else{
      return(tmp2)
    }
    
  }
  
  if(return=="cor"){ # return cross correlation matrix
    return(gene_loading_correlation)
    
  }else if(return=="hm"){ # return heatmap
    hm = Heatmap(gene_loading_correlation,
                 col = viridisLite::viridis(100),
                 cluster_rows = F,
                 cluster_columns = F,
                 heatmap_legend_param = list(legend_direction = "horizontal", title_position = "lefttop",legend_width = unit(5, "cm")),
                 row_title = paste("Components from",names[2]),
                 column_title = paste("Components from",names[1]),
                 name = paste("Gene loading",method,"correlation"),
                 row_names_max_width = unit(10, "npc"),
                 bottom_annotation = if(is.list(f1)){make_annotation(f1scores, max_cor_per_compoent_cols, "column", return)})
    
    if(is.list(f2)){
      return(draw(hm + make_annotation(f2scores, max_cor_per_compoent, "row", return), heatmap_legend_side = "bottom", annotation_legend_side = "bottom"))
    }else{
      return(draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom"))
    }
    
  }else{ # reutrn annotation df
    return(rbind(make_annotation(f1scores, max_cor_per_compoent_cols, "row", return),
                 make_annotation(f2scores, max_cor_per_compoent, "row", return)))
  }
  
}


#' Rotate SDA factorisation
#'
#' @param X SDA results list to be rotated
#' @param reference SDA results list to for X to be rotated towards
#' 
#' @return SDA results list, but with the gene loadings and scores matrices rotated by procrustes
#' 
#' @export
#' @import vegan
#' 
rotate_SDA <- function(X=results9_WT, reference=results1){
  
  if(ncol(X$loadings[[1]])>ncol(reference$loadings[[1]])){
    common_genes <- colnames(X$loadings[[1]])[colnames(X$loadings[[1]]) %in% colnames(reference$loadings[[1]])]
  }else{
    common_genes <- colnames(reference$loadings[[1]])[colnames(reference$loadings[[1]]) %in% colnames(X$loadings[[1]])]
  }
  
  rot_9WT <- vegan::procrustes(t(reference$loadings[[1]][,common_genes]), 
                               t(X$loadings[[1]][,common_genes]))
  
  colnames(rot_9WT$Yrot) <- paste0(rownames(X$loadings[[1]]),"rot")
  
  colnames(rot_9WT$rotation) <- rownames(reference$loadings[[1]])
  rownames(rot_9WT$rotation) <- rownames(X$loadings[[1]])
  
  # rotate scores too
  rot_9WT_scoresrot <- X$scores %*% rot_9WT$rotation
  colnames(rot_9WT_scoresrot) <- paste0(rownames(X$loadings[[1]]),"rot")
  str(rot_9WT_scoresrot)
  
  row_9WT_list <- list(loadings=list(t(rot_9WT$Yrot)), scores=rot_9WT_scoresrot, rotation=rot_9WT$rotation)  
  
  return(row_9WT_list)
}



#' Plot ROC curve for imputation rankings
#'
#' @param i integer or string; number or name of cell in the matrices
#' @param mt named list of  gene rank by cell matricies, cumulative sum of predicted for each cell
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ggplot2
#' 
plotCellAUC <- function(i, mt){
  
  aucg <- data.table(frac_genes = seq(1,length(mt[[1]][,i]),1)/length(mt[[1]][,i]),
                     sapply(mt, function(x) x[,i]))
  
  aucg <- melt(aucg, id.vars = "frac_genes", variable.name = "Method", value.name = "Cumulative Test data reads")
  
  return(
    ggplot(aucg, aes(frac_genes, `Cumulative Test data reads`, colour=Method)) + 
      geom_line() + 
      xlab("Fraction of Genes (Ranked high to low)") + 
      theme_minimal() +
      theme(legend.position = "bottom") +
      scale_color_brewer(palette = "Set1") #+
      # annotate("label",
      #          label = paste0("Imputed AUC: ",signif(sum(cumsum_predict[,i])/ length(cumsum_mean[,i]),3),
      #                         "\nCellwise AUC: ",signif(sum(cumsum_train[,i])/ length(cumsum_mean[,i]),3),
      #                         "\nAverage AUC: ",signif(sum(cumsum_mean[,i])/ length(cumsum_mean[,i]),3)),
      #          x = 0.6, y = 0.25)
  )
}

#' Save R Objects, as multiple files in a folder
#'
#' Like save but instead of one rds file, creates one file per object
#' This gives much more flexibility, to add, modify, rename, remove and
#' share individual objects.
#' 
#' See also load2
#'
#' @param folder string; path of folder in which to save objects
#' @param list character vector; names of objects to be saved.
#' Defaults to all objects in global environment (not functions).
#' @param compress logical; specifying whether saving to a named file is to use "gzip" compression,
#' or one of "gzip", "bzip2" or "xz" to indicate the type of compression to be used.
#' 
#' @export
#' 
save2 <- function(folder="", list = NULL, compress=FALSE){
  if(is.null(list)){
    list <- setdiff(ls(.GlobalEnv), lsf.str(.GlobalEnv))
  }
  
  for(item in list){
    saveRDS(get(item), paste0(folder,item,".rds"), compress = compress)
  }
}

#' Save R Objects, as multiple files in a folder
#'
#' Like load but instead of one rds file, loads all rds objects in a folder.
#' This gives much more flexibility, to add, modify, rename, remove and
#' share individual objects.
#' See also save2
#'
#' @param folder string; path of folder containing objects to be loaded
#' @param pattern an optional regular expression. 
#' Only files which match the regular expression will be loaded.
#' 
#' @export
#' 
load2 <- function(folder, pattern = "\\.rds$"){
  files <- list.files(path = folder, pattern = pattern, full.names = T)
  for(item in files){
    assign(gsub(".rds$","",basename(item)), readRDS(item), envir = .GlobalEnv)
  }
}

# deprecated
# e.g. a=sample(1:20000,1);cellAUC_old(a);print(a)
cellAUC_old=function(i){
  predvec=predicted_train_counts_0[i,]
  trainvec=raw_data_train[i,]
  testvec=raw_data_test[i,]
  set.seed(42)
  reorder=sample(1:length(predvec))
  predvec=predvec[reorder]
  trainvec=trainvec[reorder]
  testvec=testvec[reorder]
  order1=order(predvec,decreasing=T)
  order2=order(trainvec,decreasing=T)
  
  c4=cumsum(predvec[order1]/sum(predvec))
  c1=cumsum(testvec[order1])/sum(testvec)
  c2=cumsum(testvec[order2])/sum(testvec)
  c3=cumsum(testvec[av_exp_order])/sum(testvec) #order(reorder)
  frac1=seq(1,length(testvec),1)/length(testvec)
  
  plot(frac1,c1,type="l",xlim=c(0,1),ylim=c(0,1))
  lines(frac1,c2,type="l",col=2)
  lines(frac1,c3,lty="dotted",col="grey")
  lines(frac1,c4,lty="dotted",col="blue")
  
  
  fracs=c(sum(c1),sum(c2),sum(c3))/length(c1)
  return(fracs)
}