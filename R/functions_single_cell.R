
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
print_raw_tsne <- function(gene, expression_matrix=data, cell_metadata=datat){
  ggplot(merge(cell_metadata, expression_dt(gene, expression_matrix)), aes(Tsne1, Tsne2, color=get(gene))) +
    geom_point(size=0.2) +
    scale_color_viridis(direction=-1) +
    theme(legend.position = "bottom") +
    ggtitle(paste0(gene," - t-SNE")) + simplify
}




#' Print PCA
#'
#' @param cell_metadata data.table with columns Tsne1, Tsne2, and PCs
#' @param pc string; name of principal compopnent (column of datat), of which cell score values will be used to colour points (cells)
#'
#' @return ggplot2 object
#'
#' @export
#' @import ggplot2
print_pca <- function(cell_metadata=datat, pc){
  ggplot(cell_metadata, aes(Tsne1, Tsne2, color=get(pc))) +
    geom_point(size=0.2) +
    scale_colour_distiller(palette="YlOrRd", direction=1) +
    theme(legend.position = "bottom") +
    ggtitle(paste0(pc," - t-SNE")) + simplify
}




#' Clustered heatmap (partially deprecated)
#'
#' @param cell_metadata data.table; datat, with subset of gene columns
#' @param name string; title
#' @param annotation.col data.table; subset of datat, passed to annCol of aheatmap
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
clustered_heatmap <- function(cell_metadata=datat, name, annotation.col=NULL, colv_order=NULL, col_lab=NULL, row_text_size=0.5){
  
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

# depreceated, use print_tsne(preditc=T) instead
# print_marker2 <- function(gene){
#   gene <- paste0(gene,"")
#   ggplot(tmp[order(get(gene))], aes(Tsne1_QC1, Tsne2_QC1, color=get(gene))) +
#     geom_point(size=0.05) +
#     scale_color_viridis(direction=-1) +
#     ggtitle(gene)
# }

# depreceated, use print_tsne(preditc=T) instead
# print_marker2 <- function(gene){
#   gene <- paste0(gene,"_predict")
#   ggplot(tmp[order(get(gene))], aes(Tsne1_QC1, Tsne2_QC1, color=get(gene))) +
#     geom_point(size=0.05) +
#     scale_color_viridis(direction=-1) +
#     ggtitle(gene)
# }

# depreceated, use print_tsne() instead
# print_marker <- function(gene){
#   
#   if(!gene %in% colnames(data)){
#     return("Gene not found")
#   }
#   
#   tmp <- merge(datat, expression_dt(gene))[order(get(gene))]
#   
#   ggplot(tmp, aes(Tsne1_QC1, Tsne2_QC1, color=get(gene))) +
#     geom_point(size=0.05) +
#     scale_color_viridis(direction=-1) +
#     ggtitle(gene)
# }


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

plot_cell_scores <- function(component="V25", point_size=0.6, cell_metadata=datat){
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



#' Plot tSNE
#'
#' @param i ; either a numeric value corresponding to the component to plot, a gene name,
#' or the name of a variable stored in datat such as "group", "PseudoTime", or "library_size"
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
#'
#' @return The Y matrix from Rtsne output, rotated by angle
#' 
#' @details will add titles if you run load_component_orderings()
#' 
#' @export
#' 
#' @import ggplot2 data.table viridis ggnewscale

print_tsne <- function(i, factorisation=results, cell_metadata=datat, expression_matrix=data, princurves=principal_curves, dim1="Tsne1_QC1", dim2="Tsne2_QC1", flip=FALSE, predict=FALSE, curve=FALSE, stages=FALSE, point_size=1, log=FALSE, principal_curve="df_9", curve_width=0.5){
  
  if(i %in% colnames(cell_metadata)){ # plot feature in datat
    
    tmp <- cell_metadata[,c(i,dim1, dim2),with=FALSE]
    names(tmp)[1] <- "feature"
    
    if(log){
      tmp$feature <- log(tmp$feature)
    }
    
    p <- ggplot(tmp[order(feature)], aes(get(dim1), get(dim2))) +
      geom_point(size=point_size, shape=21, stroke=0, aes(fill=feature)) +
      ggtitle(paste(i))
    
    if(is.numeric(tmp$feature)){
      p <- p + scale_fill_viridis(guide = guide_colourbar(i))
    }else{
      p <- p + scale_fill_brewer(palette = "Paired") +
        guides(fill = guide_legend(override.aes = list(size=3, alpha=1), title = i))
    }
    
  }else if(mode(i)=="numeric"){# plot component not gene, could do this for any cell attribute in datat - future TODO
    
    tmp <- cell_metadata[,c(paste0("V",i), dim1, dim2),with=FALSE]
    names(tmp)[1] <- "score"
    
    if(flip){
      tmp[, score := score * (-1)]
    }
    
    if(tmp[,score][tmp[,which.max(abs(score))]] < 0 ){
      invert = 1
    }else{
      invert = -1
    }
    
    p <- ggplot(tmp[order(score)], aes(get(dim1), get(dim2))) +
      geom_point(size=point_size, stroke=0, aes(colour=score)) +
      scale_colour_viridis(direction = invert, guide = guide_colourbar(paste0("Cell score\n(Component ",i,")"), title.position = if(stages){"top"}else{"left"})) +
      ggtitle(paste(i,ifelse(exists("component_order_dt"),component_order_dt[component_number==i]$name,"")))
    
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
      geom_point(size=point_size, stroke=0, aes(colour=get(gene))) +
      scale_colour_viridis(direction=-1, guide = guide_colourbar(paste0(gene," Expression"), title.position = if(stages){"top"}else{"left"})) +
      ggtitle(gene)
  }
  
  curve_data <- data.table(princurves[[principal_curve]]$s[princurves[[principal_curve]]$tag, 1],
                           princurves[[principal_curve]]$s[princurves[[principal_curve]]$tag, 2])
  
  colnames(curve_data) <- c(dim1,dim2)
  
  if(curve){
    p <- p + new_scale_color() +
              geom_path(data = curve_data,
                       size = curve_width,
                       #aes(colour=Stage),
                       alpha=0.7,
                       arrow = arrow(angle = 12.5, ends = "first", type = "closed")) +
      scale_colour_brewer(palette = "Set1")
  }
  
  if(stages){
    m <- nrow(curve_data)
    
    curve_data[m:(m*0.975),Stage := "Spermatogonia"]
    curve_data[(m*0.975):(m*0.96),Stage := "Leptotene"]
    curve_data[(m*0.96):(m*0.93),Stage := "Zygotene"]
    curve_data[(m*0.93):(m*0.62),Stage := "Pachytene"]
    curve_data[(m*0.62):(m*0.5),Stage := "Division II"]
    curve_data[(m*0.5):(m*0.2),Stage := "Round Spermatid"]
    curve_data[(m*0.2):0,Stage := "Elongating \nSpermatid"]
    
    curve_data[,Stage := factor(Stage, levels=c("Spermatogonia","Leptotene","Zygotene","Pachytene","Division II","Round Spermatid","Elongating \nSpermatid"))]
    
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
    p <- p + labs(x="Umap 1", y="Umap 2")
  }
  
  return(p)

}


# for plotting mutliple tsne plots, remove duplicated extra axes
simplify <- ggplot2::theme(legend.position = "none",
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank())



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

sda_predict <- function(genes, factorisation=results, name_extension=""){
  # use SDA parameters to create posterior prediction kind of
  
  predictions <- factorisation$scores %*% factorisation$loadings[[1]][, genes,drop=FALSE]
  predictions <- data.table(predictions, keep.rownames = T)
  setnames(predictions, c("cell",paste0(names(predictions)[-1],name_extension)))
  setkey(predictions, cell)
  return(predictions)
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





#' Create list of grobs from a function and input
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

gene_expression_pseudotime <- function(genes, factorisation=results, cell_metadata=datat){
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
plot_pseudotime_expression_panel <- function(genes, factorisation=results, cell_metadata=datat, ncol=7, title="Histone Genes", gam_k=5, point_size=0.2, highlight_reigon=NULL){
  
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
#'
#' @return data.table with columns 'cell', 'PsuedoTime', 'Tsne1_QC1', 'Tsne2_QC1', 'Gene', and 'Expression'
#'
#' @export
#' @import data.table
#' 
melt_genes <- function(genes, cell_metadata=datat, expression_matrix=data, predict=FALSE){
  # generate melted version of datat on the fly for subset of genes
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
#' @param component integer or string; number or name of component in 'results' object
#' @param n integer; number of gene symbols to retrieve
#' @param values logical; if TRUE (default FALSE) named gene loadings are returned
#'
#'
#' @return charachter vector of gene symbols, or named gene loadings are returned if values=TRUE
#'
#' @export
#'
get_top_genes <- function(component, factorisation=results, n=20, values=FALSE){
  
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
#' @param rna_locations ouput of SDAtools::load_gene_locations()
#' @param i integer or string; number or name of component in 'results' object
#' @param max.items integer; number of genes to label in loadings plot
#'
#' @return ggplot2 object
#'
#' @export
#' 
#' @import ggplot2 gridExtra
print_loadings_scores <- function(i, factorisation=results, gene_locations=rna_locations, max.items=30){
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
#' @param i integer or string; number or name of component in 'results' object
#' @param n integer; number of genes to include in the list, default=100
#' @param fantom data.table with columns Gene.Name and fold_difference from FANTOM dataset
#' @param gtex data.table with column Gene.Name, Human_Orthologue, and fold_difference_gtex
#' 
#'
#' @return print of data.table
#'
#' @export
#' 
#' @import data.table
print_gene_list <- function(i, n=100, factorisation=results, infertility=infertility_genes, fantom=fantom_subset_summary, gtex=gtex_summary) {
  tmp <- data.table(as.matrix(factorisation$loadings[[1]][i,]), keep.rownames = TRUE)[order(-abs(V1))][1:n]
  setnames(tmp, c("Gene.Name","Loading"))
  setkey(tmp, Gene.Name)
  
  # Add is testis enriched annotations
  tmp <- merge(tmp, fantom_summary_subset, all.x = TRUE)
  tmp <- merge(tmp, gtex_summary[,.(Gene.Name, Human_Orthologue, fold_difference_gtex = signif(fold_difference,3))], all.x = TRUE)
  tmp$Testis_Enriched <- gsub("FALSE","",tmp$fold_difference_gtex > 2 | tmp$fold_difference > 2)
  
  # Add infertility gene annotations
  tmp$Infertility_Gene <- gsub("FALSE","",tmp$Gene.Name %in% infertility)
  
  setcolorder(tmp, c("Gene.Name","Human_Orthologue", "Loading", "Infertility_Gene", "Testis_Enriched", "fold_difference_gtex",  "fold_difference"))
  
  # Display Result
  print(tmp[order(-abs(Loading))])
}




#' Volcano plot of gene ontology enrichment resultd
#'
#' @param x R object; GO results
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
#' @param facet logical; if TRUE (default) each column is a faceted plot,
#' if FALSE just two plots are returned imputed & raw
#' 
#' @return A ggplot2 object
#' 
#' @export
#' 
#' @import ggplot2 data.table cowplot

imputed_vs_raw <- function(genes_tmp, cell_metadata=datat, factorisation=results, expression_matrix=data, facet=T){
  library(scales)
  
  tmp <- merge(cell_metadata[somatic4==FALSE], expression_dt(genes_tmp, expression_matrix))[,c(genes_tmp,"cell","PseudoTime","Tsne1_QC1", "Tsne2_QC1"), with=FALSE]
  tmp$Type = "Raw (Normalised)"
  tmp <- melt(tmp, id.vars = c("cell","PseudoTime","Tsne1_QC1", "Tsne2_QC1", "Type"))
  
  predicted_genes <- sda_predict(genes_tmp, factorisation)
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
          ggtitle("Imputed") +
          theme_minimal()  + theme(strip.text.y = element_blank())
      , ncol=2)
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

match_cells2 <- function(test_groups=NULL, cell_subset=NULL, cell_metadata=datat, swap_groups=F){
  
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
#' @param bg_genes which genes to use as the 'background' default uses all genes in the results object from SDA
#' @param test string, either binomial (from Exact::exact.test, slow) or fisher (R's fisher.test)
#'
#' @details 
#' Calculate p value of enrichemnt (fishers test) in a component, given a set of genes
#' 
#' @return
#' 
#' @export

single_component_enrichment <- function(component, genes = Mybl1_genes, factorisation=results, pos=T, threshold=200, bg_genes=NULL, test="fisher"){
  
  if(is.null(bg_genes)){
    bg_genes <- names(factorisation$loadings[[1]][1,])
  }
  
  ranks <- sort(rank((if(pos){-1}else{1}) * factorisation$loadings[[1]][component, bg_genes])[genes], na.last=T)
  
  #names(tmp) <- bg_genes
  
  yes <- sum(ranks < threshold)
  
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
#' @param bg_genes which genes to use as the 'background' default uses all genes in the results object from SDA
#'
#' @details 
#' for each component (P&N seperately), calculate p value enrichemnt (fishers test) for a set of genes
#' This is a multi-component wrapper to AZ_pvalue
#' requires the results object from SDA to be loaded as well as the component_order_all object from 'R/functions_ordering.R'
#' 
#' @return A data.table with a row for each component (+ve & -ve seperately) with columns giving the p.value of enrichment,
#' Odds Ratio, Component number and name, component type (Somatic or Meiotic), and component order
#' 
#' @export
#' 
#' @import data.table

component_enrichment <- function(genes, threshold=500, factorisation=results, bg_genes=NULL, test="fisher"){
  pos_ad_e <- lapply(1:50, function(x) single_component_enrichment(x, genes = genes, factorisation=factorisation, pos = T, threshold = threshold, bg_genes=bg_genes, test=test))
  neg_ad_e <- lapply(1:50, function(x) single_component_enrichment(x, genes = genes, factorisation=factorisation, pos = F, threshold = threshold, bg_genes=bg_genes, test=test))
  
  combined <- c(pos_ad_e,neg_ad_e)
  names(combined) <- c(paste0(1:50,"P"),paste0(1:50,"N"))
  
  tmp <- data.table(
    component = factor(names(combined), levels=names(combined)[c(rbind(component_order_all,component_order_all+50))]),
    p.value = sapply(combined, function(x) x$fisher_test$p.value),
    OR = sapply(combined, function(x) x$fisher_test$estimate),
    hits = sapply(combined, function(x) x$counts[1]),
    hit_names = I(sapply(combined, function(x) names(x$ranks[x$ranks<x$counts[5]])))
  )
  
  #t(sapply(combined, function(x) c("p.value"=x$fisher_test$p.value, "OR"=x$fisher_test$estimate[[1]], "hits"=x$counts[1]))),
  
  tmp$name <- rep(component_order_dt$name,2)
  
  tmp[order(component), rank := 1:100]
  tmp[rank >= 41, Type := "Meiotic"]
  tmp[rank < 41, Type := "Somatic"]
  
  return(tmp)
}

#' Plot Manhatten Style Enrichment for each component
#'
#' @param enrichments result of the component_enrichment() function
#' @param topn int; how many components should be labeled
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
    geom_point(aes(colour=Type)) +
    geom_label_repel(data=enrichments[order(p.value)][1:topn], force=repel_force, min.segment.length = 0) +
    geom_hline(yintercept = 0.05/100, col="red", size=0.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5),
          legend.position = legend_position,
          legend.background = element_rect(colour="black")) +
    scale_color_brewer(palette = "Set1") +
    xlab("Components by Type and Pseudotime") +
    scale_y_continuous(trans=SDAtools::reverselog_trans())
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





#' Plot ROC curve for imputation rankings
#'
#' @param i integer or string; number or name of cell
#' @param cumsum_predict gene rank by cell matrix, cumulative sum of predicted for each cell
#' @param cumsum_train gene rank by cell matrix, cumulative sum of train for each cell
#' @param cumsum_mean mean
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ggplot2
#' 
plotCellAUC <- function(i, cumsum_predict, cumsum_train, cumsum_mean){
  
  aucg <- data.table(frac_genes = seq(1,length(cumsum_predict[,i]),1)/length(cumsum_predict[,i]),
                     Imputed = cumsum_predict[,i],
                     Training = cumsum_train[,i],
                     Average = cumsum_mean[,i])
  
  aucg <- melt(aucg, id.vars = "frac_genes", variable.name = "Source of Expression Ranking", value.name = "Cumulative Test data reads")
  
  return(
    ggplot(aucg, aes(frac_genes, `Cumulative Test data reads`, colour=`Source of Expression Ranking`)) + 
      geom_line() + 
      xlab("Fraction of Genes (Ranked high to low)") + 
      theme_minimal() +
      theme(legend.position = "bottom") +
      scale_color_brewer(palette = "Set1") +
      annotate("label",
               label = paste0("Imputed AUC: ",signif(sum(cumsum_predict[,i])/ length(cumsum_mean[,i]),3),
                              "\nTraining AUC: ",signif(sum(cumsum_train[,i])/ length(cumsum_mean[,i]),3),
                              "\nAverage AUC: ",signif(sum(cumsum_mean[,i])/ length(cumsum_mean[,i]),3)),
               x = 0.6, y = 0.25)
  )
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