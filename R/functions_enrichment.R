

#' Perform Gene Ontology (Biological Process) enrichment analysis of components
#'
#' @param component string; name of component e.g. "V5N" (N meaning negative size loadings)
#' @param geneNumber integer; number of top genes to include
#' @param side char; should negative of positive values be considered the top
#' @param database OrgDb; object of class OrgDb from the package AnnotationDbi
#'
#' @return data.frame with each row a GO term, columns Description, pvalue, Enrichment, etc
#'
#' @export
#' 
#' @import clusterProfiler
GO_enrichment <- function(component, geneNumber = 250, side="N", factorisation=results, database=Musculus){
  
  if(side=="N"){
    top_genes <- data.table(as.matrix(factorisation$loadings[[1]][component, ]), keep.rownames = TRUE)[order(V1)][1:geneNumber]$rn
  }else{
    top_genes <- data.table(as.matrix(factorisation$loadings[[1]][component, ]), keep.rownames = TRUE)[order(-V1)][1:geneNumber]$rn
  }
  
  gene_universe <- data.table(as.matrix(factorisation$loadings[[1]][component,]), keep.rownames = TRUE)$rn
  
  ego <- enrichGO(gene = top_genes,
                  universe = gene_universe,
                  OrgDb = database,
                  keyType = 'SYMBOL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text=x)))
  
  ego@result$Enrichment <- frac_to_numeric(ego@result$GeneRatio)/frac_to_numeric(ego@result$BgRatio)
  ego@result$GeneOdds <- unlist(lapply(strsplit(ego@result$GeneRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
  ego@result$BgOdds <- unlist(lapply(strsplit(ego@result$BgRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
  
  return(ego@result)
  
}



motifGenes <- function(motif){
  motifRankings[[1]]@rankings$rn[order(motifRankings[[1]]@rankings[,motif,with=FALSE])]
  
  # old way  
  # motif_list <- motifRankings$`500bp`@rankings[,c(motif,"rn"), with=FALSE]
  # setnames(motif_list, c("rank","rn"))
  # motif_list <- motif_list[order(rank)]$rn
  
}


join_lists <- function(motif_list, test_list, joined=T, values="test"){
  
  test_list <- test_list[test_list %in% motif_list]
  motif_list <- motif_list[motif_list %in% test_list]
  motif_list <- motif_list[!duplicated(motif_list)]
  
  if(joined){
    if(values=="test"){
      component_rankings = seq_along(test_list)
      names(component_rankings) = test_list
      return(component_rankings[motif_list]) # test_list ranking, for each motif ranking ordered by motif rakning
    }else if(values=="motif"){
      component_rankings = seq_along(motif_list)
      names(component_rankings) = motif_list
      return(component_rankings[test_list]) # motif ranking, for each test_list ranking ordered by test_list rakning
    }
  }else{
    return(list(test_list, motif_list))
  }
  
}



plot_enrichment <- function(data, random=trace_random, xlim=2000){
  
  if(ncol(data)>3){
    rows=2
  }else{
    rows=1
  }
  
  data <- cbind(index = seq_along(random), data)
  data <- melt(data, id.vars = "index")
  data <- data[,.(index, enrichment = cumsum(value-trace_random)), by=variable]
  
  print(
    ggplot(data, aes(index, enrichment, colour=variable)) +
      geom_path() +
      facet_zoom(xy = index < xlim, zoom.size = 1, horizontal = F) +
      theme(legend.position = "bottom") +
      guides(col = guide_legend(nrow = rows))
  )
}

# QQ plot function from the internet
ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}



visualise_correlation <- function(motif_genes, component, motif_name="motif X"){
  joined_list <- join_lists(motif_genes, names(sort(results$loadings[[1]][component,])), values = "motif")
  joined_list2 <- join_lists(motif_genes, names(sort(results$loadings[[1]][component,])), values = "test")
  
  print(qplot(1:length(joined_list), joined_list, size=I(2), shape=I(21)) +
          labs(x="Genes Ranked by SDA Component Loading",
               y="Genes Ranked by iRegulon",
               title=paste("Component", component,"vs. motif:",motif_name)) +
          theme_classic())
  
  print(qplot(1:100, tapply(joined_list, cut(seq_along(joined_list), 100), mean)) +
          theme_classic() +
          labs(x="Bin of SDA Component Loading",
               y="Mean iRegulon Rank",
               title=paste("Component", component,"vs. motif:",motif_name)))
  
  print(qplot(1:100, tapply(joined_list2, cut(seq_along(joined_list2), 100), mean)) +
          theme_classic() +
          labs(x="Bin of iRegulon Rank",
               y="Mean SDA Component Loading",
               title=paste("Component", component,"vs. motif:",motif_name)))
}


# calculate cor.test for each motif for a single component
correlation_test <- function(component, positive=FALSE){
  # check first column is rn
  stopifnot(colnames(motifRankings$`500bp`@rankings)[1]=="rn")
  
  # create placeholder matrix
  correlation_statistics <- matrix(nrow=ncol(motifRankings$`500bp`@rankings)-1, ncol=2)
  rownames(correlation_statistics) <- colnames(motifRankings$`500bp`@rankings)[2:ncol(motifRankings$`500bp`@rankings)]
  
  # calc correlation for each motif
  for (motif in rownames(correlation_statistics)){
    
    if(positive){
      joined_ranking <- join_lists(motifGenes(motif), names(sort(results$loadings[[1]][component,][results$loadings[[1]][component,]>0], T)), values = "motif")
    }else{
      joined_ranking <- join_lists(motifGenes(motif), names(sort(results$loadings[[1]][component,][results$loadings[[1]][component,]<0])), values = "motif")
    }
    
    correlation_statistics[motif,] <- unlist(cor.test(joined_ranking, seq_along(joined_ranking))[c("p.value","estimate")])
  }
  
  return(correlation_statistics)
  
}

# calculate cor.test for each component for a single motif
correlation_test_motif <- function(motif){
  # check first column is rn
  stopifnot(colnames(motifRankings$`500bp`@rankings)[1]=="rn")
  
  # create placeholder matrix
  correlation_statistics <- matrix(nrow=100, ncol=2)
  rownames(correlation_statistics) <- paste0(rep(1:50, each=2), c("P","N"))
  
  # calc correlation for each motif
  for (component in 1:50){
    
    # positive
    joined_ranking <- join_lists(motifGenes(motif), names(sort(results$loadings[[1]][component,][results$loadings[[1]][component,]>0], T)), values = "motif")
    correlation_statistics[seq(1,99,2)[component],] <- unlist(cor.test(joined_ranking, seq_along(joined_ranking))[c("p.value","estimate")])
    
    # negative
    joined_ranking <- join_lists(motifGenes(motif), names(sort(results$loadings[[1]][component,][results$loadings[[1]][component,]<0])), values = "motif")
    correlation_statistics[seq(1,99,2)[component]+1,] <- unlist(cor.test(joined_ranking, seq_along(joined_ranking))[c("p.value","estimate")])
    
  }
  
  return(correlation_statistics)
  
}

# split one side of loadings into spike and slab, then test difference in mean
two_groups_test <- function(component, positive=FALSE, pip_threshold=0.5){
  # check first column is rn
  stopifnot(colnames(motifRankings$`500bp`@rankings)[1]=="rn")
  
  # create placeholder matrix
  correlation_statistics <- matrix(nrow=ncol(motifRankings$`500bp`@rankings)-1, ncol=3)
  rownames(correlation_statistics) <- colnames(motifRankings$`500bp`@rankings)[2:ncol(motifRankings$`500bp`@rankings)]
  
  # calc correlation for each motif
  for (motif in rownames(correlation_statistics)){
    
    if(positive){
      joined_ranking <- join_lists(motif, names(sort(results$loadings[[1]][component,][results$loadings[[1]][component,]>0], T)), values = "motif")
    }else{
      joined_ranking <- join_lists(motif, names(sort(results$loadings[[1]][component,][results$loadings[[1]][component,]<0])), values = "motif")
    }
    
    in_spike <- names(joined_ranking) %in% colnames(results$pips[[1]])[results$pips[[1]][component,] < pip_threshold]
    joined_ranking_spike <- joined_ranking[in_spike]
    joined_ranking_slab <- joined_ranking[!in_spike]
    
    correlation_statistics[motif,] <- c(wilcox.test(joined_ranking_spike, joined_ranking_slab, alternative="g")$p.value, mean(joined_ranking_slab), mean(joined_ranking_spike))
  }
  return(correlation_statistics)
}


# enrichment test for the top x genes (faster than all if only looking at the top)
enrichment_test2 <- function(test_list, topx=NULL){
  # check first column is rn
  stopifnot(colnames(motifRankings$`500bp`@rankings)[1]=="rn")
  
  # create placeholder matrix
  correlation_statistics <- matrix(nrow=ncol(motifRankings$`500bp`@rankings)-1, ncol=3)
  rownames(correlation_statistics) <- colnames(motifRankings$`500bp`@rankings)[2:ncol(motifRankings$`500bp`@rankings)]
  
  # calc correlation for each motif
  for (motif in rownames(correlation_statistics)){
    joined_ranking <- join_lists(motif, test_list, F) # NB joining the lists takes longer than the get trace
    
    if(is.null(topx)){topx = joined_ranking}
    
    trace <- get_trace3b(joined_ranking[[1]][1:topx], joined_ranking[[2]][1:topx])
    
    correlation_statistics[motif,] <- c(sum(cumsum(trace-trace_random[1:topx])[1:250]), sum(cumsum(trace-trace_random[1:topx])[1:500]), sum(cumsum(trace-trace_random[1:topx])[1:1000]))
  }
  
  return(correlation_statistics)
  
}



# If one list shorter:

# for (i in seq_along(ranked_list)){
#   trace[i] <- sum(duplicated(c(ranked_list[1:i], soh)))
# }


# Naive approach
get_trace <- function(ranked_list, test_list){
  
  lists <- cbind(ranked_list, test_list)
  
  counter <- numeric(length=nrow(lists))
  
  for (row in 1:nrow(lists)){
    counter[row] <- sum(duplicated(c(lists[1:row,])))
  }
  
  return(counter)
}

# Slightly better, only calculate for each additional row
get_trace2 <- function(ranked_list, test_list){
  
  lists <- cbind(ranked_list, test_list)
  
  counter <- numeric(length=nrow(lists))
  
  for (row in 1:nrow(lists)){
    counter[row] <- lists[row,1] %in% lists[1:row,2] + lists[row,2] %in% lists[1:row,1]
  }
  
  return(counter)
}

# Hashtable (data.table) version, many times faster
get_trace3b <- function(ranked_list, test_list){
  
  ranked_list <- data.table(rank1=seq_along(ranked_list),gene=ranked_list)
  test_list <- data.table(rank2=seq_along(test_list),gene=test_list)
  
  setkey(ranked_list, gene)
  setkey(test_list, gene)
  
  lists <- merge(ranked_list, test_list)
  
  # find position where duplicate occours
  lists[, dup_pos := max(rank1, rank2), by=gene]
  
  # tally duplicates per position
  table_lists <- table(lists$dup_pos)
  table_lists <- cbind(as.numeric(rownames(table_lists)), table_lists)
  
  # create duplicate count per position
  #lists[, trace := 0]
  trace <- rep(0,nrow(ranked_list))
  trace[table_lists[,1]] <- table_lists[,2]
  
  return(trace)
  
}

# Simon's code

# aa=(cumsum(joined_ranking_homer_myb))
# plot(diff(aa[seq(1,14000,14)]))
# plot(diff(aa[seq(1,14000,140)])/100)
# 
# plot(cumsum(sort(joined_ranking_homer_myb)[trace_42]))
# plot(cumsum(sort(joined_ranking_homer_myb)[sample(trace_42)]))

# sm.density(tt[motif_list], h=c(5,5), nbins=361,
# xlim=c(-180,180), ylim=c(-180,180),
# xlab="phi", ylab="psi", zlab="")