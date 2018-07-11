
#' Download sequences flanking gene
#'
#' @details 
#' This function is deprecated, use bedtools getfasta instead
#'
#' @export
#' 
#' @import biomaRt

download_flank_sequences <- function(gene_names, length, flank_type="gene_flank"){
  # Download Sequences
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = "mmusculus_gene_ensembl")
  #data.table(listAttributes(ensembl))[grepl(pattern = "Flank",x = description)]
  #coding_gene_flank
  
  flank_sequences <- getBM(attributes = c("external_gene_name",flank_type),
                           filter = c("external_gene_name","upstream_flank"),
                           values = list(gene_names,length),
                           mart = ensembl, uniqueRows=TRUE, checkFilters = FALSE)
  
  return(flank_sequences)
}




#' Download location of Transcripion Start Sites from Ensembl Biomart
#'
#' @param gene_names charachter vector; vector of gene names (symbols) to download TSS for
#'
#' @details 
#' Uses ensembl biomart to download location of Transcripion Start Sites
#' If there are multiple tss the best is determined by tsl, then by appris, then by length.
#' 
#' @return
#' 
#' @export
#' 
#' @import biomaRt data.table

download_tss <- function(gene_names){
  # Download Sequences
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = "mmusculus_gene_ensembl")
  #data.table(listAttributes(ensembl))[grepl(pattern = "Flank",x = description)]
  
  tss <- getBM(attributes = c("external_gene_name","strand","transcript_tsl","transcription_start_site","transcript_appris","transcript_length","ensembl_transcript_id","chromosome_name"),
               filter = c("external_gene_name","transcript_biotype"),
               values = list(gene_names,"protein_coding"),
               mart = ensembl, uniqueRows=TRUE, checkFilters = FALSE)
  
  tss <- data.table(tss)
  
  # tidy up names
  #tss[,.N,by=transcript_tsl][order(-N)]
  tss[, transcript_tsl := gsub(" \\(assigned to previous version [0-9]+\\)","",transcript_tsl)]
  tss[, transcript_tsl := as.numeric(gsub("tsl","",tss$transcript_tsl))]
  
  #tss[,.N,by=transcript_appris][order(-N)]
  tss[, transcript_appris := gsub("alternative1","6",transcript_appris)]
  tss[, transcript_appris := gsub("alternative2","7",transcript_appris)]
  tss[, transcript_appris := as.numeric(gsub("principal","",tss$transcript_appris))]
  
  # order such that best first
  setorder(tss, transcript_tsl, transcript_appris, transcript_length, na.last = TRUE)
  
  # remove lower ranked duplicate transcripts
  tss <- tss[!duplicated(tss, by="external_gene_name")]
  
  return(tss)
}




#' Export sequences to file in FASTA format
#'
#' @param sequences data.frame; data.table / data.frame
#'
#' @details 
#' Deprecated
#' 
#' @return
#' 
#' @export

export_FASTA <- function(sequences, file){
  sequences$external_gene_name <- paste0(">",sequences$external_gene_name)
  D <- do.call(rbind, lapply(seq(nrow(sequences)), function(i) t(sequences[i, ])))
  write.table(D, file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}




#' Find mode of values
#'
#' @param x numeric vector; values to find Mode of
#'
#' @details 
#' Finds most common value in vector
#' 
#' @return numeric value
#' 
#' @export

Mode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}




#' Create matrix representation of sequences with highlighted motif
#'
#' @param motif output of MotifFinder
#'
#' @details 
#' From the output of MotifFinder, creates matrix of sequences, with the motif found highlighted using lowercase nucleotide letters
#' 
#' @return Charachter matrix (one sequence per row, one column per nucleotide) with rownames
#' 
#' @export
#' @import stringr

sequence_matrix <- function(motif){
  
  # remove sequences with unusual lengths
  #odd_lengths <- which(nchar(motif$seqs) != Mode(nchar(motif$seqs)))
  
  seq_names <- names(motif$seqs)
  
  mean(nchar(motif$seqs))
  
  motif$seqs <- str_pad(motif$seqs, Mode(nchar(motif$seqs)), "left",pad="N")
  mean(nchar(motif$seqs))
  
  # motif$seqs <- motif$seqs[-odd_lengths]
  # 
  # reg_odd <- which(motif$whichregs %in% odd_lengths)
  # str(reg_odd)
  # str(motif$whichregs)
  # motif$whichregs <- motif$whichregs[-reg_odd]
  # str(motif$whichregs
  # str(motif$whichpos)
  # motif$whichpos <- motif$whichpos[-reg_odd]
  # str(motif$whichpos)
  # 
  # seq_names <- seq_names[-odd_lengths]
  # 
  mean(nchar(motif$seqs))
  
  # create matrix from sequence character vector
  ntmatrix <- matrix(unlist(strsplit(motif$seqs,"")), ncol = mean(nchar(motif$seqs)), byrow = TRUE)
  ntmatrix <- toupper(ntmatrix)
  rownames(ntmatrix) <- seq_names
  
  #highlight motif location using lowercase
  for(x in seq_along(motif$whichregs)){
    ntmatrix[motif$whichregs[x], motif$whichpos[x]:(motif$whichpos[x]+motif$scorematdim-1)] <- tolower(ntmatrix[motif$whichregs[x], motif$whichpos[x]:(motif$whichpos[x]+motif$scorematdim-1)])
  }
  
  return(ntmatrix)
}

# create colour pallete for creating heatmap of sequence_matrix
#dna_colours <- c('#109648', '#255C99', '#F7B32B', '#D62839')
#dna_colours2 <- c('#6fce97', '#92c3fc', '#fcda97', '#ffa5ae')
#dna_colours <- c(rbind(dna_colours, dna_colours2))
#names(dna_colours) <- c("a","A","c","C","g","G","t","T")
#brewer.pal(8, "Paired")[8:1], or brewer.pal(8, "Paired")[c(4,3,2,1,8,7,6,5)]




#' Extract matrix of motif occourences from sequence_matrix
#'
#' @param sequence_matrix output of sequence_matrix function
#' @param motif output of MotifFinder
#'
#' @details 
#' From a sequence matrix, extract occourences of the motif into a new sequence matrix
#' 
#' @return Charachter matrix (one motif occourence per row, one column per nucleotide) with rownames
#' 
#' @export

motif_consensus <- function(sequence_matrix, motif){
  # subset for sequences containing the motif
  newmatrix <- matrix(nrow = length(motif$whichregs), ncol = motif$scorematdim)
  
  for(x in seq_along(motif$whichregs)){
    newmatrix[x,] <- sequence_matrix[motif$whichregs[x], motif$whichpos[x]:(motif$whichpos[x]+motif$scorematdim-1)]
    
    # reverse complement if motif on opposite strand
    if(motif$whichstrand[x]==0){
      newmatrix[x,] <- chartr("atgc","tacg",newmatrix[x, ncol(newmatrix):1])
    }
    
  }
  
  rownames(newmatrix) <- rownames(sequence_matrix)[motif$whichregs]
  
  return(newmatrix)
}




#' For an SDA component, find motifs /de novo/
#'
#' @param component numeric; component number in results$loadings[[1]]
#' @param positive logical; should we look for motifs in the genes with positive loadings or negative
#' @param random logical; if TRUE (default: FALSE) 250 random genes are used
#'
#' @details 
#' Requires tss_seqs_s and results (from SDA load_results) to be loaded
#' For an SDA
#' 
#' @return Will save pdf of motifs found, and return list of motifs found
#'
#' @export
#' @import MotifFinder ggseqlogo

get_component_motifs <- function(component, positive, random=FALSE){
  
  gene_names17 <- head(names(sort(results$loadings[[1]][component,], decreasing = positive)),550)
  
  regs_spermio <- tss_seqs_s[gene_names17]
  regs_spermio <- regs_spermio[!is.na(regs_spermio)][1:250]
  stopifnot(sum(is.na(regs_spermio))==0)
  
  reigons <- rbind(
    c(0,250),
    c(0,200),
    c(0,150),
    c(-250,0),
    c(-200,0),
    c(-150,0),
    c(-100,100),
    c(-150,150),
    c(-200,200),
    c(-400,400)
  )
  
  sp1_motifs <- c("CCGCCC","CCCGCC","CGCCCC","CCCCGC")
  
  auto_motifs <- vector("list", nrow(reigons)*3*3)
  count <- 1
  
  if(positive){
    side<-"P"
  }else{
    side<-"N"
  }
  
  if(random){ # to get "Null" Hypothesis Motifs
    set.seed(42+component)
    gene_names_random <- sample(names(results$loadings[[1]][42,]),1500)
    
    regs_random <- tss_seqs_s[gene_names_random]
    regs_random <- regs_random[!is.na(regs_random)][1:250]
    stopifnot(sum(is.na(regs_random))==0)
    regs_spermio <- regs_random
    
    side <- "random"
  }
  
  
  pdf(file = paste0(component,side,".pdf"))
  
  for (i in 1:nrow(reigons)){
    for (j in 1:3){ # different motif ranks
      for (k in 1:3){ # repeat for multiple seeds
        
        print(paste("reigon",i))
        print(paste("motif rank",j))
        print(paste("seed",k))
        
        tmp <- findamotif(gsub("a|t|c|g","N",substr(regs_spermio, 500+reigons[i,1], 500+reigons[i,2])),
                          len=6, nits=100,
                          motif_blacklist = sp1_motifs,
                          motif_rank = j)
        
        if(!is.null(tmp)){
          auto_motifs[[count]] <- tmp
          print(
            ggseqlogo(list(get_PWM(auto_motifs[[count]]),
                           get_PWM(auto_motifs[[count]],T)),
                      ncol=1) +
              ggtitle(paste0("id: ",count,
                             "reigon: c(",reigons[i,1],",",reigons[i,2],
                             "), motif_rank:",j,
                             ", seed:",auto_motifs[[count]]$seed)) + ylim(0,2)
          )
          
          auto_motifs[[count]]$command <- paste(reigons[i,1],reigons[i,2],j,k,sep = ",")
          
          print(plot(auto_motifs[[count]]$alphas))
          print(plot(auto_motifs[[count]]$prior))
          print(plot_motif_location(auto_motifs[[count]]))
        }
        
        count <- count+1
        
      }
    }
  }
  
  dev.off()
  
  return(auto_motifs)
  
}



#' Export PWMs in MEME format
#'
#' @param results_all list of outputs from MotifFinder
#' @param file string; name of file to save PWMs to
#'
#' @details Export list of MotifFinder results to MEME format file, e.g. for input into TOMTOM
#' 
#' @return nothing, results written to file
#'
#' @export

export_all_pwms <- function(results_all, file){
  
  meme_header <- paste0("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n")
  write.table(meme_header, file,row.names=F, col.names = F, quote = F) # write header
  
  for(i in 1:length(results_all)){
    print(names(results_all[i]))
    
    if(!is.null(results_all[[i]])){ #  & mode(results[[j]][[i]])!="character"
      pwm <- get_PWM(results_all[[i]])
      
      motif_header <- paste0("\nMOTIF ",names(results_all[i]),
                             "\n\nletter-probability matrix: alength= 4 w= ",
                             ncol(pwm)," nsites= 20 E= 0")
      
      write.table(motif_header, file,row.names=F, col.names = F, quote = F, append = T)
      
      write.table(format(t(pwm), digits=5), file, row.names=F, col.names = F, quote = F, append = T)
    }
    
  }
  
}



#' For denovo motif and known motif, plot sequence logos
#'
#' @param best_denovo integer; motif id in object 'hocomoco'
#' @param titles logical; if TRUE (default) the logos titles are shown using $Target.Name and $X.Query.ID from best_denovo
#' @param yaxis logical; if TRUE (default: FALSE) yaxis is shown
#'
#' @details 
#' requires the object results_all
#' best_denovo is a data.table with column X.Query.ID with value e.g. C13N.I47 matching names of results_all list
#' also needs orientation etc.
#' 
#' This function is a wrapper for the plot_tomtom_match function from MotifFinder
#' 
#' @return list of ggplot2 objects
#'
#' @export
#' @import MotifFinder

plot_motif_matches <- function(best_denovo, titles=TRUE, yaxis=FALSE){
  
  matched_pwms <- list()
  
  for(i in 1:nrow(best_denovo)){
    print(i)
    print(best_denovo[i]$X.Query.ID)
    
    if(titles==T){
      title <- paste0(best_denovo[i]$Target.Name," - (",best_denovo[i]$X.Query.ID,")")
    }else{ # if(!is.null(titles))
      title <- titles[i]
    }
    
    matched_pwms[[i]] <- plot_tomtom_match(query_motif = results_all[[best_denovo[i]$X.Query.ID]],
                                           tomtom_match = best_denovo[i],
                                           yaxis=yaxis,
                                           titles=title)
  }
  
  return(matched_pwms)
  
}



# from ggseqlogo
information_content <- function(x){
  if(!is.null(results_all[[x]])){
    return(sum(computeBits(t(exp(results_all[[x]]$scoremat)))))
  }else{
    return(NA)
  }
}

# from ggseqlogo
computeBits <- function(pwm, N=4, Nseqs=NULL){
  Nseqs = attr(pwm, 'nongapped')
  H_i = - apply(pwm, 2, function(col) sum(col * log2(col), na.rm=T))
  e_n = 0
  if(!is.null(Nseqs)) e_n = (1/logb(2)) * (N-1)/(2*Nseqs) 
  
  R_i = log2(N) - (H_i  + e_n)
  R_i = pmax(R_i, 0)
  # Set any negatives to 0
  return(R_i)
}




#' For known motifs, run MotifFinder to find motif presence probability
#'
#' @param motif integer; motif id in object 'hocomoco'
#' @param database string; either "hocomoco" (default), or "jaspar"
#' @param custom logical; if TRUE (default: FALSE) the motif parameter should be a PWM
#' @param seqs charachter vector; DNA sequences (ALL UPPERCASE) to search motif presence in
#'
#' @details 
#' requires the object hocomoco or MotifDb to be loaded (depending on the database set)
#' 
#' @return
#'
#' @export
#' @import MotifFinder

find_motifs_parallel <- function(motif, database="hocomoco", custom=F, seqs=tss_seqs_s){
  
  if(database=="jaspar"){
    single_motif <- as.list(subset(MotifDb,
                                   dataSource=='jaspar2018' &
                                     organism %in% c("Hsapiens","Mmusculus","NA")
    )[motif])
    
    tmp <- single_motif[[1]]
    tmp <- t((tmp+0.01)/colSums(tmp+0.01))
    
  }else if(custom){
    single_motif <- motif
    tmp <- t(single_motif)
    str(tmp)
  }else { # hocomoco
    single_motif <- hocomoco[motif]
    tmp <- pcm2pwm(t(hocomoco[[motif]]))
  }
  
  start_time <- proc.time()
  
  found_motif_tmp <- getmotifs(log(tmp),
                               nrow(tmp),
                               seqs, # should be upper already
                               maxwidth=max(nchar(seqs)),
                               alpha=0.6, # prev 0.2 for jaspar run
                               maxits=20,
                               updatemot=0,
                               ourprior=rep(0.1,10),
                               seed=671747)
  
  found_motif_tmp$time <- proc.time() - start_time
  
  found_motif_tmp$motifID <- names(single_motif)
  
  found_motif_tmp$seqs_digest <- digest::digest(found_motif_tmp$seqs)
  
  found_motif_tmp$seqs <- NULL
  found_motif_tmp$trimmedseqs <- NULL
  
  return(found_motif_tmp)
  
}




#' For a given component, display binned associations
#'
#' @param component numeric or string; rowname of `results$loadings[[1]]` corresponding to a component name or number
#'
#' @details 
#' Requires regprobs_matrix to be loaded.
#' Calculate mean regulation probability of motif for binned gene loadings from SDA. Plot is faceted by motif.
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ggplot2

print_binned_meanprob <- function(component){
  tmp <- names(sort(results$loadings[[1]][component,]))
  tmp <- tmp[tmp %in% rownames(regprobs_matrix)]
  
  # for each motif, calculate binned mean regprob
  bin_factor <- cut(seq_along(regprobs_matrix[,1]), 50)
  binned_meanprob <- matrix(unlist(lapply(1:ncol(regprobs_matrix),
                                          function(i) tapply(regprobs_matrix[,i][tmp], bin_factor, function(x) mean(x,na.rm=T)))),
                            ncol=ncol(regprobs_matrix))
  colnames(binned_meanprob) <- colnames(regprobs_matrix)
  binned_meanprob <- data.table(binned_meanprob)
  binned_meanprob$bin <- seq_along(levels(bin_factor))
  binned_meanprob <- melt(binned_meanprob, id.vars = "bin", variable.name = "TF", value.name = "mean_regprob")
  
  return(
    ggplot(binned_meanprob, aes(bin, mean_regprob)) + geom_point() + geom_line() + facet_wrap(~TF, scales = "free_y") + xlab("Bin of SDA Gene Loading")
  )
}




#' For a matrix of genes by motifs, calculate correlations
#'
#' @param regprobs matrix; a matrix of genes by motifs
#'
#' @details 
#' Calculates two-sided Pearson's product-moment correlation test. Positive and negative loadings are split and calculated seperately.
#' 
#' @return
#' Matrix of components by motifs, with entries for the test statistic t for each combination
#'
#' @export

correlate_regprobs <- function(regprobs=regprobs_matrix){
  
  pred2 <- t(results$loadings[[1]][,rownames(regprobs_matrix)])
  
  vvv2 = matrix(nrow=100, ncol=ncol(regprobs))
  for(j in 1:ncol(pred2)){
    # positive side
    vvv2[j,] = unlist(lapply(1:ncol(regprobs), function(i) cor.test(regprobs[which(pred2[,j]>0),i],pred2[which(pred2[,j]>0),j])$statistic))
    
    # negative side
    vvv2[j+50,] = unlist(lapply(1:ncol(regprobs), function(i) cor.test(regprobs[which(pred2[,j]<0),i],abs(pred2[which(pred2[,j]<0),j]))$statistic))
  }
  
  dimnames(vvv2) <- dimnames(t(cor(regprobs,pred)))
  colnames(vvv2) <- sapply(strsplit(colnames(vvv2), "_"),function(x) x[1])
  
  return(vvv2)
}
