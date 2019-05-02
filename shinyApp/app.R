
##
# for running this app on a memory constrained machine (1GB RAM) set the following option to TRUE
# In low memory mode, the data will be stored on disk and read when requested by the end user
# To use low memory mode, the script setup_low_memory_mode.R must be run

low_memory_mode = TRUE

##

library(data.table)
library(ggplot2)
library(viridis)
library(SDAtools)
library(testisAtlas)
library(shinycssloaders)


if(low_memory_mode){
  library(bigmemory)
  options(bigmemory.allow.dimnames=TRUE)

  load("cached_objects_small.rds", envir = .GlobalEnv)
  assign('data', attach.big.matrix("data.big.desc"), envir=.GlobalEnv)
  #data <- attach.big.matrix("data.big.desc")

}else{
  load("cached_objects.rds", envir = .GlobalEnv)
  load_component_orderings()
}

results$loadings_ranks <- t(sapply(rownames(results$loadings[[1]]), function(i) rank(abs(results$loadings[[1]][i,]))))

gene_symbols <- readRDS("Ensembl_92_gene_symbols.rds")

# clip outlying cell
datat[msci_ratio>0.15,msci_ratio := 0.12]

gene_or_component <- function(tmp, allways_component=F){
  gene <<- NULL
  if(!is.na(as.numeric(tmp))){
    tmp <- as.numeric(tmp)
  }else{
    tmp <- strsplit(tmp, " ")[[1]]
    
    if(length(tmp)==0){
      tmp <- c("Prdm9")
    }
    
    if(!tmp %in% colnames(datat)){
      tmp <- paste0(toupper(substring(tmp, 1, 1)), substring(tmp, 2))
    }
    
    if(length(which(!tmp %in% c(gene_symbols$external_gene_name,colnames(datat))))!=0){
      stop(paste("Gene Symbol not recognised:", paste(tmp[which(!tmp %in% colnames(data))], collapse = "",sep = "")))
    }else if(length(which(!tmp %in% c(colnames(data),colnames(datat))))!=0){
      stop(paste("Gene not detected:", paste(tmp[which(!tmp %in% colnames(data))], collapse = "",sep = "")))
    }
    
    if(allways_component){
      gene <<- tmp
      tmp <- which.max(results$loadings_ranks[,gene])
    }
    
  }
  
  return(tmp)
}

options(shiny.sanitize.errors = FALSE)

server <- function(input, output, session) {
  
  library(testisAtlas)
  
  ranges <- reactiveValues(x = NULL, y = NULL)

  output$distPlot <- renderPlot({
    
    tmp <- gene_or_component(trimws(input$genes))
    
    if(!input$umap){
    p <- print_tsne(tmp,
               predict = input$show_predict,
               curve = input$show_arrow,
               stages = input$show_stages,
               point_size = input$decimal) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = T)
    }else{
      p <- print_tsne(tmp,
                      dim1="Umap1", dim2="Umap2",
                      predict = input$show_predict,
                      curve = FALSE,
                      stages = FALSE,
                      point_size = input$decimal) +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = T)
    }
    
    p <- p + annotate("text", x=-Inf, y=-Inf, hjust=0, vjust=-1, label="Jung & Wells et. al. 2018", colour='grey')
    
    return(p)
    
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$distPlot_dblclick, {
    brush <- input$distPlot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  # select (pseudo)random gene
  observeEvent(input$rand_gene, {
    # Exclude boring components
    #dput(component_order_dt[QC_fail==FALSE]$component_number)
    good_components <- c(2, 3, 5, 7, 10, 11, 13, 15, 16, 17, 18, 19, 20, 21, 23, 24, 
                         26, 27, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 44, 45, 
                         47, 48, 49, 50)
    # weight component choice by mean cell score of component
    c <- sample(good_components, size = 1,
                prob = results$component_statistics[good_components]$mean_score /
                  sum(results$component_statistics[good_components]$mean_score))
    new_gene <- names(sort(abs(results$loadings[[1]][c,]), T)[sample(1:200,1)])

    updateTextInput(session, inputId="genes", value = new_gene)
  })
  
  output$loadings <- renderPlot({
    
    tmp <- gene_or_component(input$genes, allways_component = T)
    
    if(!is.null(gene)){
      title <- paste("Component",tmp,"(",component_order_dt[component_number==tmp]$name,")","- highest ranked loading for",gene)
    }else{
      title <- paste("Component",tmp,"(",component_order_dt[component_number==tmp]$name,")")
    }
    
	genome_loadings(results$loadings[[1]][tmp,], label_both = TRUE, max.items = input$n_genes, label.size = 4, gene_locations=rna_locations, highlight_genes = gene, label_genes = gene) +
	ylab(paste("Gene Loading (Component",tmp,")")) + theme_minimal() + theme(legend.position = "none") + ggtitle(title)
    
  })
  
  output$whichComp <- renderPlot({
    
    tmp <- trimws(input$genes)
    
    if(!is.na(as.numeric(tmp))){
      stop("Please enter a gene name")
    }else{
      tmp <- strsplit(tmp, " ")[[1]]
      tmp <- paste0(toupper(substring(tmp, 1, 1)), substring(tmp, 2))
      
      if(length(tmp)==0){
        tmp <- c("Prdm9")
      }else if(length(which(!tmp %in% gene_symbols$external_gene_name))!=0){
        stop(paste("Gene Symbol not recognised:", paste(tmp[which(!tmp %in% colnames(data))], collapse = "",sep = "")))
      }else if(length(which(!tmp %in% colnames(data)))!=0){
        stop(paste("Gene not detected:", paste(tmp[which(!tmp %in% colnames(data))], collapse = "",sep = "")))
      }
    }
    
    highest_components(results, tmp)
    
  })

  
}

ui <- fluidPage(
  titlePanel("Mouse Testis Single Cell RNAseq Atlas"),
  sidebarLayout(
    sidebarPanel(
      textInput("genes", "Input gene name / component number:", placeholder = "Prdm9"),
      tags$div(class="label", style="text-align:centered; color:black; font-size:100%;", checked=NA, tags$strong("OR")),
      actionButton("rand_gene", "Show me a random gene!"),
      tags$hr(),
      tags$div(class="header", checked=NA, tags$h4("tSNE Settings:")),
      helpText("To Zoom, select a reigon and double click, to reset double click."),
      checkboxInput("show_predict", "Show imputed expression?", value=TRUE),
      checkboxInput("show_arrow", "Display Pseudotime Arrow?"),
      checkboxInput("show_stages", "Annotate Stages?", value=TRUE),
      checkboxInput("umap", "Umap Projection?"),
      sliderInput("decimal", "Point size:",
                  min = 0.1, max = 3, value = 1.5, step = 0.1),
      tags$hr(),
      tags$div(class="header", checked=NA, tags$h4("Genome Loadings Settings:")),
      sliderInput("n_genes", "Number of genes to label:",
                  min = 0, max = 100, value = 20, step = 1),
      tags$head(tags$style("#distPlot{height:85vh !important;}")),
      tags$head(tags$style("#loadings{height:85vh !important;}")),
      tags$head(tags$style("#whichComp{height:85vh !important;}")),
      tags$hr(),
      tags$div(class="header", checked=NA, tags$h4("About:")),
      tags$div(p(HTML(paste0('For more information see our ',a(target="_blank", href = 'https://doi.org/10.1101/393769', 'paper')))))
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("tSNE",
                           withSpinner(plotOutput("distPlot",
                                      height = "auto",
                                      dblclick = "distPlot_dblclick",
                                      brush = brushOpts(
                                        id = "distPlot_brush",
                                        resetOnNew = TRUE)
                           ))
                  ),
                  tabPanel("Gene Loadings",
                           plotOutput("loadings")),
                  tabPanel("Which Component",
                           plotOutput("whichComp"))
      )
    )
  )
)

shinyApp(ui = ui, server = server)