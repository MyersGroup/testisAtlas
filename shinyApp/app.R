
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

gene_symbols <- readRDS("Ensembl_92_gene_symbols.rds")

options(shiny.sanitize.errors = FALSE)

server <- function(input, output) {
  
  library(testisAtlas)
  
  ranges <- reactiveValues(x = NULL, y = NULL)

  output$distPlot <- renderPlot({
    
    tmp <- trimws(input$genes)
    
    if(!is.na(as.numeric(tmp))){
      tmp <- as.numeric(tmp)
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
    
    if(input$diverging_colour){
      return(p + scale_fill_gradient2())
    }else{
      return(p) 
    }
    
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
  
  output$loadings <- renderPlot({
    
    tmp <- input$genes
    
    if(!is.na(as.numeric(tmp))){
      tmp <- as.numeric(tmp)
    }else{
		stop("Please enter a component number")
    }
	
	genome_loadings(results$loadings[[1]][tmp,], label_both = TRUE, max.items = input$n_genes, label.size = 4, gene_locations=rna_locations) +
	ylab(paste("Gene Loading (Component",tmp,")")) + theme_minimal() + theme(legend.position = "none") +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    
  })

  
}

ui <- fluidPage(
  titlePanel("Mouse Testis Single Cell RNAseq Explorer"),
  sidebarLayout(
    sidebarPanel(
      textInput("genes", "Input gene name or component number:", placeholder = "Prdm9"),
      tags$hr(),
      tags$div(class="header", checked=NA, tags$h4("tSNE Settings:")),
      helpText("To Zoom, select a reigon and double click, to reset double click."),
      checkboxInput("show_predict", "Show imputed expression?", value=TRUE),
      checkboxInput("show_arrow", "Display Pseudotime Arrow?"),
      checkboxInput("show_stages", "Annotate Stages?", value=TRUE),
      checkboxInput("diverging_colour", "Diverging Colour Scale?"),
      checkboxInput("umap", "Umap Projection?"),
      sliderInput("decimal", "Point size:",
                  min = 0.1, max = 2, value = 1, step = 0.05),
      tags$hr(),
      tags$div(class="header", checked=NA, tags$h4("Genome Loadings Settings:")),
      sliderInput("n_genes", "Number of genes to label:",
                  min = 0, max = 100, value = 20, step = 1),
      tags$head(tags$style("#distPlot{height:88vh !important;}")),
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
                           plotOutput("loadings"))
      )
    )
  )
)

shinyApp(ui = ui, server = server)