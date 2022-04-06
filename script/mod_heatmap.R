mod_heatmap_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3("The main options:"),
      fileInput(ns("filename"),"Choose File to Upload:", accept = c(".txt")),
      radioButtons(ns('filetype'), 'File type', c(xlsx='excel', txt='table', csv='CSV'), 'excel'),
      selectInput(ns("color"), label = "Select color for data(*)",
                  c("default", "green yellow red", "blue yellow red")),
      
      selectInput(ns("scale"), label = "Select scale for data",
                  c("row", "column", "none")),
      selectInput(ns("clusterrows"), label = "Select cluster rows for cluster(*)",
                  c("TRUE", "FALSE")),
      selectInput(ns("clustercols"), label = "Select cluster columns for cluster(*)",
                  c("TRUE", "FALSE")),
      selectInput(ns("algorithm"), label = "Select algorithm to cluster",
                  c("complete", "single", "average", "centroid", "median")),
      selectInput(ns("rownames"), label = "Select show row name(*)",
                  c("TRUE", "FALSE")),
      selectInput(ns("colnames"), label = "Select show column name(*)",
                  c("TRUE", "FALSE")),
      selectInput(ns("colnames_angle"), label = "Select angle of column name",
                  c(270, 0, 45, 90, 315)),
      selectInput(ns("displaynumber"), label = "Select show number in cell(*)",
                  c("FALSE", "TRUE")),
      numericInput(ns("cellheight"),label = "Cell heigh value",value = 5),
      numericInput(ns("cellwidth"),label = "Cell width value",value = 30),
      br(),
      actionButton(ns("pheatmap_submit"), strong("Submit"), styleclass = "success"),
      helpText("After re-selecting the options with '*', you need to click submit!!!"),
      br(),
      br(),
      h3(strong("Download options:")),
      numericInput(ns("heatmapheigh"),label = "Graph heigh value",value = 20),
      numericInput(ns("heatmapwidth"),label = "Graph width value",value = 10),
      downloadButton(ns("downloadheatmap"),label = "Download heatmap")
    ),
    
    mainPanel(
      h3("heatmap plot:"),
      withSpinner(plotOutput(ns("heatmapplot"), width='100%', height='2000px'))
    )
  )
}


mod_heatmap_server <- function(input, output, session){
  ns <- session$ns
  
  filedata <- eventReactive(input$pheatmap_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }
    if (input$filetype == "table"){
      read.table(infile$datapath,sep = "\t", header = T, row.names = 1)
    }else{
      if(input$filetype == "excel"){
        output<-read.xlsx(infile$datapath,1, header=T)
        row.names(output)<-output[,1]
        output[,-1]
      }else{
        read.csv(infile$datapath,sep=',', header=T, row.names = 1)
      }
    }
    
  })
  
  show_clusterrows <- eventReactive(input$pheatmap_submit,{
    clusterrowname <- input$clusterrows
    if (clusterrowname == "FALSE")
      return( FALSE )
    else
      return( TRUE )
  })
  
  show_clustercols <- eventReactive(input$pheatmap_submit,{
    clustercolname <- input$clustercols
    if (clustercolname == "FALSE")
      return( FALSE )
    else
      return( TRUE )
  })
  
  show_rowname <- eventReactive(input$pheatmap_submit,{
    rowname <- input$rownames
    if (rowname == "FALSE")
      return( FALSE )
    else
      return( TRUE )
  })
  
  show_colname <- eventReactive(input$pheatmap_submit,{
    colname <- input$colnames
    if (colname == "FALSE")
      return( FALSE )
    else
      return( TRUE )
  })
  
  show_number <- eventReactive(input$pheatmap_submit,{
    number <- input$displaynumber
    if (number == "FALSE")
      return( FALSE )
    else
      return( TRUE )
  })
  

  show_color <- eventReactive(input$pheatmap_submit,{
    colorname <- input$color
    if (colorname == "green yellow red")
      return( colorRampPalette(c("green4", "yellow", "red"))(500) )
    if (colorname == "blue yellow red")
      return( colorRampPalette(c("navy", "yellow", "red"))(500) )
  })
  
  heatmap_plot <- eventReactive(input$pheatmap_submit,{
    df <- filedata()
    
    if(input$color == "default"){
      pheatmap(df, 
               scale= input$scale, angle_col= input$colnames_angle,
               cluster_rows = show_clusterrows(), cluster_cols = show_clustercols(),
               cellwidth=input$cellwidth, cellheight= input$cellheight, 
               show_rownames = show_rowname(), show_colnames = show_colname(), display_numbers=show_number(),
               clustering_method_rows = input$algorithm
      )
    }else{
      pheatmap(df, color = show_color(),
               scale= input$scale, angle_col= input$colnames_angle,
               cluster_rows = show_clusterrows(), cluster_cols = show_clustercols(),
               cellwidth=input$cellwidth, cellheight= input$cellheight, 
               show_rownames = show_rowname(), show_colnames = show_colname(), display_numbers=show_number(),
               clustering_method_rows = input$algorithm
      )
    }
  })
  
  output$heatmapplot <- renderPlot({
    heatmap_plot()
  })

  output$downloadheatmap <- downloadHandler(
    filename = function() { 
      paste0("heatmap_plot", '.pdf')
    },
    contentType = "image/pdf",
    content = function(file) {
      pdf(file, width = input$heatmapwidth, height = input$heatmapheigh)
      print(heatmap_plot())
      dev.off()
    }
  )
  
}
