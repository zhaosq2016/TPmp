mod_target_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3("The main options:"),
      fileInput(ns("filename1"),"Choose Gene Exon File to Upload(*):", accept = c(".txt")),
      radioButtons(ns('filetype1'), 'File type', c(xlsx='excel', txt='table', csv='CSV'), 'excel'),
      
      fileInput(ns("filename2"),"Choose miRNA_target File to Upload:", accept = c(".txt")),
      radioButtons(ns('filetype2'), 'File type', c(xlsx='excel', txt='table', csv='CSV'), 'excel'),
      
      numericInput(ns("splicesite"),label = "The Splice Site Number",value = 10),
      
      actionButton(ns("target_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      numericInput(ns("targetheigh"),label = "Graph heigh value",value = 9.5),
      numericInput(ns("targetwidth"),label = "Graph width value",value = 9.5),
      downloadButton(ns("downloadtarget"),label = "Download picture")
    ),
    
    mainPanel(
      h3("Target plot:"),
      withSpinner(plotOutput(ns("targetplot"), width='80%', height='800px'))
    )
  )
}


mod_target_server <- function(input, output, session){
  ns <- session$ns
  
  filedata1 <- eventReactive(input$target_submit,{
    infile1 <- input$filename1
    if (is.null(infile1)){
      return(NULL)
    }
    if (input$filetype1 == "table"){
      read.table(infile1$datapath,sep = "\t", header = T)
    }else{
      if(input$filetype1 == "excel"){
        output1<-read.xlsx(infile1$datapath,1, header=T)
        row.names(output1)<-output1[,1]
        output1[,-1]
      }else{
        read.csv(infile1$datapath,sep=',', header=T)
      }
    }
    
  })
  
  filedata2 <- eventReactive(input$target_submit,{
    infile2 <- input$filename2
    if (is.null(infile2)){
      return(NULL)
    }
    if (input$filetype2 == "table"){
      read.table(infile2$datapath,sep = "\t", header = T, quote = "")
    }else{
      if(input$filetype2 == "excel"){
        output2<-read.xlsx(infile2$datapath,1, header=T)
        row.names(output2)<-output2[,1]
        output2[,-1]
      }else{
        read.csv(infile2$datapath,sep=',', header=T)
      }
    }
    
  })
  
  mirtarget_plot <- function(){
    target_plot(filedata1(), filedata2(), input$splicesite)
  }
  
  
  output$targetplot <- renderPlot({
    mirtarget_plot()
  })

  output$downloadtarget <- downloadHandler(
    filename = function() { 
      paste0("target_plot", '.pdf')
    },
    contentType = "image/pdf",
    content = function(file) {
      pdf(file, width = input$targetwidth, height = input$targetheigh)
      print(mirtarget_plot())
      dev.off()
    }
  )
  
}
