mod_mirna_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3("The main options of miRNA structure:"),
      selectInput(ns("mirnaid"), "Please select a MIRNA:", choices = mirna_list),
      numericInput(ns("pointsize"),label = "Point size value",value = 1),
      numericInput(ns("textsize"),label = "Nucleotide size value",value = 0.6),
      actionButton(ns("mirna_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3("The main options of download:"),
      numericInput(ns("mirnaheigh"),label = "Graph heigh value",value = 10),
      numericInput(ns("mirnawidth"),label = "Graph width value",value = 10),
      br(),
      br(),
      downloadButton(ns("downloadmirna"),label = "Download miRNA Secondary Structure"),
      br(),
      br(),
      downloadButton(ns("downloadmirnaseq"),label = "Download miRNA Sequence")
      ),
    
    mainPanel(
      h3("miRNA Secondary Structure:"),
      withSpinner(plotOutput(ns("mirnaplot"), width='60%', height='800px'))
      
    )
  )
}


mod_mirna_server <- function(input, output, session){
  ns <- session$ns
  
  mirna_plot <- eventReactive(input$mirna_submit,{
    ct <- makeCt(mirna_info[input$mirnaid,9],mirna_info[input$mirnaid,8])
    dat <- ct2coord(ct)
    RNAPlot2(dat,pointSize=input$pointsize)
    RNAPlot2(dat,nt=TRUE,tsize=input$textsize, add=TRUE,dp=1,hl=c(as.character(mirna_info[input$mirnaid,6]),as.character(mirna_info[input$mirnaid,7])),seqcol=c(2,4),labTF=TRUE)
  })
  
  output$mirnaplot <- renderPlot({
    mirna_plot()
  })
  
  output$downloadmirna <- downloadHandler(
    filename = function() { 
      a<-paste0("mirnaplot", '.pdf')
      paste0(input$mirnaid, a)
    },
    contentType = "image/pdf",
    content = function(file) {
      pdf(file, width = input$mirnawidth, height = input$mirnaheigh)
      print(mirna_plot())
      dev.off()
    }
  )
  
  output$downloadmirnaseq <- downloadHandler(
    filename = function(){
      a<-paste0("mirnaplot", '.csv')
      paste0(input$mirnaid, a)
    },
    content=function(file){
      write.csv(mirna_info[input$mirnaid,],file)
    }
  )
}
