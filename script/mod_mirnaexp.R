mod_mirnaexp_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3("The main options:"),
      selectInput(ns("mirnaidexp"), "Please select a MIRNA:", choices = mirna_list3),
      
      selectInput(ns("plottype"), "Please select plot type:", choices = c("line chart", "bar chart")),
      br(),
      h4(strong("Select samples:")),
      chooserInput(ns("sample_accession"), "Available frobs", "Selected frobs", colnames(allmiRNAexp), c(), size = 10, multiple = TRUE),
      br(),
      actionButton(ns("id_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      numericInput(ns("mirnaexpheigh"),label = "Graph heigh value",value = 10),
      numericInput(ns("mirnaexpwidth"),label = "Graph width value",value = 10),
      downloadButton(ns("downloadmirnaexp"),label = "Download miRNA expression plot"),
      br(),
      br(),
      downloadButton(ns("downloadmirnaexpdata"),label = "Download miRNA expression data")
    ),
    
    mainPanel(
      h3("miRNA expression plot:"),
      withSpinner(plotOutput(ns("mirnaexpplot"), width='60%', height='800px'))
    )
  )
}


mod_mirnaexp_server <- function(input, output, session){
  ns <- session$ns
  
  miRNAexp_data <- eventReactive(input$id_submit,{
    data_exp <- as.data.frame(t(allmiRNAexp[input$mirnaidexp, input$sample_accession$selected]))
    data_exp[,2] <- row.names(data_exp)
    colnames(data_exp) <- c("expression", "sample")
    data_exp$sample <- factor(data_exp$sample, levels=data_exp$sample)
    data_exp
  })
  
  mirnaexp_plot <- function(){
    if(input$plottype == "line chart"){
      ggplot(miRNAexp_data(), aes(x=sample, y=expression, group=1)) + 
        geom_line() + 
        geom_point(size=4, shape=20) + 
        xlab("Sample") + 
        ylab("Expression of miRNA (TPM)") + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              axis.text=element_text(size=12,face = "bold"),
              axis.title=element_text(size=15,face = "bold")
              )
    }else{
      ggplot(miRNAexp_data(), aes(x=sample, y=expression, group=factor(1))) + 
        geom_bar(stat="identity") + 
        xlab("Sample") + 
        ylab("Expression of miRNA (TPM)") + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              axis.text=element_text(size=12,face = "bold"),
              axis.title=element_text(size=15,face = "bold")
              )
    }
    
  }
  
  output$mirnaexpplot <- renderPlot({
    mirnaexp_plot()
  })

  output$downloadmirnaexp <- downloadHandler(
    filename = function() { 
      a<-paste0("_expplot", '.pdf')
      paste0(input$mirnaidexp, a)
    },
    contentType = "image/pdf",
    content = function(file) {
      pdf(file, width = input$mirnaexpwidth, height = input$mirnaexpheigh)
      print(mirnaexp_plot())
      dev.off()
    }
  )
  
  output$downloadmirnaexpdata <- downloadHandler(
    filename = function(){
      a<-paste0("_data", '.csv')
      paste0(input$mirnaidexp, a)
    },
    content=function(file){
      write.csv(miRNAexp_data(),file)
    }
  )
}
