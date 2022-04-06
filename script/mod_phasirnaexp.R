mod_phasirnaexp_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3("The main options:"),
      selectInput(ns("phasirnalength"), "Please select PHAS type:", choices = c("21nt", "24nt")),
      selectInput(ns("phasirnaid"), "Please select a MIRNA:", choices = phasirna21_list),
      
      selectInput(ns("plottype"), "Please select plot type:", choices = c("line chart", "bar chart")),
      br(),
      h4(strong("Select samples:")),
      chooserInput(ns("sample_accession"), "Available frobs", "Selected frobs", phasiRNAsample, c(), size = 10, multiple = TRUE),
      br(),
      actionButton(ns("id_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      numericInput(ns("phasirnaexpheigh"),label = "Graph heigh value",value = 10),
      numericInput(ns("phasirnaexpwidth"),label = "Graph width value",value = 10),
      downloadButton(ns("downloadphasirnaexp"),label = "Download phasiRNA expression plot"),
      br(),
      br(),
      downloadButton(ns("downloadphasirnaexpdata"),label = "Download phasiRNA expression data")
    ),
    
    mainPanel(
      h3("phasiRNA expression plot:"),
      withSpinner(plotOutput(ns("phasirnaexpplot"), width='60%', height='800px'))
    )
  )
}


mod_phasirnaexp_server <- function(input, output, session){
  ns <- session$ns
  
  phasirna_list<- reactive({
    if(input$phasirnalength == "21nt"){
      phasirna21_list
    }else{
      phasirna24_list
    }
  })
  
  observe({
    x <- phasirna_list()
    updateSelectInput(session, "phasirnaid",
                      label = "Please select a PHAS: ",
                      choices = x,
                      selected = head(x, 1)
    )
  })
  
  phasiRNAexp_data <- eventReactive(input$id_submit,{
    
    if(input$phasirnalength == "21nt"){
      phasiRNAexp <- phasiRNA21exp
    }else{
      phasiRNAexp <- phasiRNA24exp
    }
    
    data_exp <- as.data.frame(t(phasiRNAexp[input$phasirnaid, input$sample_accession$selected]))
    data_exp[,2] <- row.names(data_exp)
    colnames(data_exp) <- c("expression", "sample")
    data_exp
  })
  
  phasirnaexp_plot <- function(){
    if(input$plottype == "line chart"){
      ggplot(phasiRNAexp_data(), aes(x=sample, y=expression, group=1)) + 
        geom_line() + 
        geom_point(size=4, shape=20) + 
        xlab("Sample") + 
        ylab("Expression of phasiRNA") + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              axis.text=element_text(size=12,face = "bold"),
              axis.title=element_text(size=15,face = "bold")
              )
    }else{
      ggplot(phasiRNAexp_data(), aes(x=sample, y=expression, group=factor(1))) + 
        geom_bar(stat="identity") + 
        xlab("Sample") + 
        ylab("Expression of phasiRNA") + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              axis.text=element_text(size=12,face = "bold"),
              axis.title=element_text(size=15,face = "bold")
              )
    }
    
  }
  
  output$phasirnaexpplot <- renderPlot({
    phasirnaexp_plot()
  })

  output$downloadphasirnaexp <- downloadHandler(
    filename = function() { 
      a<-paste0("_expplot", '.pdf')
      paste0(input$phasirnaid, a)
    },
    contentType = "image/pdf",
    content = function(file) {
      pdf(file, width = input$phasirnaexpwidth, height = input$phasirnaexpheigh)
      print(phasirnaexp_plot())
      dev.off()
    }
  )
  
  output$downloadphasirnaexpdata <- downloadHandler(
    filename = function(){
      a<-paste0("_data", '.csv')
      paste0(input$phasirnaid, a)
    },
    content=function(file){
      write.csv(phasiRNAexp_data(),file)
    }
  )
}
