mod_mirtar_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3("The main options of miRNA targets:"),
      selectInput(ns("mirnaidseq"), "Please select a miRNA:", choices = mirna_list2),
      selectInput(ns("targetid"), "Please select a target id:", choices = c("CSS1", "CSS2")),
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


mod_mirtar_server <- function(input, output, session){
  ns <- session$ns
  
  miRNA_target <- reactive({
    m <- dplyr::filter(miRNAtarget, miRNAtarget[,2] == input$mirnaidseq)
  })
  
  observe({
    x <- miRNA_target()[,4]
    x <- unique(x)
    updateSelectInput(session, "targetid",
                      label = "Please select a target:",
                      choices = x,
                      selected = head(x, 1)
    )
  })
  
  targetgff <- reactive({
    targetid2 <- gsub("_.*", "", input$targetid)
    dplyr::filter(GFF_info, GFF_info[,9] == targetid2)
  })
  
  miRNA_target2 <- reactive({
    a <- dplyr::filter(miRNAtarget, miRNAtarget[,2] == input$mirnaidseq)
    dplyr::filter(a,a[,4] == input$targetid)
  })
  
  mirtarget_plot <- eventReactive(input$target_submit,{
    splicesite <- miRNA_target2()[1,5]
    target_plot(targetgff(), miRNA_target2(), splicesite)
  })

  output$targetplot <- renderPlot({
    showtext_auto()
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
