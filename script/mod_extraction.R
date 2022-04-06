mod_extraction_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 3,
      h3("The main options:"),
      selectInput(ns("srnatype"), "sRNA type:", choices = c("miRNA", "phasiRNA21", "phasiRNA24", "hcsiRNA")),
      selectInput(ns("chrom"), "Chromosome:", choices = allresult),
      numericInput(ns("chromstart"),label = "Start:",value = 10000000),
      numericInput(ns("chromend"),label = "End:",value = 20000000),
      actionButton(ns("id_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      downloadButton(ns("downloadsrnadata"),label = "Download sRNA information:")
      ),
    mainPanel(
      h2("sRNA information:"),
      withSpinner(DTOutput(ns("srnadata")))
    )
  )
}

mod_extraction_server <- function(input,output,session){
  ns <- session$ns
  srnadata1 <- eventReactive(input$id_submit,{
    if(input$srnatype == "miRNA"){
      mirna_info%>% dplyr::filter(Chr_ID==input$chrom, 
                                                Precursor_Start >= input$chromstart,
                                                Precursor_End <= input$chromend
                                                )
    }else if(input$srnatype == "phasiRNA21"){
      phasiRNA21 %>% dplyr::filter(Chr_ID==input$chrom, 
                                                Start_Pos >= input$chromstart,
                                                End_Pos <= input$chromend
      )
    }else if(input$srnatype == "phasiRNA24"){
      phasiRNA24 %>% dplyr::filter(Chr_ID==input$chrom, 
                                                Start_Pos >= input$chromstart,
                                                End_Pos <= input$chromend
      )
    }else if(input$srnatype == "hcsiRNA"){
      hcsirna %>% dplyr::filter(Chr_ID==input$chrom, 
                                                Start_Pos >= input$chromstart,
                                                End_Pos <= input$chromend
      )
    }
  })
  output$srnadata <- renderDT({
    srnadata1()
  })
  
  output$downloadsrnadata <- downloadHandler(
    filename = function(){
      paste0("sRNA_data", '.csv')
    },
    content=function(file){
      write.csv(srnadata1(),file)
    }
  )
  
}