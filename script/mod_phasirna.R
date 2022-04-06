mod_phasirna_ui <- function(id){
  warnings('off')
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 3,
      h3("The main options:"),
      selectInput(ns("phasirnalength"), "Select PHAS type:", choices = c("21nt", "24nt")),
      selectInput(ns("phasirnasample"), "Select a sample:", choices = phasiRNAsample),
      selectInput(ns("phasirnaid"), "Select a PHAS:", choices = phasirna21_list),
      br(),
      actionButton(ns("phasirna_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h4("Download the align and phasing scores of PHAS:"),
      downloadButton(ns("downloadphasirna"), "Download")
      ),
    mainPanel(
      h3("PHAS information:"),
      withSpinner(tableOutput(ns("table"))),
      h3("Align and Phasing scores of PHAS:"),
      withSpinner(plotOutput(ns("plotphas"), width='60%', height='800px'))
    )
  )
}

mod_phasirna_server <- function(input,output,session){
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
  
  phasirna_table<- reactive({
    if(input$phasirnalength == "21nt"){
      phasiRNA21
    }else{
      phasiRNA24
    }
  })


  plotphasirna <- eventReactive(input$phasirna_submit,{
    url <- paste0("https://gitee.com/zhaoshiqi89/phasirna/raw/master/pdf/",input$phasirnasample,"/PHASList_",input$phasirnalength,"_ToGraph/",input$phasirnaid, ".pdf")
    download.file(url,destfile = "data/pdf/phasiRNA.pdf", quiet = TRUE, mode="wb")
    m <- "data/pdf/phasiRNA.pdf"
    
    if(file_test("-f",m)== TRUE){
      bitmap <- pdftools::pdf_render_page(m)
      png::writePNG(bitmap, "www/img/test.png", dpi=500)
    }else{
      bitmap <- pdftools::pdf_render_page("data/pdf/no.pdf")
      png::writePNG(bitmap, "www/img/test.png", dpi=500)
    }
    imagefile <- "www/img/test.png"
    list(src=imagefile, width=1000, height=700)
  })
  
  phasirnatable <- eventReactive(input$phasirna_submit,{
    phasirna_table()[input$phasirnaid,1:6]
  })
  
  output$table <- renderTable({
    phasirnatable()
  })
  
  
  output$plotphas <- renderImage({
    plotphasirna()
  })
  
  output$downloadphasirna <- downloadHandler(
    filename = function() { 
      paste0(input$phasirnaid, '.pdf')
    },
    content = function(file){
      m <- "data/pdf/phasiRNA.pdf"
      if(file_test("-f",m)== TRUE){
        file.copy("data/pdf/phasiRNA.pdf", file)
      }else{
        file.copy("data/pdf/no.pdf", file)
      }
    }
  )
}