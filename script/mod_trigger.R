mod_trigger_ui <- function(id){
  warnings('off')
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 3,
      h3("The main options:"),
      selectInput(ns("phasirnalength"), "Please select PHAS type:", choices = c("21nt", "24nt")),
      selectInput(ns("phasirnaid"), "Please select a PHAS:", choices = phasirna21_list),
      selectInput(ns("triggertype"), "Please select the type of trigger:", choices = c("Trigger", "Potential")),
      br(),
      actionButton(ns("trigger_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      ),
    mainPanel(
      h3("Trigger information:"),
      withSpinner(tableOutput(ns("table")))
    )
  )
}

mod_trigger_server <- function(input,output,session){
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

  trigger_table<- reactive({
    if(input$phasirnalength == "21nt"){
      trigger_21
    }else{
      trigger_24
    }
  })
  
  triggertable <- eventReactive(input$trigger_submit,{
    a <- trigger_table()
    b <- a[input$phasirnaid, input$triggertype]
    d <- a[input$phasirnaid,1]
    d <- paste0(d, "\t")
    b <- gsub("),", "\n", b)
    b <- gsub("\\(", d, b)
    b <- gsub(" ", "\t", b)
    b <- gsub("\\)", "", b)
    if(!is.na(b) & !is.null(b) & b != ""){
      e <- read.table(text=b)
      colnames(e) <- c("Trigger", "phasiRNA_position", "Trigger_site", "Chain", "Trigger_sequence")
      e
    }else{
      "No trigger was identified!"
    }
  })
  
  output$table <- renderTable({
    triggertable()
  })
}