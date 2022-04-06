mod_resource_ui <- function(id){
  ns <- NS(id)
  tagList(
    mainPanel(
      width = 12,
      h2("High-throughput sequencing data used in this study:"),
      withSpinner(DTOutput(ns("rnaseqinfo")))
    )
  )
}

mod_resource_server <- function(input, output, session) {
  ns <- session$ns
  output$rnaseqinfo <- renderDT({
    allsample
  })
}
