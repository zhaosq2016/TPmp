mod_hcsiRNA_ui <- function(id){
  ns <- NS(id)
  tagList(
    mainPanel(
      width = 12,
      h2("hcsiRNA information:"),
      withSpinner(DTOutput(ns("hcsirnainfo")))
    )
  )
}

mod_hcsiRNA_server <- function(input, output, session) {
  ns <- session$ns
  output$hcsirnainfo <- renderDT({
    hcsirna
  })
}
