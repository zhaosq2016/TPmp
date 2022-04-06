source("global.R")

ui <- fluidPage(
  tags$script(src="css/addhash.js"),
  div(img(src = "img/TPmp_web.png")),
  includeCSS("www/css/custom.css"),
  includeCSS("www/css/footer.css"),
  disconnectMessage(
    text = "Your session timed out, reload the application!",
    refresh = "Reload now",
    background = "#f89f43",
    colour = "white",
    overlayColour = "grey",
    overlayOpacity = 0.75,
    top = 250,
    refreshColour = "brown"
  ),
  navbarPage(
    title = "",
    windowTitle = "TPmp",
    theme = shinytheme("flatly"),
    tabPanel("Home", homepage, icon = icon("home")),
    navbarMenu(
      title = "miRNA", icon = icon("wave-square"),
      tabPanel("miRNA structure", mod_mirna_ui("mirna"), icon = icon("wave-square"),),
      tabPanel("miRNA target", mod_mirtar_ui("mirtar"), icon = icon("wave-square"),),
      tabPanel("miRNA expression", mod_mirnaexp_ui("mirnaexp"), icon = icon("chart-line"),),
      tabPanel("regulation element", mod_promoter_ui("mirnapro"), icon = icon("sliders-h"),)
    ),
    
    navbarMenu(
      title = "phasiRNA",  icon = icon("wave-square"),
      tabPanel("phasiRNA", mod_phasirna_ui("phasirna"),icon = icon("align-left"),),
      tabPanel("trigger", mod_trigger_ui("trigger"),icon = icon("align-left"),),
      tabPanel("phasiRNA expression", mod_phasirnaexp_ui("phasirnaexp"), icon = icon("chart-line"),)
    ),
    
    tabPanel("hcsiRNA", mod_hcsiRNA_ui("hcsirna"),icon = icon("align-left"),),
    
    tabPanel("Extraction", mod_extraction_ui("extraction"),icon = icon("search"),),
    
    navbarMenu(
      title = "Small tools", icon = icon("tools"),
      tabPanel("miRNA_target", mod_target_ui("mirtarget"),icon = icon("th", lib="glyphicon"),),
      tabPanel("Heatmap", mod_heatmap_ui("heatmap"),icon = icon("th", lib="glyphicon"),)
    ),
    
    navbarMenu(
      title = "Documentation", icon = icon("book"),
      tabPanel("Resources", mod_resource_ui("resource"),icon = icon("sourcetree"),),
      tabPanel("Tea Genome", mod_genome_ui("teagenome"),icon = icon("stumbleupon-circle"),),
      tabPanel("References", mod_reference_ui("references"),icon = icon("redo"),),
      tabPanel("FAQs", mod_faqs_ui("faqs"),icon = icon("question-circle"),)
    ),
    tabPanel("About", mod_about_ui("about"), icon = icon("info-circle")),
    footer = footerTagList
  )
)

server <- function(input, output, session) {
  callModule(mod_mirna_server, "mirna")
  callModule(mod_mirtar_server, "mirtar")
  callModule(mod_mirnaexp_server, "mirnaexp")
  callModule(mod_promoter_server, "mirnapro")
  
  callModule(mod_phasirna_server, "phasirna")
  callModule(mod_trigger_server, "trigger")
  callModule(mod_phasirnaexp_server, "phasirnaexp")
  
  callModule(mod_hcsiRNA_server, "hcsirna")
  
  callModule(mod_extraction_server, "extraction")
  
  callModule(mod_target_server, "mirtarget")
  callModule(mod_heatmap_server, "heatmap")
  
  
  callModule(mod_resource_server, "resource")
  #callModule(mod_search_server, "diversity")
  callModule(mod_download_server, "downloads")
  observeEvent(input$disconnect, {
    session$close()
  })
}

shinyApp(ui, server)