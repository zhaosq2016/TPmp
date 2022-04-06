mod_reference_ui <- function(id){
  ns <- NS(id)
  tagList(
    mainPanel(
      width = 12,
      h2("Articles referenced in this research:"),
      h3("(1) ", tags$a(href="https://www.sciencedirect.com/science/article/pii/S2468014119301992?utm_medium=cpc&utm_source=TrendMD", "Comprehensive Characterization of miRNA and PHAS Loci in the Diploid Strawberry (Fragaria vesca) Genome")),
      h3("(2) ", tags$a(href="http://www.plantcell.org/content/25/5/1555", "MicroRNA Superfamilies Descended from miR390 and Their Roles in Secondary Small Interfering RNA Biogenesis in Eudicots")),
      br(),
      br(),
      h2("Softwares used in this research: "),
      h3("(1) ", tags$a(href="http://hannonlab.cshl.edu/fastx_toolkit", "Fastx-toolkit")),
      h3("(2) ", tags$a(href="https://www.yuque.com/books/share/169db1d9-bd85-4d92-b3d3-9ff6139ee29b?#", "sRNAminer")),
      h3("(3) ", tags$a(href="https://github.com/MikeAxtell/CleaveLand4", "CleaveLand4")),
      br(),
      br(),
      br()
    )
  )
}