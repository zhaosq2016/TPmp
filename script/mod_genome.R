mod_genome_ui <- function(id){
  ns <- NS(id)
  tagList(
    mainPanel(
      width = 12,
      h2("Information about Tea Genome (", tags$em("Camellia Sinensis"), ")"),
      h3("(1) The Tea Genome of Shuchazao (AAU):", tags$a(href="http://tpia.teaplant.org/", "  Click here!")),
      h3("(2) The Tea Genome of Shuchazao (TRI):", tags$a(href="https://github.com/JiedanChen/TeaGenomeData", "  Click here!")),
      h3("(3) The Tea Genome of Longjing43:",tags$a(href="https://bigd.big.ac.cn/gwh/Assembly/1086/show", "  Click here!")),
      h3("(4) The Tea Genome of Biyun:",tags$a(href="https://bigd.big.ac.cn/gwh/Assembly/8796/show", "  Click here!")),
      h3("(5) The Tea Genome of Yunkang10:",tags$a(href="www.plantkingdomgdb.com/tea_tree/", "  Click here!")),
      h3("(6) The Tea Genome of DASZ(an ancient tea tree):",tags$a(href="https://figshare.com/articles/journal_contribution/Assembly_and_annotation_of_DASZ_genome/12560462", "  Click here!")),
      h3("(7) The Tea Genome of TV-1:",tags$a(href="https://www.biorxiv.org/content/10.1101/762161v1.full", "  Click here!")),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br()
    )
  )
}