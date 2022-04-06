mod_download_ui <- function(id){
  ns <- NS(id)
  tagList(
    mainPanel(
      width = 12,
      h3("(1) The table of sRNA sample ID adn degradome ID from NCBI:"),
      downloadButton(ns("downloadresource"),label = "Download"),
      br(),
      h3("(2) The table of information of all miRNAs:"),
      downloadButton(ns("downloadmiRNAinfor"),label = "Download"),
      br(),
      h3("(3) The table of regulatory elements on miRNA promoter:"),
      downloadButton(ns("downloadcare"),label = "Download"),
      br(),
      h3("(4) The table of TF binding sites on miRNA promoter:"),
      downloadButton(ns("downloadtfbind"),label = "Download"),
      br(),
      br()
    )
  )
}

mod_download_server <- function(input, output, session) {
  ns <- session$ns
 # all_resource<-read.table("data/all_sample.txt",head=T,sep='\t')
    
  output$downloadresource <- downloadHandler(
    filename = function(){
      paste0("resources", ".txt")
    },
    content=function(file){
      write.table(allsample,file,row.names = TRUE, col.names = TRUE, quote=FALSE)
    }
  )
  
  output$downloadmiRNAinfor <- downloadHandler(
    filename = function(){
      paste0("mirna_table", ".txt")
    },
    content=function(file){
      write.table(mirna_info,file,row.names = TRUE, col.names = TRUE, quote=FALSE)
    }
  )
  output$downloadmiRNAinfor <- downloadHandler(
    filename = function(){
      paste0("mirna_table", ".txt")
    },
    content=function(file){
      write.table(mirna_info,file,row.names = TRUE, col.names = TRUE, quote=FALSE)
    }
  )
  
  output$downloadcare <- downloadHandler(
    filename = function(){
      paste0("plantcare", ".txt")
    },
    content=function(file){
      write.table(miRNAcare,file,row.names = TRUE, col.names = TRUE, quote=FALSE)
    }
  )
  
  output$downloadtfbind <- downloadHandler(
    filename = function(){
      paste0("tfbind", ".txt")
    },
    content=function(file){
      write.table(TFs,file,row.names = TRUE, col.names = TRUE, quote=FALSE)
    }
  )
}
