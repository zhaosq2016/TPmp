mod_promoter_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 3,
      h3("The main options of miRNA promoter:"),
      chooserInput(ns("mirnaidpro"), "Available IDs", "Selected IDs", MIRNA[2:638], MIRNA[1], size = 10, multiple = TRUE),
      
      h3("The main options:"),
      selectInput(ns("caretype"), "Please select regulatory element type:", choices = c("all Plantcare","Plantcare html","phytohormone responsive", "abiotic responsive","box","TFs")),
      br(),
      actionButton(ns("promoter_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      downloadButton(ns("downloadelement"),label = "Download element table"),
      br()
      ),
    
    mainPanel(
      h3("regulatory element on miRNA promoter"),
      withSpinner(DTOutput(ns("carehtml"))),
      br(),
      br(),
      br()
    )
  )
}


mod_promoter_server <- function(input, output, session){
  ns <- session$ns

  carehtml_table <- eventReactive(input$promoter_submit,{
    miRNAid_select <- input$mirnaidpro$selected
    if(input$caretype =="all Plantcare"){
      filter1<- miRNAcare[1,] 
      for(i in 1:length(miRNAid_select)){
        filter1<- rbind( filter1, dplyr::filter(miRNAcare, miRNAcare[,1] == miRNAid_select[i]))
      }
      filter1<-filter1[-1,]
      rownames(filter1) <- seq(1,nrow(filter1),1)
      filter1
    }else if(input$caretype =="Plantcare html"){
      url <- paste0("https://gitee.com/zhaoshiqi89/mp/raw/master/",miRNAid_select,"_plantcare.html")
      download.file(url,destfile = "www/mp/promoter_mir.html", quiet = TRUE)
      linktext <- paste0("<a href='mp/promoter_mir.html' target='blank' >",miRNAid_select,"</a>")
      data.frame(ID = miRNAid_select, PlantCARE = linktext)
    }else if(input$caretype =="phytohormone responsive"){
      filter1<- data.frame(ID= character(), position= numeric(), type= character()) 
      for(i in 1:length(miRNAid_select)){
        filter1<- rbind( filter1, dplyr::filter(phytohormone, phytohormone[,1] == miRNAid_select[i]))
      }
      filter1
    }else if(input$caretype == "abiotic responsive"){
      filter1<- data.frame(ID= character(), position= numeric(), type= character())
      for(i in 1:length(miRNAid_select)){
        filter1<- rbind( filter1, dplyr::filter(abiotic, abiotic[,1] == miRNAid_select[i]))
      }
      filter1
    }else if(input$caretype == "box"){
      filter1<- data.frame(ID= character(), type= character(), position= numeric())
      for(i in 1:length(miRNAid_select)){
        filter1<- rbind( filter1, dplyr::filter(care_box, care_box[,1] == miRNAid_select[i]))
      }
      filter1
    }else if(input$caretype == "TFs"){
      filter1<- TFs[1,]
      for(i in 1:length(miRNAid_select)){
        filter1<- rbind( filter1, dplyr::filter(TFs, TFs[,3] == miRNAid_select[i]))
      }
      filter1<-filter1[-1,]
      rownames(filter1) <- seq(1,nrow(filter1),1)
      filter1
    }
  }
  )
  
  output$carehtml <- renderDataTable({
    carehtml_table()
  }, escape = FALSE)
  
  output$downloadelement <- downloadHandler(
    filename = function(){
      paste0("element", '.txt')
    },
    content=function(file){
      write.table(carehtml_table(),file,row.names = FALSE, col.names = FALSE, quote=FALSE)
    }
  )
  
}
