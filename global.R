options(warn = -1)

if(!require(shiny)){
  install.packages("shiny")
}

if(!require(shinythemes)){
  install.packages("shinythemes")
}

if(!require(shinydisconnect)){
  install.packages("shinydisconnect")
}

if(!require(shinydashboard)){
  install.packages("shinydashboard")
}

if(!require(shinycssloaders)){
  install.packages("shinycssloaders")
}

if(!require(shinyWidgets)){
  install.packages("shinyWidgets")
}

if(!require(pdftools)){
  install.packages("pdftools")
}

if(!require(DT)){
  install.packages("DT")
}

if(!require(tibble)){
  install.packages("tibble")
}

if(!require(dplyr)){
  install.packages("dplyr")
}

if(!require(tidyverse)){
  install.packages("tidyverse")
}


if(!require(tableHTML)){
  install.packages("tableHTML")
}

if(!require(markdown)){
  install.packages("markdown")
}

if(!require(pheatmap)){
  install.packages("pheatmap")
}

if(!require(ggplot2)){
  install.packages("ggplot2")
}

if(!require(showtext)){
  install.packages("showtext")
}

if(!require(sysfonts)){
  install.packages("sysfonts")
}

if(!require(tidyr)){
  install.packages("tidyr")
}

if(!require(stringr)){
  install.packages("stringr")
}

library(stringr)
library(tidyr)
library(sysfonts)
library(dplyr)
library(tibble)
library(ggplot2)
library(showtext)
library(shiny)
library(shinythemes)
library(shinydisconnect)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(pdftools)
library(markdown)
library(DT)
library(tidyverse)
library(tableHTML)
library(pheatmap)
library(ggplot2)

source("script/RNAPlot2.R")
source("script/chooser.R")
source("script/mod_mirna.R")
source("script/mod_mirtar.R")
source("script/mod_mirnaexp.R")
source("script/mod_heatmap.R")
source("script/mod_phasirna.R")
source("script/mod_trigger.R")
source("script/mod_phasirnaexp.R")
source("script/mod_hcsiRNA.R")
source("script/mod_extraction.R")
source("script/boxformate.R")
source("script/homepage.R")
source("script/mod_about.R")
source("script/mod_reference.R")
source("script/mod_resources.R")
source("script/mod_faqs.R")
source("script/mod_genome.R")
source("script/mod_download.R")
source("script/mod_promoter.R")
source("script/mod_target.R")
source("script/target_plot.R")

footerTagList <- list(
  tags$footer(id = "myFooter",
              shiny::includeHTML("www/md/footer.html")
  )
)

font_add("simsun",regular = "www/font/simsunb.ttf")

allsample <- read.table("data/all_sample.txt",head=T,sep='\t')

mirna_info <- read.table("data/Css_miRNA_info.txt", header = T, sep='\t')
mirna_list <- mirna_info[,1]
rownames(mirna_info) <- mirna_info[,1]

miRNAtarget <- read.table("data/miRNA_target.txt", head=F, sep='\t', quot="")
mirna_list2 <- unique(miRNAtarget[,2])
GFF_info <- read.table("data/CSS.gff3", header = F, sep='\t')

allmiRNAexp <- read.table("data/Css_miRNA_exp.txt", head=T, sep='\t', row.names = 1)
mirna_list3 <- row.names(allmiRNAexp)

phasiRNA21 <- read.table("data/merged.21.PHAS.list.filter.txt", head=T, sep='\t', row.names = 1)
phasiRNA24 <- read.table("data/merged.24.PHAS.list.filter.txt", head=T, sep='\t', row.names = 1)
phasirna21_list <- row.names(phasiRNA21)
phasirna24_list <- row.names(phasiRNA24)

phasiRNA21exp <- read.table("data/merged.21.PHAS.exp.txt", head=T, sep='\t', row.names = 1)
phasiRNA24exp <- read.table("data/merged.24.PHAS.exp.txt", head=T, sep='\t', row.names = 1)
phasiRNAsample <- colnames(phasiRNA21exp)

trigger_21 <- read.table("data/phasiRNA_21nt_trigger.result.txt", head = T, sep = '\t', row.names = 1)
trigger_24 <- read.table("data/phasiRNA_24nt_trigger.result.txt", head = T, sep = '\t', row.names = 1)

hcsirna <- read.table("data/hcsiRNA.txt", head = T, sep='\t', row.names = 1)

allresult <- unique(c(mirna_info[,2], phasiRNA21[,1], phasiRNA24[,1], hcsirna[,1]))

miRNAcare <- read.table("data/miRNA_CARE.txt", head=T, sep='\t')
MIRNA <- as.character(unique(miRNAcare[,1]))

phytohormone<-read.table("data/phytohormone.txt",head=T,sep='\t')
abiotic<-read.table("data/abiotic.txt",head=T,sep='\t')
care_box<-read.table("data/box.txt",head=T,sep='\t')
TFs <- read.table("data/TFs.txt",head=T,sep='\t')
