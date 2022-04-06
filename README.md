# TPmp

&nbsp;**TPmp** (microRNAs and phasiRNAs of tea plant ) is an easy-to-use, interactive web application with multiple analysis modules developed to visualize and explore miRNA and phasiRNA among 81 sRNA-seq samples  and 10 degradome-seq samples based on the data in NCBI. 

## Online

The app is hosted on Shinyapps.io here:  https://zhaosq89.shinyapps.io/tpmp for online use.

## Locally

### Requirements

* R: https://www.r-project.org v4.0.0+
* RStudio: https://rstudio.com/products/rstudio/download

### Run

1. Download this repository
    &nbsp;Code&rarr;Download ZIP
    
2. Unzip "TPmp-master.zip" file

3. Open the app.R file in the folder with Rstudio

4. Start the app

```
install.packages("shiny")
```
  then:

   Click Rstudio's Run App icon to run

   or
```
shiny::runApp("app.R")
```

> During the running process, a series of R packages will be installed, please be patient.

