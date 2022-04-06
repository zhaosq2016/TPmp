target_plot <- function(GFF_file, target_file, splice_number){

  GFF_info <- as_tibble(GFF_file)
  colnames(GFF_info) <- c("Chr", "info", "element", "start", "end", "info1", "info2", "info3", "ID")
  GFF_info <- GFF_info %>% dplyr::filter(element %in% c('exon', 'five_prime_UTR', 'three_prime_UTR'))
  GFF_info <- GFF_info[,-c(2,8)]
  GFF_info[,5] <- GFF_info[,4] - GFF_info[3] + 1
  gene_exon <- GFF_info %>% dplyr::filter(element == "exon")
  gene_exon <- gene_exon %>% select(-element)
  gene_exon <- arrange(gene_exon, start)
  gene_exon <- as.data.frame(gene_exon)
  
  gene_exon[,6] <- gsub("Parent=(.*);", "\\1", gene_exon[,6])
  
  gene_utr <- GFF_info %>% dplyr::filter(grepl("UTR",element))
  gene_utr <- gene_utr %>% select(-element)
  gene_utr <- arrange(gene_utr, start)
  gene_utr <- as.data.frame(gene_utr)
  
  target_info <- target_file
  
  if(str_detect(target_info[1,1], "Transcript") == TRUE){
    gene_name <- gsub("^.* Transcript: ", "", target_info[1,1])
  }else{
    gene_name <- gene_exon[1,6]
  }
  
  target_info[1,1] <- gsub(" Transcript.*", "", target_info[1,1])
  target_info[3,1] <- gsub(" Query.*", "", target_info[3,1])
  
  splice_number <- splice_number
  
  splice_site <- function(x,y){
    exons_data <- x
    exons <- exons_data[,4]
    exon_num <- nrow(exons_data)
    number <- y
    if(exons_data[1,5] == "+"){
      for(i in 1:exon_num){
        if(sum(exons[1:i]) < number){
          next
        }
        b <- i
        break
      }
      if(b == 1){
        site_num <- exons_data[b,2] + number
      }else{
        site_num <- exons_data[b,2] + (number - sum(exons[1:(b-1)]))
      }
    }else{
      for(i in exon_num:1){
        if(sum(exons[i:exon_num]) < number){
          next
        }
        b <- i
        break
      }
      if(b == exon_num){
        site_num <- exons_data[b,3] - number
      }else{
        site_num <- exons_data[b,3] - (number - sum(exons[(b+1):exon_num]))
      }
    }
    return(site_num)
  }
  
  if(!is.na(splice_number)){
    splice_sitenum <- splice_site(gene_exon, splice_number)
  }else{
    splice_sitenum <- 10
  }
  
  
  a <- (gene_exon[nrow(gene_exon),3]- gene_exon[1,2])/20
  b <- (gene_exon[nrow(gene_exon),3]- gene_exon[1,2])/2
  
  if(nrow(gene_utr) != 0){
    if(gene_exon[1,5] == "+"){
      p <- ggplot() + 
        geom_segment(aes(x=gene_exon[1,2],xend=gene_exon[nrow(gene_exon),3],y=0,yend=0),colour="black") +
        geom_rect(data=gene_exon,aes(xmin=start, xmax=end,ymin=-0.2,ymax=0.2),
                  colour="white", fill="orange") +
        geom_rect(data=gene_utr,aes(xmin=start, xmax=end,ymin=-0.2,ymax=0.2),
                  colour="white", fill="white") +
        geom_rect(data=gene_utr,aes(xmin=start, xmax=end,ymin=-0.17,ymax=0.17),
                  colour="blue", fill="blue") +
        ylim(c(-3,1)) + 
        geom_segment(aes(x=splice_sitenum,xend=splice_sitenum,y=-0.4,yend=-0.2),colour="black", arrow = arrow(length=unit(0.2, "cm"), type="closed"))+
        geom_segment(aes(x=splice_sitenum,xend=gene_exon[1,2]+a,y=-0.4,yend=-1.1),colour="black")+
        geom_segment(aes(x=splice_sitenum,xend=gene_exon[nrow(gene_exon),3]-a,y=-0.4,yend=-1.1),colour="black")+
        geom_text(aes(x=gene_exon[1,2]+a, y = -1.3,label=target_info[1,1]),family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=gene_exon[1,2]+a, y = -1.7,label=target_info[2,1]),family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=gene_exon[1,2]+a, y = -2.1,label=target_info[3,1]),family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=gene_exon[nrow(gene_exon),3]-b, y = 0.7,label=gene_name), size=10) +
        geom_text(aes(x=gene_exon[nrow(gene_exon),3]-b, y = -2.5,label=target_info[1,2]), size=10) +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
    }else{
      p <- ggplot(gene_exon) + 
        geom_segment(aes(x=-gene_exon[1,2],xend=-gene_exon[nrow(gene_exon),3],y=0,yend=0),colour="black") +
        geom_rect(data=gene_exon,aes(xmin=-start, xmax=-end,ymin=-0.2,ymax=0.2),
                  colour="white", fill="orange") +
        geom_rect(data=gene_utr,aes(xmin=-start, xmax=-end,ymin=-0.2,ymax=0.2),
                  colour="white", fill="white") +
        geom_rect(data=gene_utr,aes(xmin=-start, xmax=-end,ymin=-0.17,ymax=0.17),
                  colour="blue", fill="blue") +
        ylim(c(-3,1)) + 
        geom_segment(aes(x=-splice_sitenum,xend=-splice_sitenum,y=-0.4,yend=-0.2),colour="black", arrow = arrow(length=unit(0.2, "cm"), type="closed"))+
        geom_segment(aes(x=-splice_sitenum,xend=-gene_exon[1,2]-a,y=-0.4,yend=-1.1),colour="black")+
        geom_segment(aes(x=-splice_sitenum,xend=-gene_exon[nrow(gene_exon),3]+a,y=-0.4,yend=-1.1),colour="black")+
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+a, y = -1.3,label=target_info[1,1]), family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+a, y = -1.6,label=target_info[2,1]), family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+a, y = -1.9,label=target_info[3,1]), family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+b, y = 0.7,label=gene_name), size=10) +
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+b, y = -2.5,label=target_info[1,2]), size=10) +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
    }
  }else{
    if(gene_exon[1,5] == "+"){
      p <- ggplot() + 
        geom_segment(aes(x=gene_exon[1,2],xend=gene_exon[nrow(gene_exon),3],y=0,yend=0),colour="black") +
        geom_rect(data=gene_exon,aes(xmin=start, xmax=end,ymin=-0.2,ymax=0.2),
                  colour="white", fill="orange") +
        ylim(c(-3,1)) + 
        geom_segment(aes(x=splice_sitenum,xend=splice_sitenum,y=-0.4,yend=-0.2),colour="black", arrow = arrow(length=unit(0.2, "cm"), type="closed"))+
        geom_segment(aes(x=splice_sitenum,xend=gene_exon[1,2]+a,y=-0.4,yend=-1.1),colour="black")+
        geom_segment(aes(x=splice_sitenum,xend=gene_exon[nrow(gene_exon),3]-a,y=-0.4,yend=-1.1),colour="black")+
        geom_text(aes(x=gene_exon[1,2]+a, y = -1.3,label=target_info[1,1]),family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=gene_exon[1,2]+a, y = -1.7,label=target_info[2,1]),family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=gene_exon[1,2]+a, y = -2.1,label=target_info[3,1]),family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=gene_exon[nrow(gene_exon),3]-b, y = 0.7,label=gene_name), size=10) +
        geom_text(aes(x=gene_exon[nrow(gene_exon),3]-b, y = -2.5,label=target_info[1,2]), size=10) +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
    }else{
      p <- ggplot(gene_exon) + 
        geom_segment(aes(x=-gene_exon[1,2],xend=-gene_exon[nrow(gene_exon),3],y=0,yend=0),colour="black") +
        geom_rect(data=gene_exon,aes(xmin=-start, xmax=-end,ymin=-0.2,ymax=0.2),
                  colour="white", fill="orange") +
        ylim(c(-3,1)) + 
        geom_segment(aes(x=-splice_sitenum,xend=-splice_sitenum,y=-0.4,yend=-0.2),colour="black", arrow = arrow(length=unit(0.2, "cm"), type="closed"))+
        geom_segment(aes(x=-splice_sitenum,xend=-gene_exon[1,2]-a,y=-0.4,yend=-1.1),colour="black")+
        geom_segment(aes(x=-splice_sitenum,xend=-gene_exon[nrow(gene_exon),3]+a,y=-0.4,yend=-1.1),colour="black")+
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+a, y = -1.3,label=target_info[1,1]), family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+a, y = -1.6,label=target_info[2,1]), family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+a, y = -1.9,label=target_info[3,1]), family="simsun", hjust = 0, size=14) +
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+b, y = 0.7,label=gene_name), size=10) +
        geom_text(aes(x=-gene_exon[nrow(gene_exon),3]+b, y = -2.5,label=target_info[1,2]), size=10) +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
    }
  }
  showtext_auto()
  p
}
