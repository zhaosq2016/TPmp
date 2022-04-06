loadCoords<-function(filename){
  data<-read.table(file=filename,sep=",",header=FALSE)
  if ( dim(data)[2]==6 ) {
    names(data)<-c("id","x","y","seq","num","bound")
  } else {
    if ( dim(data)[2]==7 ){
      names(data)<-c("id","x","y","seq","num","bound","energy")
    } else {
      data=NULL
      print(paste("The file ",filename," does not have the correct number of columns",sep=""))
      print(paste("Expected 6 or 7, found ",dim(data)[2],sep=""))
    }
  }
  data
}

alignCoord<-function(data,n1,n2,x,y){
  ids<-data$id
  ids<-table(ids)
  ids<-as.numeric(names(ids))
  
  ## Loop through all ids in the dataset ##
  for(i in ids){
    tmp<-NULL
    #print(i)
    tmp<-data[data$id==i,]
    x1<-tmp$x[tmp$num==n1]
    x2<-tmp$x[tmp$num==n2]
    
    y1<-tmp$y[tmp$num==n1]
    y2<-tmp$y[tmp$num==n1]
    
    ## Translate x1,y1 to x,y ##
    tmp$x<-tmp$x-(x1-x)
    tmp$y<-tmp$y-(y1-y)
    dy<-tmp$y[tmp$num==n2]-tmp$y[tmp$num==n1]
    dx<-tmp$x[tmp$num==n2]-tmp$x[tmp$num==n1]
    ang<-atan(abs(dy)/abs(dx))
    if(dx<0 & dy < 0){
      #print("2")
      ang<-ang+pi
    }else{
      if(dx< 0 & dy >0){
        #print("1")
        ang<-pi-ang
      }else{
        if(dx>0 & dy <0){
          #print("3")
          ang<-2*pi-ang
        }
      }
    }
    
    ## Rotate x2,y2 around x1,y1 until y1==y2
    for(j in 1:length(tmp$x)){
      pt<-rotateS(tmp$x[j],tmp$y[j],x,y,-ang)
      tmp$x[j]<-pt[1]
      tmp$y[j]<-pt[2]
    }
    data[data$id==i,]<-tmp
  }
  data
}

RNAPlot2<-function(data,ranges=0,add=FALSE,hl=NULL,seqcols=NULL,seqTF=FALSE,labTF=FALSE,nt=FALSE,dp=0.5,modspec=FALSE,modp=NULL,mod=NULL,modcol=NULL,tsize=0.5,main="",pointSize=2,lineWd=2){
  l=1
  if(seqTF){
    s<-data$seq
    sequence<-paste(s,collapse="")
  }
  if(is.null(dim(ranges)) & length(seqcols)<1){l=0}
  if(is.null(dim(ranges))){
    ranges<-NULL
    ranges$min[1]<-0
    ranges$max[1]<-dim(data)[1]
    ranges$desc[1]<-c("RNA Molecule")
    ranges$col[1]<--1
    ranges<-as.data.frame(ranges)
  }
  
  s2<-dim(data)[1]
  if(length(hl)<1){}else{
    hlranges<-NULL
    for(i in 1:length(hl)){
      t<-strsplit(hl[i],"")
      s1<-length(t[[1]])
      
      groups<-NULL
      if(s2<=s1){
        groups[[1]]<-as.character(data$seq)
      }else{
        for(j in 1:(s2-s1+1)){
          groups[[j]]<-as.character(data$seq[j:(j+(s1-1))])
        }
      }
      groups<-lapply(groups,paste,collapse="")
      groups<-unlist(groups)
      ind<-c(1:length(groups))[groups==hl[i]]
      ind<-data$num[ind]
      
      hlranges$min<-c(hlranges$min,ind)
      hlranges$max<-c(hlranges$max,(ind+(s1-1)))
      hlranges$desc<-c(hlranges$desc,rep(hl[i],length(ind)))
      hlranges$col=c(hlranges$col,rep(seqcols[i],length(ind)))
      
    }
    
    hlranges<-as.data.frame(hlranges)
    ranges<-rbind(ranges,hlranges)
  }
  s1<-max(c(data$x,data$y))
  s2<-min(c(data$x,data$y))
  if(add==FALSE){
    plot(data$x,data$y,type="n",xlab="",ylab="",axes=FALSE,xlim=c(s2,s1),ylim=c(s2-10,s1))
    for(i in data$num){
      if(data$bound[data$num==i]!=-1){	
        lines(c(data$x[data$num==i],data$x[data$num==data$bound[data$num==i]]),c(data$y[data$num==i],data$y[data$num==(data$bound[data$num==i])]),lwd=lineWd)	
      }	
      
    }
  }else{
    if(add==TRUE){
      points(data$x,data$y,type="n")
      for(i in data$num){
        if(data$bound[data$num==i]!=-1){	
          lines(c(data$x[data$num==i],data$x[data$num==data$bound[data$num==i]]),c(data$y[data$num==i],data$y[data$num==(data$bound[data$num==i])]),lwd=lineWd)	
        }	
        
      }
    }
  }

  x<-data$x
  y<-data$y
  if(nt){}else{
    lines(x,y,lw=lineWd)
    points(x,y,pch=16,cex=pointSize)
  }
  if(l==1){
    if(labTF){
      legend(s2[1],s2[1]-2,
             legend=ranges$desc[ranges$col!=0],
             col=abs(ranges$col[ranges$col!=0]),bty="n",
             lty=1			
      )
    }
  }
  
  if(nt==TRUE){
    x1<-data$x[1]
    x2<-data$x[2]
    y1<-data$y[1]
    y2<-data$y[2]
    vx<-(x1-x2)
    vy<-(y1-y2)
    dist<-sqrt(vx^2+vy^2)
    vx<-vx/dist
    vy<-vy/dist
    
    text((data$x[1]+(vx*dp)),(data$y[1]+(vy*dp)),"5'",font=2,cex=tsize)
    
    x3<-(-1*(y1-y2))
    y3<-(x1-x2)
    vx<-x3/dist
    vy<-y3/dist
    text((data$x[1]+(vx*dp)),(data$y[1]+(vy*dp)),data$seq[1],cex=tsize,font=2)
    
    t1<-length(data$x)
    t2<-(length(data$x)-1)
    x1<-data$x[t1]
    x2<-data$x[t2]
    y1<-data$y[t1]
    y2<-data$y[t2]
    vx<-(x1-x2)
    vy<-(y1-y2)
    dist<-sqrt(vx^2+vy^2)
    vx<-vx/dist
    vy<-vy/dist
    
    text((data$x[t1]+(vx*dp)),(data$y[t1]+(vy*dp)),"3'",cex=tsize,font=2)
    
    x3<-(-1*(y1-y2))
    y3<-(x1-x2)
    vx<-x3/dist
    vy<-y3/dist
    text((data$x[t1]-(vx*dp)),(data$y[t1]-(vy*dp)),data$seq[t1],font=2,cex=tsize)
    
    
    for(i in 2:(length(data$x)-1)){
      x1<-data$x[(i-1)]
      x2<-data$x[(i+1)]
      y1<-data$y[(i-1)]
      y2<-data$y[(i+1)]
      
      x3<-(-1*(y1-y2))
      y3<-(x1-x2)
      dist<-sqrt(x3^2+y3^2)
      x3<-x3/dist
      y3<-y3/dist
      
      text(data$x[i]+(x3*dp),data$y[i]+(y3*dp),data$seq[i],font=2,cex=tsize)
    }
    
    
    
  }
  if(modspec==TRUE){
    for(i in 1:length(modp)){
      points(data$x[modp[i]],data$y[modp[i]],pch=mod[i],col=modcol[i])
      
    }
    
  }
  
  
  for(i in 1:dim(ranges)[1]){
    
    keep<-data$num >= ranges$min[i] & data$num <= ranges$max[i]
    if(ranges$col[i]!=0){
      if(ranges$col[i]==-1){
        
      }else{	
        if(ranges$col[i]==3){
          points(data$x[keep],data$y[keep],col=ranges$col[i],pch="x")
        }else{
          points(data$x[keep],data$y[keep],col=ranges$col[i],pch=19)
        }
      }
    }
    
    
  }
  title(main)
}



bplfile<-function(dat,name){
  seq<-paste(dat$seq,collapse="")
  fileo=paste(">",name,sep="")
  fileo=paste(fileo,"\n",seq,"\n",sep="")
  fileo=paste(fileo,"\n$ list",sep="")
  for(i in 1:length(dat$num)){
    if(dat$bound[i]==-1){}else{
      fileo=paste(fileo,"\n",(dat$num[i]+1)," ",(dat$bound[i]+1))}
  }
  write(fileo,file=name)
}


transformFold<-function(dat,x0,y0,ang){
  x<-dat$x
  y<-dat$y
  new<-rotateV(x,y,x0,y0,ang)
  dat$x<-new$x
  dat$y<-new$y
  dat
}
aptPlotCT<-function(file,ranges=0,add=FALSE,hl=NULL,seqcols=NULL,seqTF=FALSE,labTF=FALSE,nt=FALSE,dp=0.5,modspec=FALSE,modp=NULL,mod=NULL,modcol=NULL,tsize=0.5,main="",pseudoTF=FALSE,pseudo_nums=NULL,ticks=NULL,ticksTF=FALSE){
  
  dat<-loadCt(file)
  pseudo<-pseudoKnot(dat)
  data<-ct2coord(pseudo[[2]])
  names(data)[1]="num"	
  data<-as.data.frame(data)	
  l=1
  if(seqTF){
    s<-dat$seq
    sequence<-paste(s,collapse="")
  }
  if(ranges==0 & length(seqcols)<1){l=0}
  if(ranges==0){
    ranges<-NULL
    ranges$min[1]<-0
    ranges$max[1]<-dim(data)[1]
    ranges$desc[1]<-c("Aptamer")
    ranges$col[1]<--1
    ranges<-as.data.frame(ranges)
  }
  s2<-dim(data)[1]
  if(length(hl)<1){print("no sequences to match")}else{
    hlranges<-NULL
    for(i in 1:length(hl)){
      t<-strsplit(hl[i],"")
      s1<-length(t[[1]])
      groups<-NULL
      if(s2<=s1){
        groups[[1]]<-as.character(data$seq)
      }else{
        for(j in 1:(s2-s1)){
          groups[[j]]<-as.character(data$seq[j:(j+(s1-1))])
        }
      }
      groups<-lapply(groups,paste,collapse="")
      groups<-unlist(groups)
      ind<-c(1:length(groups))[groups==hl[i]]
      ind<-data$num[ind]
      hlranges$min<-c(hlranges$min,ind)
      hlranges$max<-c(hlranges$max,(ind+(s1-1)))
      hlranges$desc<-c(hlranges$desc,rep(hl[i],length(ind)))
      hlranges$col=c(hlranges$col,rep(seqcols[i],length(ind)))
      
    }
    hlranges<-as.data.frame(hlranges)
    ranges<-rbind(ranges,hlranges)
  }
  
  s1<-max(c(data$x,data$y))
  s1x<-max(c(data$x))+3
  s1y<-max(c(data$y))+3
  s2<-min(c(data$x,data$y))
  s2x<-min(c(data$x))-3
  s2y<-min(c(data$y))-3
  
  if(add==FALSE){
    plot(data$x,data$y,type="n",xlab="",ylab="",axes=FALSE,xlim=c(s2x,s1x),ylim=c(s2y,s1y))
  }else{
    if(add==TRUE){
      points(data$x,data$y,type="n")
    }
  }
  x<-data$x
  y<-data$y
  if(nt){}else{
    lines(x,y,lw=4)}
  if(l==1){
    if(labTF){
      legend(s2[1],s1[1],
             legend=ranges$desc[ranges$col!=0],
             col=abs(ranges$col[ranges$col!=0]),bty="n",
             lty=1			
      )
    }
  }
  
  if(nt==TRUE){
    x1<-data$x[1]
    x2<-data$x[2]
    y1<-data$y[1]
    y2<-data$y[2]
    vx<-(x1-x2)
    vy<-(y1-y2)
    dist<-sqrt(vx^2+vy^2)
    vx<-vx/dist
    vy<-vy/dist
    text((data$x[1]+(vx*dp)),(data$y[1]+(vy*dp)),"5'",cex=tsize)
    
    x3<-(-1*(y1-y2))
    y3<-(x1-x2)
    vx<-x3/dist
    vy<-y3/dist
    text((data$x[1]+(vx*dp)),(data$y[1]+(vy*dp)),data$seq[1],cex=tsize)
    
    t1<-length(data$x)
    t2<-(length(data$x)-1)
    x1<-data$x[t1]
    x2<-data$x[t2]
    y1<-data$y[t1]
    y2<-data$y[t2]
    vx<-(x1-x2)
    vy<-(y1-y2)
    dist<-sqrt(vx^2+vy^2)
    vx<-vx/dist
    vy<-vy/dist
    
    text((data$x[t1]+(vx*dp)),(data$y[t1]+(vy*dp)),"3'",cex=tsize)
    
    x3<-(-1*(y1-y2))
    y3<-(x1-x2)
    vx<-x3/dist
    vy<-y3/dist
    text((data$x[t1]-(vx*dp)),(data$y[t1]-(vy*dp)),data$seq[t1],cex=tsize)
    
    if(modspec==TRUE){
      for(i in 1:length(modp)){
        if(mod[i]==1){
        }else{
          if(mod[i]==2){
          }else{
            if(mod[i]==3){
            }else{
              if(mod[i]==4){
                points(data$x[modp[i]],data$y[modp[i]],pch=1,col=modcol[i])
              }
            }
          }
        }
        
      }
      
    }
    
    
    for(i in 2:(length(data$x)-1)){
      x1<-data$x[(i-1)]
      x2<-data$x[(i+1)]
      y1<-data$y[(i-1)]
      y2<-data$y[(i+1)]
      ## Rotate 90 degrees and create unit vector ##
      x3<-(-1*(y1-y2))
      y3<-(x1-x2)
      dist<-sqrt(x3^2+y3^2)
      x3<-x3/dist
      y3<-y3/dist
      ## Add text to data$seq[i] ###
      text(data$x[i]+(x3*dp),data$y[i]+(y3*dp),data$seq[i],cex=tsize)
    }
    
    
    
  }
  
  if(ticksTF){
    for(i in 2:(length(data$x)-1)){
      k<-data$num[i]==ticks
      t<-table(k)
      if(length(t)>1){
        x1<-data$x[(i-1)]
        x2<-data$x[(i+1)]
        y1<-data$y[(i-1)]
        y2<-data$y[(i+1)]
        
        x3<-(-1*(y1-y2))
        y3<-(x1-x2)
        dist<-sqrt(x3^2+y3^2)
        x3<-x3/dist
        y3<-y3/dist
        
        points(c(data$x[i],(data$x[i]+x3*dp)),c(data$y[i],(data$y[i]+y3*dp)),type="l")
        dp2<-dp*1.8
        text(data$x[i]+(x3*dp2),data$y[i]+(y3*dp2),data$pos[i],cex=tsize)
      }
    }
    
  }
  if(pseudoTF){
    
    for(i in 2:(length(data$x)-1)){
      k<-data$num[i]==pseudo_nums
      t<-table(k)
      if(length(t)>1){
        x1<-data$x[(i-1)]
        x2<-data$x[(i+1)]
        y1<-data$y[(i-1)]
        y2<-data$y[(i+1)]
        
        x3<-(-1*(y1-y2))
        y3<-(x1-x2)
        dist<-sqrt(x3^2+y3^2)
        x3<-x3/dist
        y3<-y3/dist
        
        text(data$x[i]+(x3*dp),data$y[i]+(y3*dp),data$seq[i],cex=tsize)
        print(c(data$x[i]+(x3*dp),data$y[i]+(y3*dp)))
        
      }
    }
    
  }
  
  
  for(i in 1:dim(ranges)[1]){
    keep<-data$num >= ranges$min[i] & data$num <= ranges$max[i]
    if(ranges$col[i]!=0){
      if(ranges$col[i]==-1){
      }else{	
        if(ranges$col[i]==3){
          points(data$x[keep],data$y[keep],col=ranges$col[i],pch="x")
        }else{
          points(data$x[keep],data$y[keep],col=ranges$col[i],pch=19)
        }
      }
    }
    for(i in data$num){
      if(data$bound[data$num==i]!=-1){	
        lines(c(data$x[data$num==i],data$x[data$num==data$bound[data$num==i]]),c(data$y[data$num==i],data$y[data$num==(data$bound[data$num==i])]))	
        
      }	
    }
  }
  title(main)
}

pseudoKnot<-function(ctDat){
  
  pseudo<-NULL
  pk<-1
  for(i in ctDat$pos){
    bound=ctDat$bound[ctDat$pos==i]
    if(bound!=0 & i<bound){
      for(j in i:bound){
        nbound=ctDat$bound[ctDat$pos==j]
        if(nbound!=0){
          if(nbound>bound | nbound <i){
            input$bound[input$pos==j]=0
            input$bound[input$pos==nbound]=0
            pseudo$pos[pk]=j
            pseudo$bound[pk]=nbound
            pk<-pk+1
          }
        }	
      }
    }
  }
  out<-NULL
  pseudo<-as.data.frame(pseudo)
  psk<-unique(pseudo$pos)
  outpk<-NULL
  sk<-1
  for(pk in psk){
    outpk$pos[sk]=pk
    bound<-pseudo$bound[pseudo$pos==pk]
    bound<-bound[1]
    outpk$bound[sk]=bound
    sk<-sk+1
  }
  
  if(length(psk)>0){
    
    stems<-NULL
    d1<-diff(outpk$bound)
    d2<-diff(outpk$pos)
    ind<-1:length(d1)
    r1<-ind[d1!=-1]
    r2<-ind[d2!=1]
    r<-c(r1,r2)
    r<-unique(r)
    r<-sort(r)
    s<-1
    for(i in 1:length(r)){
      tmp<-NULL
      if(i==1){
        tmp$pos<-outpk$pos[1:r[1]]
        tmp$bound<-outpk$bound[1:r[1]]
      }else{
        tmp$pos<-outpk$pos[(r[(i-1)]+1):r[i]]
        tmp$bound<-outpk$bound[(r[(i-1)]+1):r[i]]
      }
      stems[[s]]<-tmp
      s<-s+1
    }
    l<-length(r)
    lp<-length(outpk$pos)
    tmp$pos<-outpk$pos[(r[l]+1):lp]
    tmp$bound<-outpk$bound[(r[l]+1):lp]
    stems[[s]]<-tmp
    
    compm<-NULL
    m1t=1
    for(i in 1:s){
      for(j in 1:s){
        comp<-0
        st1<-stems[[i]]
        st2<-stems[[j]]
        mb1<-min(c(st1$bound))
        mp1<-min(c(st1$pos))
        mb2<-min(c(st2$bound))
        mp2<-min(c(st2$pos))
        mm2<-min(c(mp2,mb2))
        mx2<-max(c(mp2,mb2))
        
        if(mp1 < mx2 & mp1 > mm2 & mb1 < mm2 ){
          comp<-1
        }else{
          if(mp1 < mx2 & mp1 > mm2 & mb1 > mx2 ){
            comp<-1
          }else{
            if(mb1 < mx2 & mb1 > mm2 & mp1 > mx2 ){
              comp<-1
            }else{
              if(mb1 < mx2 & mb1 > mm2 & mp1 < mm2 ){
                comp<-1
              }else{
                
                
              }
            }
          }
        }
        compm$i[m1t]=i
        compm$j[m1t]=j
        compm$c[m1t]=comp
        m1t<-m1t+1
        
      }
    }
    comp<-as.data.frame(compm)
    
    mySets<-NULL
    l<-c(1:length(stems))
    
    mys<-1
    while(length(l)>0){
      s1<-l[1]
      s<-comp[comp$c==0 & comp$i==s1,]
      
      v<-s$j[1]
      for(j in s$j){
        v<-c(v,j)
        l<-l[l!=j]
      }
      mySets[[mys]]=v
      mys<-mys+1
    }
    
    lst<-NULL
    for(i in 1:length(stems)){
      tmp<-stems[[i]]
      tmp<-tmp$pos
      lst[i]<-length(tmp)
    }
    lst2<-NULL
    for(i in 1:length(mySets)){
      tmp<-mySets[[i]]
      tmp<-tmp[2:length(tmp)]
      lst2[i]<-sum(lst[tmp])
    }
    
    largest<-c(1:length(lst2))[lst2==max(lst2)]
    largest<-mySets[[largest]]
    largest<-largest[2:length(largest)]
    
    addBack<-NULL
    tmp<-stems[[largest[1]]]
    addBack$pos<-tmp$pos
    addBack$bound<-tmp$bound
    if(length(largest)>1){
      for(i in largest[2:length(largest)]){
        tmp<-stems[[i]]
        addBack$pos<-c(addBack$pos,tmp$pos)
        addBack$bound<-c(addBack$bound,tmp$bound)
      }
    }
    for(pos in addBack$pos){
      bnd<-addBack$bound[addBack$pos==pos]
      input$bound[input$pos==pos]=bnd
      input$bound[input$pos==bnd]=pos
    }
    addBack<-as.data.frame(addBack)
    outpk<-as.data.frame(outpk)
    
    for(pos in addBack$pos){
      outpk<-outpk[outpk$pos!=pos,]
    }
  }else{
    outpk=NULL
    
  }
  out[[1]]<-outpk
  out[[2]]<-input
  out
}

loadCt<-function(file){
  input<-read.table(file=file,skip=1)
  names(input)=c("pos","seq","before","after","bound","pos2")
  input$seq<-as.character(input$seq)
  input
}
makeCt<-function(struct,seq){
  nts<-strsplit(seq,"")
  nts<-nts[[1]]
  stc<-strsplit(struct,"")
  stc<-stc[[1]]
  out<-NULL
  for(i in 1:length(stc)){
    out$pos[i]=i
    out$before[i]=(i-1)
    out$after[i]=(i+1)
    out$seq[i]=nts[i]
    out$pos2[i]=i
    if(stc[i]=="."){
      out$bound[i]=0
    }else{
      if(stc[i]==")"){
        b<-backward(stc,i)
        out$bound[i]=b
      }else{
        if(stc[i]=="("){
          b<-forward(stc,i)
          out$bound[i]=b
        }
      }
    }
    
  }
  
  out<-as.data.frame(out)
  out
}

forward<-function(stc,i){
  k<-1
  go<-1
  out<-0
  if(stc[i]=="("){
    i<-i+1
    for(j in i:length(stc)){
      if(go==1){
        if(stc[j]==")"){
          k<-k-1
          if(k==0){
            out=j
            go=0
          }
        }else{
          if(stc[j]=="("){
            k=k+1
          }
          
        }
      }
    }
  }
  out
}
backward<-function(stc,i){
  k<-1
  go<-1
  out<-0
  if(stc[i]==")"){
    i<-i-1
    for(j in i:1){
      if(go==1){
        if(stc[j]=="("){
          k<-k-1
          if(k==0){
            out=j
            go=0
          }
        }else{
          if(stc[j]==")"){
            k=k+1
          }
          
        }
      }
    }
  }
  out
}


ct2knet<-function(file,ind=0){
  dat<-loadCt(file)
  out<-""
  seq<-dat$seq[(ind+1):length(dat$seq)]
  seq<-paste(seq,collapse="")
  out=paste(">",file,"\n",sep="")
  out=paste(out,seq,"\n\n",sep="")
  out=paste(out,"$ list\n",sep="")
  for(i in 1:dim(dat)[1]){
    if(dat$bound[i]!=0){
      out=paste(out," ",(dat$pos[i]-ind),"   ",(dat$bound[i]-ind),"\n",sep="")
    }
  }
  out
}


ct2coord<-function(input){
  group=0;
  firstrun=1
  mnp<-min(input$pos)
  map<-max(input$pos)
  a1<-map+1
  a2<-map+2
  b1<-mnp-1
  b2<-mnp-2
  
  new1<-NULL
  new1$pos[1]=a1
  new1$before[1]=(a1-1)
  new1$after[1]=(a1+1)
  new1$seq[1]="A"
  new1$pos2[1]=a1
  new1$bound[1]=b1
  
  new2<-NULL
  new2$pos[1]=a2
  new2$before[1]=(a2-1)
  new2$after[1]=(a2+1)
  new2$seq[1]="A"
  new2$pos2[1]=a2
  new2$bound[1]=b2
  
  new3<-NULL
  new3$pos[1]=b1
  new3$before[1]=(b1-1)
  new3$after[1]=(b1+1)
  new3$seq[1]="A"
  new3$pos2[1]=b1
  new3$bound[1]=a1
  
  new4<-NULL
  new4$pos[1]=b2
  new4$before[1]=(b2-1)
  new4$after[1]=(b2+1)
  new4$seq[1]="A"
  new4$pos2[1]=b2
  new4$bound[1]=a2
  
  input<-rbind(new1,new2,input,new3,new4)
  input<-as.data.frame(input)
  
  mnp<-min(input$pos)
  
  output<-NULL
  nextNT=input[input$pos==mnp,]
  if(nextNT$bound!=0){
    
    output$pos[1]=mnp
    output$x[1]=0
    output$y[1]=0
    output$pos[2]=nextNT$bound
    output$x[2]=0
    output$y[2]=1
  }
  stems<-NULL
  stems[[1]]<-c(output$pos[1],output$pos[2])
  
  j<-3
  prev<-mnp
  npos<-input$after[input$pos==mnp]
  nextNT=input[input$pos==npos,]
  mp<-max(input$pos)
  while(length(stems)>0){
    newstems<-NULL
    newloops<-NULL
    ns<-1
    nl<-1
    for(i in 1:length(stems)){
      s1<-stems[[i]]
      
      p1<-s1[1]
      p2<-s1[2]
      if(firstrun==1){
        x3=-1
        y3=0
        firstrun=0
      }else{
        
        prev<-input$before[input$pos==p1]
        x3<-output$x[output$pos==prev]
        y3<-output$y[output$pos==prev]
      }
      x1<-output$x[output$pos==p1]
      y1<-output$y[output$pos==p1]
      x2<-output$x[output$pos==p2]
      y2<-output$y[output$pos==p2]
      
      sout<-stemCords(input,p1,p2,x1,y1,x2,y2,x3,y3)
      sdat<-sout[[1]]
      sdat<-as.data.frame(sdat)
      sdat<-sdat[sdat$pos!=p1,]
      sdat<-sdat[sdat$pos!=p2,]
      l<-dim(sdat)[1]
      if(l>0){
        s<-length(output$pos)
        s<-s+1
        for(sd in 1:length(sdat$pos)){
          output$pos[s]<-sdat$pos[sd]
          output$x[s]<-sdat$x[sd]
          output$y[s]<-sdat$y[sd]
          s<-s+1
        }
      }
      newloops[[nl]]=sout[[2]]
      nl<-nl+1
    }

    newstems<-NULL
    for( i in 1:length(newloops)){
      lps<-newloops[[i]]
      
      lp1<-loopLength(input,lps[1])
      
      if(length(lp1)>1){
        nstm<-lp1[[2]]
        for(nt in 1:length(nstm)){
          newstems[[ns]]<-nstm[[nt]]
          ns<-ns+1
        }
      }
      pp<-input$before[input$pos==lps[1]]
      if(output$x[output$pos==lps[1]] < output$x[output$pos==pp]){
        p3=1
      }else{
        p3=1
      }
      
      lout<-genCords(lp1,lps[1],lps[2],output,p3)
      pt<-lp1[[1]]
      tout<-lout
      lpl<-length(pt$pos)
      pt$pos=pt$pos[lpl:1]
      lout$pos=pt$pos
      lout<-as.data.frame(lout)
      lout<-lout[lout$pos!=lps[1],]
      lout<-lout[lout$pos!=lps[2],]
      s<-length(output$pos)
      s<-s+1
      for(lo in 1:length(lout$pos)){
        output$pos[s]<-lout$pos[lo]
        output$x[s]<-lout$x[lo]
        output$y[s]<-lout$y[lo]
        output$group[s]=group
        s<-s+1
      }
      group=group+1
    }
    stems<-newstems
    
  }		
  my<-min(output$y)-5
  xy<-max(output$y)+5
  mx<-min(output$x)-5
  xx<-max(output$x)+5
  
  output<-as.data.frame(output)
  input<-as.data.frame(input)
  all<-merge(output,input,by="pos")
  names(all)[1]="num"
  
  rmax<-dim(all)[1]
  rmax<-rmax-2
  all<-all[3:rmax,]
  return(all)
}

stemCords<-function(input,p1,p2,x1,y1,x2,y2,x3,y3){
  
  output<-NULL
  output$pos[1]=p1
  output$x[1]=x1
  output$y[1]=y1
  output$pos[2]=p2
  output$x[2]=x2
  output$y[2]=y2
  ends<- -1
  
  vperp<-NULL
  vperp$x=x3-x1
  vperp$y=y3-y1
  
  v<-rotateS((x2-x1),(y2-y1),0,0,(pi/2))
  m<-sqrt(v[1]^2+v[2]^2)+sqrt(vperp$x^2+vperp$y^2)
  ang<-acos((v[1]*vperp$x+v[2]*vperp$y)/m)
  
  if(ang<pi/2){
    v<-rotateS((x2-x1),(y2-y1),0,0,(-1*pi/2))
  }
  vect<-NULL
  vect$x=v[1]
  vect$y=v[2]
  nv<-sqrt((vect$x)^2+(vect$y)^2)
  vect$x<-vect$x/nv
  vect$y<-vect$y/nv
  prev1=p1
  j<-3
  mf<-length(input$pos)+5
  mp<-max(input$pos)
  while(j<mf){
    npos=input$after[input$pos==prev1]
    if(npos<(mp+1)){
      nextNT=input[input$pos==npos,]
      if(dim(nextNT)[1] < 1){
        j=mf+5
      }else{
        if(nextNT$bound==0){
          
          bound<-input$bound[input$pos==prev1]
          ends<-c(prev1,bound)
          j=mf+5
        }else{
          
          pbpos=input$bound[input$pos==prev1]
          pbafter=input$before[input$pos==pbpos]
          bpos=input$bound[input$pos==npos]
          
          if(pbafter==bpos){	
            
            output$pos[j]=npos
            output$x[j]=output$x[output$pos==prev1]+vect$x[1]
            output$y[j]=output$y[output$pos==prev1]+vect$y[1]
            j<-j+1
            if(j>3){
              pbpos=input$bound[input$pos==prev1]
              bpos=input$bound[input$pos==npos]
              output$pos[j]=bpos
              output$x[j]=output$x[output$pos==pbpos]+vect$x[1]
              output$y[j]=output$y[output$pos==pbpos]+vect$y[1]
              j<-j+1
            }
          }else{
            ends<-c(prev1,pbpos)
            j=mf+5
            
          }
          prev1=npos
          
        }
      }
    }
    else{
      j=mf+5
    }
    
  }
  l<-NULL
  l[[1]]<-output
  l[[2]]<-ends
  l
}

loopLength<-function(input,start){
  stems<-NULL
  s<-1
  first<-input[input$pos==start,]
  mf<-length(input$pos)
  mp<-max(input$pos)
  j<-0
  ntot=1
  output<-NULL
  output$pos[1]=first$pos
  i<-2
  nextNT<-NULL
  nextNT=input[input$pos==first$pos,]
  bound=1
  while(j< mf){
    if(nextNT$bound!=0 & bound==0){
      npos=nextNT$bound
      bound=1
      
      stems[[s]]=c(nextNT$pos,npos)
      s<-s+1
    }else{
      bound=0
      npos=input$after[input$pos==nextNT$pos]
    }
    nextNT=input[input$pos==npos,]
    if(nextNT$pos==first$bound){
      j=mf
      output$pos[i]=nextNT$pos
    }else{
      if(nextNT$pos==mp){
        ntot=ntot+1
        j=mf
        output$pos[i]=nextNT$pos
      }else{
        ntot=ntot+1
        output$pos[i]=nextNT$pos
        j<-j+1
        i<-i+1
      }
      
    }
    
  }
  
  out<-NULL
  out[[1]]<-output
  out[[2]]<-stems
  out
}

rotateS<-function(x2,y2,x0,y0,ang){
  pt<-NULL
  x1<-x2-x0
  y1<-y2-y0
  pt[1]<-x1*cos(ang)-y1*sin(ang)+x0
  pt[2]<-y1*cos(ang)+x1*sin(ang)+y0
  pt
}

circleCoord<-function(n){
  ang=2*pi/(n)
  output<-NULL
  output$y[1]=0
  output$x[1]=1
  for(i in 1:(n-1)){
    v<-rotateS(1,0,0,0,(i*ang))
    j<-i+1
    output$x[j]<-v[1]+output$x[(j-1)]
    output$y[j]<-v[2]+output$y[(j-1)]
  }
  ### Rotate ##
  output
}

genCords<-function(loop,p1,p2,input,vn){
  
  dat<-loop[[1]]
  n=length(dat$pos)
  
  c<-circleCoord(n)
  
  pt1<-NULL
  pt2<-NULL
  pt1$x<-input$x[input$pos==p1]
  pt1$y<-input$y[input$pos==p1]
  pt2$x<-input$x[input$pos==p2]
  pt2$y<-input$y[input$pos==p2]
  
  pt2<-translate(pt2$x,pt2$y,pt1$x,pt1$y)
  ang<-atan2(pt2$y,pt2$x)
 
  if(vn==1){
    
    c$y=-1*c$y
  }
  c<-rotateV(c$x,c$y,0,0,ang)
  c<-translate(c$x,c$y,(-1*pt1$x),(-1*pt1$y))	
  return(c)
}

translate<-function(x1,y1,x2,y2){
  out<-NULL
  out$x<-x1-x2
  out$y<-y1-y2
  out
}

rotateV<-function(x2,y2,x0,y0,ang){
  ## Rotates and translates a vector ###
  pt<-NULL
  x1<-x2-x0
  y1<-y2-y0
  pt$x<-x1*cos(ang)-y1*sin(ang)+x0
  pt$y<-y1*cos(ang)+x1*sin(ang)+y0
  pt
}

