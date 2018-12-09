###########################################################################################
############### Directions on Regular Grid ################################################

direction <- function(x, y, r=1, l=1)
{
  nc=dim(x)[2]
  nr=dim(x)[1]
  tt=min(nc,nr)
  h1=seq(0,ceiling(tt/l), by=r)
  h2=h1
  h1=c(-h1[length(h1):2],h1)
  H=expand.grid(h1,h2)
  H<-H[H[,1]^2+H[,2]^2<=(tt/l)^2,]
  print(dim(H)[1])
  invisible(list(xcirc = H[,1], ycirc = H[,2]))
}

###########################################################################################
############ Codispersion Map for Regular Lattice #########################################

hCodismap <- function(x, y, r=1, l=1, plot.it=TRUE, plotname="hcodismap.jpg", border=TRUE, ncell=40,
                      color=colorRampPalette(c("#0000FF","#FFFFFF","#FF6666"))(128), semicirc=TRUE, 
                      zlim=c(-1,1), cex.lab=1.7, cex.main=2.0, cex.axis=1.5, lwd=8, width=512, height=480,
                      xlab = "Spatial Lag X", ylab = "Spatial Lag Y", main="Codispersion Map")
{
  require("jpeg")
  require("fields")
  
  nc=dim(x)[2]
  nr=dim(x)[1]
  tt=min(nc,nr)
  hh1=seq(0,ceiling(tt/l), by=r)
  hh2=hh1
  hh1=c(-hh1[length(hh1):2],hh1)
  H=expand.grid(hh1,hh2)
  H<-H[H[,1]^2+H[,2]^2<=(tt/l)^2,]
  xcirc=H[,1]
  ycirc=H[,2]
  H3<-H
  H3[,1]<-(H[,1]+max(abs(H[,1]))+r)/r
  H3[,2]<-(H[,2]+r)/r

  aux <-  matrix(NA,ncol=max(H3[,2]),nrow=max(H3[,1]))
  
  for(k in 1:dim(H)[1])
  {
    aux[H3[k,1],H3[k,2]]=rhoh.c(x,y,hx=H[k,1],hy=H[k,2])
  }

  if(semicirc)
  {
    if(border)
    {
      jpeg(plotname, width=width, height=height)
      par(pty = "s")
      image.plot(x=hh1, y=hh2, z=aux, col = color, cex.lab=cex.lab, cex.axis=cex.axis,
                 main=main, xlab=xlab, ylab=ylab, cex.main=cex.main, zlim=zlim,
                 xlim=c(min(hh1)-r-3,max(hh1)+r+3), ylim=c(-3,max(hh2)+r+3) ) 
      xx=seq(min(hh1)-r,max(hh1)+r,length=100000)
      yy=sqrt((max(hh1)+r)^2-xx^2)
      points(c(min(hh1)-r,max(hh1)+r),c(-0.05,-0.05),type="l",lwd=lwd)
      points(xx,yy,type="l",lwd=lwd)
      dev.off()
    }else{
      jpeg(plotname, width=width, height=height)
      par(pty = "s")
      image.plot(x=hh1, y=hh2, z=aux, col = color, cex.lab=cex.lab, cex.axis=cex.axis,
                 main=main, xlab=xlab, ylab=ylab, cex.main=cex.main, zlim=zlim,
                 xlim=c(min(hh1)-r-3,max(hh1)+r+3), ylim=c(-3,max(hh2)+r+3) ) 
      dev.off()
    }
    if(plot.it)
    {
      if(border)
      {
        par(pty = "s")
        image.plot(x=hh1, y=hh2, z=aux, col = color, cex.lab=cex.lab, cex.axis=cex.axis,
                   main=main, xlab=xlab, ylab=ylab, cex.main=cex.main, zlim=zlim,
                   xlim=c(min(hh1)-r-3,max(hh1)+r+3), ylim=c(-3,max(hh2)+r+3) ) 
        xx=seq(min(hh1)-r,max(hh1)+r,length=100000)
        yy=sqrt((max(hh1)+r)^2-xx^2)
        points(c(min(hh1)-r,max(hh1)+r),c(-0.05,-0.05),type="l",lwd=lwd)
        points(xx,yy,type="l",lwd=lwd)
      }else{
        par(pty = "s")
        image.plot(x=hh1, y=hh2, z=aux, col = color, cex.lab=cex.lab, cex.axis=cex.axis,
                   main=main, xlab=xlab, ylab=ylab, cex.main=cex.main, zlim=zlim,
                   xlim=c(min(hh1)-r-3,max(hh1)+r+3), ylim=c(-3,max(hh2)+r+3) ) 
      }
    }
    z=aux
  }else{
    aux2 <- matrix(NA,ncol=2*max(H3[,2])-1,nrow=max(H3[,1]))
    aux2[,max(H3[,2]):(2*max(H3[,2])-1)] <- aux
    aux2[,1:(max(H3[,2])-1)] <- aux[max(H3[,1]):1,max(H3[,2]):2]

    if(border)
    {
      jpeg(plotname, width=width, height=height)
      par(pty = "s")
      image.plot(x=hh1, y=hh1, z=aux2, col = color, cex.main=cex.main, 
                 cex.lab=cex.lab, zlim=zlim, cex.axis=cex.axis, xlab=xlab, ylab=ylab, main=main,
                 xlim=c(min(hh1)-r-3,max(hh1)+r+3), ylim=c(min(hh1)-r-3,max(hh1)+r+3) ) 
      xx=seq(min(xcirc)-r,max(xcirc)+r,length=100000)
      yy=sqrt((max(xcirc)+r)^2-xx^2)
      points(xx,yy,type="l",lwd=lwd)
      points(xx,-yy,type="l",lwd=lwd)
      dev.off()
    }else{
      jpeg(plotname, width=width, height=height)
      image.plot(x=hh1, y=hh1, z=aux2, col = color, cex.main=cex.main, 
                 cex.lab=cex.lab, zlim=zlim, cex.axis=cex.axis, xlab=xlab, ylab=ylab, main=main,
                 xlim=c(min(hh1)-r-3,max(hh1)+r+3), ylim=c(min(hh1)-r-3,max(hh1)+r+3) ) 
      dev.off()
    }
    if(plot.it)
    {
      if(border)
      {
        par(pty = "s")
        image.plot(x=hh1, y=hh1, z=aux2, col = color, cex.main=cex.main, 
                   cex.lab=cex.lab, zlim=zlim, cex.axis=cex.axis, xlab=xlab, ylab=ylab, main=main,
                   xlim=c(min(hh1)-r-3,max(hh1)+r+3), ylim=c(min(hh1)-r-3,max(hh1)+r+3) ) 
        xx=seq(min(xcirc)-r,max(xcirc)+r,length=100000)
        yy=sqrt((max(xcirc)+r)^2-xx^2)
        points(xx,yy,type="l",lwd=lwd)
        points(xx,-yy,type="l",lwd=lwd)
      }else{
        image.plot(x=hh1, y=hh1, z=aux2, col = color, cex.main=cex.main, 
                   cex.lab=cex.lab, zlim=zlim, cex.axis=cex.axis, xlab=xlab, ylab=ylab, main=main,
                   xlim=c(min(hh1)-r-3,max(hh1)+r+3), ylim=c(min(hh1)-r-3,max(hh1)+r+3) ) 
      }
    }
    z=aux2
    H=expand.grid(hh1,hh1)
    H<-H[H[,1]^2+H[,2]^2<=(tt/l)^2,]
    xcirc=H[,1]
    ycirc=H[,2]
  }
  invisible(list(xcirc = xcirc, ycirc = ycirc, z = z))
}


####################################################################################
##################### Call functions ############################################### 
rhoh.c <- function(xM,yM,hx,hy)
{
  ncol=dim(xM)[2]
  nrow=dim(xM)[1]
  rhoh3 <- .C("hcod_coef", as.double(xM), as.double(yM), as.integer(nrow), as.integer(ncol),
              as.integer(hx), as.integer(hy), rhoh = double(1) ) 
  rhoh3[["rhoh"]]
}

