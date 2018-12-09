#########################################################################################
###################### Codispersion between X and Y #####################################

rho <- function(x, y, uvec, max.dist, angle)
{
  z <- as.geodata(cbind(x$coords, x$data + y$data))
  nz <- variog(z, uvec = uvec, max.dist = max.dist, direction = angle, messages = FALSE)
  dx <- variog(x, uvec = uvec, max.dist = max.dist, direction = angle, messages = FALSE)
  dy <- variog(y, uvec = uvec, max.dist = max.dist, direction = angle, messages = FALSE)
  rho <- .5 * (nz$v - dx$v - dy$v) / sqrt(dx$v * dy$v)
}

#########################################################################################
############### Codispersion Map for no Regular Lattice #################################

codisp.map <-
  function(x, y, coords, nclass=13, ncell=40, plot.it=TRUE, plotname="codismap.jpg", border=TRUE, dmax0=0.5,
           color=colorRampPalette(c("#0000FF","#FFFFFF","#FF6666"))(128), semicirc=TRUE, nangles=200, 
           zlim=c(-1,1), cex.lab=1.7, cex.main=2.0, cex.axis=1.5, lwd=8, width=512, height=480,
           xlab = "Spatial Lag X", ylab = "Spatial Lag Y", main="Codispersion Map")
  {
    require("akima")
    require("fields")
    require("geoR")
    require("jpeg")

    x <- as.geodata(cbind(coords, x))
    y <- as.geodata(cbind(coords, y))
    if(dmax0>1|dmax0<=0){ dmax <- 0.5 * max(dist(coords))}
    else{ dmax <- dmax0 * max(dist(coords)) }
    angles <- seq(from = 0, to = pi, l = nangles)
    uvec <- seq(from = 0, to = dmax, length = nclass + 1)[-1]
    
    xcirc <- 0
    ycirc <- 0
    
    for (i in seq_len(nclass))
    {
      xcirc[(nangles*(i-1)+1):(nangles*i)] <- seq(-uvec[i], uvec[i], length = nangles)
      ycirc[(nangles*(i-1)+1):(nangles*i)] <- sqrt(uvec[i]^2 - xcirc[(nangles*(i-1)+1):(nangles*i)]^2)
    }
    z <- matrix(0, nrow = nangles, ncol = nclass)
    for (i in seq_len(nangles)){z[i,] <- rho(x, y, uvec = uvec, max.dist = dmax, angle = angles[i])}
    z <- as.vector(z)
    
    if(semicirc)
    {
      xl <- seq(min(xcirc), max(xcirc), length=ncell)
      yl <- seq(min(ycirc), max(ycirc), length=ncell)
    
      if (plot.it) 
      {
        if(border)
        {
         par(pty = "s")
         image.plot(interp(xcirc, ycirc, z, xo = xl,yo = yl), col = color, cex.main=cex.main, 
                   cex.lab=cex.lab, zlim=zlim, cex.axis=cex.axis, ylim=c(-5,max(ycirc)+10),
                   xlab=xlab, ylab=ylab, main=main)
         xx=seq(min(xcirc)-1,max(xcirc)+1,length=100000)
         yy=sqrt((max(xcirc)+1)^2-xx^2)
         points(xx,yy,type="l",lwd=lwd)
         points(range(xx),c(0,0),type="l",lwd=lwd)
        }else{
          par(pty = "s")
          image.plot(interp(xcirc, ycirc, z, xo = xl,yo = yl), col = color, cex.main=cex.main, 
                     cex.lab=cex.lab, zlim=zlim, cex.axis=cex.axis, ylim=c(-5,max(ycirc)+10),
                     xlab=xlab, ylab=ylab, main=main)
        }
      }
      if(border)
      {
      jpeg(plotname, width=width, height=height)
        image.plot(interp(xcirc, ycirc, z, xo = xl,yo = yl), col = color, cex.main=cex.main, 
               cex.lab=cex.lab, cex.axis=cex.axis, zlim=zlim, ylim=c(-5,max(ycirc)+10),
                xlab=xlab, ylab=ylab, main=main)
        xx=seq(min(xcirc)-1,max(xcirc)+1,length=100000)
        yy=sqrt((max(xcirc)+1)^2-xx^2)
        points(xx,yy,type="l",lwd=lwd)
        points(range(xx),c(0,0),type="l",lwd=lwd)
      dev.off()
      }else{
        jpeg(plotname, width=width, height=height)
          image.plot(interp(xcirc, ycirc, z, xo = xl,yo = yl), col = color, cex.main=cex.main, 
                   cex.lab=cex.lab, cex.axis=cex.axis, zlim=zlim, ylim=c(-5,max(ycirc)+10),
                   xlab=xlab, ylab=ylab, main=main)
      }
    }else{
      xcirc2=c(xcirc,-xcirc)
      ycirc2=ycirc
      z2=z
      for(i in 1:length(z))
      {
        if(ycirc[i]==0){ycirc2[length(z)+i]<-NA}else{ycirc2[length(z)+i]<-(-ycirc[i])}
        z2[length(z)+i]<-z[i]
      }
      
      aux=na.omit(cbind(xcirc2,ycirc2,z2))
      xcirc=aux[,1]
      ycirc=aux[,2]
      z=aux[,3]
      
      xl <- seq(min(xcirc), max(xcirc), length=ncell)
      yl <- seq(min(ycirc), max(ycirc), length=ncell)
      
      if (plot.it) 
      {
        if(border)
        {
          par(pty = "s")
          image.plot(interp(xcirc, ycirc, z, xo = xl,yo = yl), col = color, cex.main=cex.main, 
                     cex.lab=cex.lab, zlim=zlim, cex.axis=cex.axis,
                     xlab=xlab, ylab=ylab, main=main)
          xx=seq(min(xcirc)-1,max(xcirc)+1,length=100000)
          yy=sqrt((max(xcirc)+1)^2-xx^2)
          points(xx,yy,type="l",lwd=lwd)
          points(xx,-yy,type="l",lwd=lwd)
        }else{
          par(pty = "s")
          image.plot(interp(xcirc, ycirc, z, xo = xl,yo = yl), col = color, cex.main=cex.main, 
                     cex.lab=cex.lab, zlim=zlim, cex.axis=cex.axis,
                     xlab=xlab, ylab=ylab, main=main)
        }
      }
      if(border)
      {
        jpeg(plotname, width = width, height = height)
        image.plot(interp(xcirc, ycirc, z, xo = xl,yo = yl), col = color, cex.main=cex.main, 
                 cex.lab=cex.lab, cex.axis=cex.axis, zlim=zlim,
                 xlab=xlab, ylab=ylab, main=main)
        xx=seq(min(xcirc)-1,max(xcirc)+1,length=100000)
        yy=sqrt((max(xcirc)+1)^2-xx^2)
        points(xx,yy,type="l",lwd=lwd)
        points(xx,-yy,type="l",lwd=lwd)
        dev.off()
      }else{
        jpeg(plotname, width = width, height = height)
        image.plot(interp(xcirc, ycirc, z, xo = xl,yo = yl), col = color, cex.main=cex.main, 
                   cex.lab=cex.lab, cex.axis=cex.axis, zlim=zlim,
                   xlab=xlab, ylab=ylab, main=main)
        dev.off()
      }
    }
    invisible(list(xcirc = xcirc, ycirc = ycirc, z = z))
  }
