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

hCodismap <- function(x, y, r=1, l=1)  {
  nc=dim(x)[2] # number of columns
  nr=dim(x)[1] # numbner of rows
  tt=min(nc,nr) # find the minimum plot dimension
  hh1=seq(0,ceiling(tt/l), by=r) # divide by the value to set quarter of the plot size 
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
  z=aux
  return(list(xcirc = xcirc, ycirc = ycirc, z = z))
  }

####################################################################################
##################### Call functions ############################################### 

rhoh.c <- function(xM,yM,hx,hy)

  {
  ncol=dim(xM)[2]
  nrow=dim(xM)[1]
  rhoh3 <- .C("hcod_coef", as.double(xM), as.double(yM), as.integer(nrow), as.integer(ncol), as.integer(hx), as.integer(hy), rhoh = double(1) ) 
  rhoh3[["rhoh"]]
}

