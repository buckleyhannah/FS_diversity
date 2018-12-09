##########################################################################
###################### Salt and Pepper ###################################

saltnpepper <- function(x, epsilon = 0.10)
{
  xrow <- nrow(x)
  xcol <- ncol(x)
  eps <- rbinom(xrow * xcol, size = 1, prob = epsilon)
  eps <- matrix(eps, nrow = xrow, ncol = xcol)
  
  xbar <- mean(x)
  xsd <- sd(x)
  
  noise <- rnorm(xrow * xcol, mean = xbar, sd = xsd)
  y <- (1 - eps) * x + eps * noise
  y <- ifelse(y > 1, 1, y)
  y <- ifelse(y < 0, 0, y)
  y
}


##########################################################################
####### Missing Blocks Observations at Random Locations ##################

NAcont <- function(x, epsilon=0.01, k=1)
{
  xrow <- nrow(x)
  xcol <- ncol(x)
  bin <- rbinom(xrow * xcol, size = 1, prob = epsilon)
  eps <- matrix(bin, nrow = xrow, ncol = xcol)
  eps2 <- eps
  for(j in 1:(xcol))
  {
    for(i in 1:(xrow))
    {
      if(eps[i,j]==1) eps2[max(1,(i-k)):min((i+k),xrow),max(1,(j-k)):min(xcol,(j+k))]<-1
    }
  }
  y <- matrix(xrow * xcol, nrow = xrow, ncol = xcol)
  y <- (1 - eps2) * x + eps2
  y
  
}

##########################################################################
####### Gaps Resulting from Clusters of Missing Observations #############

NAblock<- function(x, k=10, m=1)
{
  xrow <- nrow(x)
  xcol <- ncol(x)
  eps <- x
  if(length(m)==1){m=c(m,m)}
  rx=index(xrow,k,m[1])
  ry=index(xcol,k,m[2])
  eps[rx,ry]<-1
  eps
}

## Auxiliary functions for Gaps function.

Par<-function(x)
{
  ceiling(x/2)==x/2
}

index<-function(a,k,m)
{
  aux=c()
  if(Par(k)==FALSE)
  {
    for(i in 1:m)
    {
      mx=round(i*a/(m+1))
      aux=c(aux, (mx-(k+1)/2+1):(mx+(k+1)/2-1))
    }
  }else
  {
    for(i in 1:m)
    {
      mx=round(i*a/(m+1))
      aux=c(aux, (mx-k/2+1):(mx+k/2))
    }
  }
  aux
}

##########################################################################
#################### Bivariate Matern Simulation #########################

simmatern<-function(N, nu1=0.5, nu2=0.5, varx=1, vary=1, meanx=0.5, meany=0.5, rho12=0.1, plot=TRUE)
{
  require("RandomFields")
  
  xpos <- ypos <- seq(1, N, 1) 												  
  model <- RMbiwm(nudiag = c(nu1,nu2), c = c(varx, rho12, vary))			
  data <- RFsimulate(model, xpos, ypos)										
  X <- matrix(data$variable1, ncol = N) + meanx								
  Y <- matrix(data$variable2, ncol = N) + meany	
  
  X <- (X-min(X))/(max(X)-min(X))
  Y <- (Y-min(Y))/(max(Y)-min(Y))
  if (plot)
    {
     par(mfrow=c(1,2), pty="s",mar=c(2,2,1,1))
     image(X,main="Image X", col=gray((0:32)/32))
     image(Y,main="Image Y", col=gray((0:32)/32))
     r=cor(c(X),c(Y))
    }
  return(list(X,Y))
}


