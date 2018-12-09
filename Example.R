#######################################################
################## Sources ############################

source("Noises.R")
source("Input.R")
source("Codismap.R")
source("hCodismap.R")
dyn.load("hCodismap.so")

#######################################################
##### Example Regular Grid ############################

sim1=simmatern(256, 0.5,0.5,1,1,0.5,0.5,0.4)

X1 <- sim1[[1]]
Y1 <- sim1[[2]]

SP_Y1 <- saltnpepper(Y1, epsilon = 0.25)
NAcont_Y1 <- NAcont(Y1, epsilon=0.001, k=5)
NAblock_Y1<- NAblock(Y1, k=20, m=1)
Input_Y1 <- Input(y=NAblock_Y1, k=20)

par(mfrow=c(2,3), pty="s",mar=c(2,2,1,1))
image(X1,main="Image X1", col=gray((0:32)/32))
image(Y1,main="Image Y1", col=gray((0:32)/32))
image(SP_Y1,main="Salt and Pepper", col=gray((0:32)/32))
image(NAcont_Y1,main="Missing Block", col=gray((0:32)/32))
image(NAblock_Y1,main="Gaps", col=gray((0:32)/32))
image(Input_Y1,main="Input", col=gray((0:32)/32))

hCodismap(x=X1, y=Y1, r=3, l=2, lwd=4, plotname="example1.jpg")
hCodismap(x=X1, y=SP_Y1, r=3, l=2, lwd=4, semicirc=FALSE, plotname="example2.jpg", zlim=c(-0.5,1))
hCodismap(x=X1, y=NAcont_Y1, r=3, l=2, lwd=4, semicirc=FALSE, plotname="example3.jpg", zlim=c(-0.5,1))
hCodismap(x=X1, y=NAblock_Y1, r=2, l=2, lwd=4, plotname="example4.jpg", col=topo.colors(256))
hCodismap(x=X1, y=Input_Y1, r=3, l=2, lwd=4, plotname="example5.jpg", col=gray((0:256)/256), zlim=c(-0.5,1))

color=colorRampPalette(c("#0000FF","#FFFFFF","#FF6666"))(128)
X11()
par(mfrow=c(2,3), pty="s",mar=c(2,2,1,1))
a<-hCodismap(x=X1, y=Y1, r=1, l=1, lwd=4, plotname="example1.jpg", border=FALSE)
plot(a$xcirc,a$ycirc)
image(x=unique(a$xcirc), y=unique(a$ycirc), z=a$z, col = color, zlim=c(-0.5,1))

b<-hCodismap(x=X1, y=SP_Y1, r=1, l=1, lwd=4, semicirc=FALSE, plotname="example2.jpg", zlim=c(-0.5,1), border=FALSE)
plot(b$xcirc,b$ycirc)
x1=unique(b$xcirc); length(x1)
y1=unique(b$ycirc); length(y1)
image(x=x1[order(x1)], y=y1[order(y1)], z=b$z, col = color, zlim=c(-0.5,1))

##########################################################
##### Example No-Regular Grid ############################

coord=matrix(round(runif(500,0,200),2),nc=2)
plot(coord)
data1=rnorm(250)
data2=rnorm(250)

aa<-codisp.map(x=data1, y=data2, coords=coord,lwd=4, dmax0=0.3)
bb<-codisp.map(x=data1, y=data2, coords=coord,lwd=4, semicirc=FALSE, plotname="example_codismap.jpg")



