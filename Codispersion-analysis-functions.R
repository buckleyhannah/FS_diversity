#####################################
### Codispersion analysis functions
#####################################

#### Modified codispersion function (modified from Cuevas et al. 2013)
#### See 'Box 1' in Buckley et al. (2016) for a detailed explanation.
#### Buckley, HL., Case, BS., Ellison, AM. 2016. Using codispersion analysis to characterize spatial patterns in species co-occurrences. Ecology 97(1): 32–39.

Codisp.Kern <- function(X, Y, h, k, gamma = 1)
{
  Kernel<-function(u,gamma)
  {
    v=0
    v=ifelse(abs(u)<=1,(1/beta(0.5,gamma+1))*(1-u^2)^gamma,0)
  }
  ifelse(X$coords==Y$coords,1,
         {
           break
           print("The coordinates of X and Y are different")
         })
  
  n=length(X$data)
  mX <- matrix(X$data,nrow=n,ncol=n,byrow=FALSE)
  mY <- matrix(Y$data,nrow=n,ncol=n,byrow=FALSE)
  MatriXX <- (mX - t(mX))^2
  MatriYY <- (mY - t(mY))^2
  MatriXY <- (mX - t(mX))*(mY - t(mY))
  mX <- matrix(X$coords[,1],nrow=n,ncol=n,byrow=FALSE)
  DesignX <- mX - t(mX)
  mY <- matrix(X$coords[,2],nrow=n,ncol=n,byrow=FALSE)
  DesignY <- mY - t(mY)
  
  KERNMATRIXX=Kernel((h[1]-DesignX)/k[1],gamma)*Kernel((h[2]-DesignY)/k[1],gamma)
  
  if(k[1]==k[2]&k[1]==k[3]){
    KERNMATRIYY=KERNMATRIXX
    KERNMATRIXY=KERNMATRIXX } else{
      KERNMATRIYY=Kernel((h[1]-DesignX)/k[2],gamma)*Kernel((h[2]-DesignY)/k[2],gamma)
      KERNMATRIXY=Kernel((h[1]-DesignX)/k[3],gamma)*Kernel((h[2]-DesignY)/k[3],gamma) 
    }
  
  Numerador=sum(KERNMATRIXY*MatriXY)/(2*sum(KERNMATRIXY))
  Denominador1=sum(KERNMATRIYY*MatriYY)/(2*sum(KERNMATRIYY))
  Denominador2=sum(KERNMATRIXX*MatriXX)/(2*sum(KERNMATRIXX))
  v1=Denominador1
  v2=Denominador2
  v3=Numerador
  v4=Numerador/sqrt(Denominador1*Denominador2)
  print(c(v1,v2,v3,v4))
}  



#### Function to run codispersion window analysis (modified from Cuevas et al. 2013)

# geodata1 = first input data object (a geoR geodata object)
# geodata2 = second input object
# k = c(k1, k2, k3) = a vector of three bandwidth values for X, Y and XY
# max.window.size = the maximum lag distance
# lx = is the number of divisions in the lags in x (up to the max.window.size) that the kernal is applied over. Half of these divisions are in the 'left', or positive direction, and half are in the 'right', or negative x direction.
# ly = is the number of divisions in the lags in y (up to the max.window.size) that the kernal is applied over in the 'up' direction of the plot

codisp.fn <- function(geodata1, geodata2, k = k, max.window.size = max.window.size, lx = 20, ly = 10) {
  X = geodata1  # input data process 1
  Y = geodata2  # input data process 2
  k = c(k[1],k[2],k[3]) # Set the bandwith for the kernel
  
  h_range <- max.window.size     # set the spatial lags over which to calculate codisp
  h1 = seq(-h_range, h_range, l = lx)  # x-axis values for codispersion graph (lags)
  h2 = seq(min(k), h_range, l = ly)    # y-axis values for codispersion graph (lags)
  
  
  MCodisp = matrix(0, ncol = ly, nrow = lx) # loop through the lags
  for(i in 1:lx)     # 'left-right' lags
  {
    for(j in 1:ly)   # 'up' lags
    {
      MCodisp[i,j] = Codisp.Kern(X, Y, c(h1[i], h2[j]), k)[4]; # calculate codisp using the kernel function given above
    }
  }
  Codispersion <- as.numeric(MCodisp) # save codisp object as output
  X <- rep(h1, length(h2))            # write out values for x-axis
  Y <- rep(h2, each = length(h1))       # write out values for y-axis
  graphing.data <- data.frame(X, Y, Codispersion) # graphing object
  
  # output the dataframe
  return(graphing.data) 
  
}  # end function


#### Plotting variograms and cross variograms

cross.variog <- function(geoR.v1, geoR.v2, v1.name = "var1", v2.name = "var2") {
  
  max.x <- ceiling(max(geoR.v1$coords[,1])) # determine plot dimensions
  max.y <- ceiling(max(geoR.v1$coords[,2]))
  
  temp <- data.frame(geoR.v1$coords, var1 = scale(geoR.v1$data), var2 = scale(geoR.v2$data)) # create dataframe with coordinates and two data vectors
  names(temp)[1] <- "gx" # rename the variables for variogram calculations 
  names(temp)[2] <- "gy"
  names(temp)[3] <- "var1" 
  names(temp)[4] <- "var2"
  
  # compute the variogram and cross variogram
  g <- gstat(id = v1.name, formula = var1 ~ 1, locations = ~ gx + gy, data = temp)    
  g <- gstat(g, id = v2.name, formula = var2 ~ 1, locations = ~ gx + gy, data = temp)
  v <- variogram(g, cutoff = (min(max.x, max.y) * 0.67), cross = TRUE)
  
  # plot the figure
  print(ggplot(v, aes(x = dist, y = gamma, group = id, colour = id)) + 
    geom_line(lwd = 2) + 
    labs(x = "Distance (m)", y = "Semivariance") +
    theme_bw(base_size = 14) +
    theme(legend.title = element_blank())  
      ) # end print plot
  
} # end function


#### Plotting codispersion graphs

plot.codisp <- function(codisp.out, v1.name = "var1", v2.name = "var2", scaled = TRUE, binwidth = 0.1) {
  p <- ggplot(codisp.out, aes(x = X, y = Y, fill = Codispersion)) +
  coord_fixed(ratio = 1) + 
  geom_tile() + 
  stat_contour(aes(x = X, y = Y, z = Codispersion), binwidth = binwidth) +
  xlab('Spatial lag in X (m)') + ylab('Spatial lag in Y (m)') +
  theme_bw() +
  ggtitle(paste(v1.name, "vs.", v2.name))
  if(scaled == TRUE) { p <- p + scale_fill_gradientn(colours = c('#0000FF', '#FFFFFF', '#FF6666'), limits = c(-1, 1)) }
  if(scaled == FALSE) { p <- p + scale_fill_gradientn(colours = rainbow(20)) }
  print(p)
  return(p)
} # end function

plot.null.model.result <- function(null.model.result, v1.name = "var1", v2.name = "var2", null.model = null.model, graph = c("OE", "PV")) {
  
  # Observed minus expected CoDispersion value graph
  if(graph == "OE") {
    print(ggplot(null.model.result, aes(X, Y)) + 
      geom_tile(aes(fill = Difference)) + 
      scale_fill_gradientn(colours = c("#0000FF", "#FFFFFF", "#FF6666"), limits = c(-1, 1)) + 
      coord_fixed(ratio = 1) +
      stat_contour(aes(x = X, y = Y, z = Difference), binwidth = 0.1) +
      xlab('Spatial lag in X (m)') + ylab('Spatial lag in Y (m)') +
      theme_bw() + theme(legend.title = element_blank()) +
      ggtitle(paste(v1.name, "vs.", v2.name, ":", null.model)) )
      }
  
  # P-value category graph
  if(graph == "PV") {
  my.cols <- c("steelblue3", "firebrick3")  # select colours for graph
  if(levels(null.model.result$P.value.cat)[1] == "Sig.") {my.cols <- c("firebrick3") }
  print( ggplot(null.model.result, aes(X, Y)) +
    geom_tile(aes(fill = P.value.cat)) + 
    scale_fill_manual(values = my.cols) +
    scale_colour_discrete(name = "P.value.cat", limits = c(0,1)) +
    coord_fixed(ratio = 1) +
    xlab('Spatial lag in X (m)') + ylab('Spatial lag in Y (m)') +
    theme_bw() + theme(legend.title = element_blank()) +
    ggtitle(paste(v1.name, "vs.", v2.name, ":", null.model)) )
    }
  } # end function


#### Simulating covarying raster patterns

#   grid.points = 20     # scale: grain of grid
#   xdim = 500           # Dimensions of the plot area, e.g. 500 x 500m 
#   ydim = 500
#   sp1.pattern = "CSR"  # The desired pattern for species 1: CSR, decreasing.x, increasing.x, decreasing.xy, bivariate.normal
#   sp2.pattern = "CSR"  # The desired pattern for species 2: CSR, decreasing.x, increasing.x, decreasing.xy, bivariate.normal
#   Print = TRUE  # Whether you want plots of the spatial patterns or not

simulate.rasters <- function(grid.points = grid.points, sp1.pattern = c("CSR", "decreasing.x", "increasing.x", "decreasing.y", "decreasing.xy", "increasing.xy", "bivariate.normal"), sp2.pattern = c("CSR", "decreasing.x", "increasing.x", "increasing.y", "decreasing.xy", "increasing.xy", "bivariate.normal", "inc.x.dec.y", "dec.x.inc.y"), xdim = xdim, ydim = ydim, print = c("TRUE", "FALSE")){
  
  # 1. Create an empty list to add output geodata objects
  copp.sim <- vector("list", 2)
  
  # 2. Set up underlying grid coordinates
  X <- seq(from = 0, to = xdim - grid.points, by = grid.points)
  Y <- seq(from = 0, to = ydim - grid.points, by = grid.points)
  gridxy <- expand.grid(x = X, y = Y)
  
  # 3a. Create a set of quadrat abundance values for sp1 based on the 'sp1.pattern' argument
  
  if(sp1.pattern == "CSR") {Z <- rnorm(length(gridxy$x), mean = 30, sd = 10) }
  
  if(sp1.pattern == "decreasing.x") {Z <- 1 + (rev(2 * gridxy$x + 5)) / 10 }
  
  if(sp1.pattern == "increasing.x") {Z <- 1 + (2*gridxy$x + 5) / 10 }
  
  if(sp1.pattern == "decreasing.y") {Z <- 1 + (rev(2 * gridxy$y + 5)) / 10 }
  
  if(sp1.pattern == "decreasing.xy") { Z <- 1 + rev(((gridxy$x + 1)^2 + (gridxy$y + 1)^2) / 3000) }
  
  if(sp1.pattern == "increasing.xy") { Z <- 1 + ((gridxy$x + 2)^2 + (gridxy$y + 1)^2) / 3000 }
  
  if(sp1.pattern == "bivariate.normal") { Z <- 300 * bivariate(((gridxy$x - min(gridxy$x)) / (max(gridxy$x) - min(gridxy$x)) * 4) - 2, ((gridxy$y - min(gridxy$y)) / (max(gridxy$y) - min(gridxy$y)) * 4) - 2) } #  bivariate normal
  
  if(sp1.pattern == "inc.x.dec.y") {  Z <- 1 + ((gridxy$x + 2)^2 + (rev(gridxy$y + 1))^2) / 3000 }
  
  if(sp1.pattern == "dec.x.inc.y") { Z <- 1 + ((rev(gridxy$x + 2))^2 + (gridxy$y + 1)^2) / 3000 }
  
  # 3b. Add data to the output list as a geodata object
  copp.sp1.df <- data.frame(x = gridxy$x, y = gridxy$y, Z = jitter(Z, mean(Z) / 5))
  copp.sim[[1]] <- as.geodata(copp.sp1.df, coords.col = 1:2, data.col = 3)
  
  # 4a. Create a set of quadrat abundance values for sp2 based on the 'sp1.pattern' argument
  
  if(sp2.pattern == "CSR"){Z <- rnorm(length(gridxy$x),mean=30,sd=10) }
  
  if(sp2.pattern == "decreasing.x") { Z <- 1 + (rev(2 * gridxy$x + 5)) / 10 }
  
  if(sp2.pattern == "decreasing.x") { Z <- 1 + (rev(2 * gridxy$x + 5)) / 10 }
  
  if(sp2.pattern == "increasing.x") { Z <- 1 + (2 * gridxy$x + 5) / 10 }
  
  if(sp2.pattern == "increasing.y") { Z <- 1 + (2 * gridxy$y + 5) / 10 }
  
  if(sp2.pattern == "decreasing.xy") { Z <- 1 + rev(((gridxy$x + 1)^2 + (gridxy$y + 1)^2) / 3000) }
  
  if(sp2.pattern == "increasing.xy") { Z <- 1 + ((gridxy$x + 2)^2 + (gridxy$y + 1)^2) / 3000 }
  
  if(sp2.pattern == "bivariate.normal") { 
    Z <- 300 * bivariate(((gridxy$x - min(gridxy$x)) / (max(gridxy$x) - min(gridxy$x)) * 4) - 2, ((gridxy$y - min(gridxy$y)) / (max(gridxy$y) - min(gridxy$y)) * 4) - 2)  } 
  
  if(sp2.pattern == "inc.x.dec.y") { Z <- 1 + ((gridxy$x + 2)^2 + (rev(gridxy$y + 1))^2) / 3000 }
  
  if(sp2.pattern == "dec.x.inc.y") {  
    Z <- 1 + ((rev(gridxy$x + 2))^2 + (gridxy$y + 1)^2) / 3000 }
  
  # 4b. Add data to the output list as a geodata object
  copp.sp2.df <- data.frame(x = gridxy$x, y = gridxy$y, Z = jitter(Z, mean(Z)/10))
  copp.sim[[2]] <- as.geodata(copp.sp2.df, coords.col = 1:2, data.col = 3) 
  
  # 5. Print abundance maps if desired
  if(print == "TRUE") {
    gdat <- data.frame(x = rep(copp.sim[[1]]$coords[,1], 2), y = rep(copp.sim[[1]]$coords[,2], 2), sp = rep(c("Species 1", "Species 2"), each = length(copp.sim[[1]]$data)), ab = c(copp.sim[[1]]$data, sp2 = copp.sim[[2]]$data)) 
    
    print( ggplot(gdat, aes(x = x, y = y, fill = sp)) + geom_point(aes(size = ab), alpha = 0.6, shape = 21) + scale_fill_manual(values = c('darkseagreen3', 'darkslategray4')) + facet_grid(. ~ sp) + labs(x = "X", y = "Y") + theme_bw(base_size = 14) + theme(legend.position = "none") )
    } # end print loop    
  
  # 6. Output the list of geodata objects
  return(copp.sim)  
} # end function


#### Function to generate a geodata object (used by packages geoR and the codispersion function) from a ppp object.

# ppp.dat = input ppp object
# xmin, xmax, ymin, ymax = plot dimensions
# method = the measure that is used to generate the 'data' value for the geodata object

ppp.to.geoR.fn <- function(ppp.dat, xmin, xmax, ymin, ymax, quad.size, method = c("abundance","mean.mark","mean.ba","total.ba","sum")) {  # function to generate geoR objects with abundance and basal area in 20x20m quadrats. Note that DBH must be measured in cm. Input data = ppp object.
  x <- ppp.dat$x # extract x coordinate
  y <- ppp.dat$y # extract y coordinate
  z <- ppp.dat$marks # extract DBH values
  ba <- (pi * (z)^2) / 40000 # calculate basal area in m^2 
  xt <- cut(x, seq(xmin, xmax, quad.size)) # cut x coordinates using quad.size spacing
  yt <- cut(y, seq(ymin, ymax, quad.size)) # cut y coordinates using quad.size spacing
  coords <- dimnames(table(yt, xt)) # extract quadrat coordinate lists
  qx <- rep(seq(xmin, xmax - quad.size, length = length(coords$xt)), each = length(coords$yt)) # vector of x coordinates for the bottom left corner of the quadrat
  qy <- rep(seq(ymin, ymax - quad.size, length = length(coords$yt)), length(coords$xt)) # vector of y coordinates for the bottom left corner of the quadrat
  if(method == "abundance"){
    out.grid <- table(yt, xt) # count the trees in each quadrat
    out.grid[is.na(out.grid) == T] <- 0 # replace NAs in table with zeros for empty quadrats
  }
  if(method == "mean.mark"){
    out.grid <- tapply(z, list(yt, xt), mean) # calculate mean DBH in each quadrat
    out.grid[is.na(out.grid) == T] <- 0
  }    
  if(method == "mean.ba"){
    out.grid <- tapply(ba, list(yt, xt), mean) # calculate mean ba in each quadrat
    out.grid[is.na(out.grid) == T] <- 0 
  }    
  if(method == "total.ba"){
    out.grid <- tapply(ba, list(yt, xt), sum) # calculate total ba in each quadrat    
    out.grid[is.na(out.grid) == T] <- 0 
  }  
  if(method == "sum"){
    out.grid <- tapply(z, list(yt, xt), sum) # calculate sum of the marks in each quadrat
    out.grid[is.na(out.grid) == T] <- 0
  }    
  
  out.df <- data.frame(qx, qy, as.vector(out.grid))
  out.geo <- as.geodata(out.df, coords.col = 1:2, data.col = 3)
  return(out.geo)
} # end function

#### geoR object to matrix input object for hCodisp function
 # (test.df <- data.frame(x = rep(1:4, 3), y = rep(1:3, each = 4), z = 21:32))
 # (mat.obs.tab.sp1 <- matrix(test.df$z, ncol = length(unique(test.df$x)), nrow = length(unique(test.df$y)), byrow = TRUE))
 # (mat.obs.tab.sp1 <- mat.obs.tab.sp1[ nrow(mat.obs.tab.sp1):1, ])
 # plot(mat.geo <- as.geodata(test.df, coords.col = 1:2, data.col = 3))
 # mat.geo
 # geoR.to.matrix.fn(mat.geo)
 # rand <- sample(nrow(test.df))
 # new.test.df <- test.df[rand, ]
 
geoR.to.matrix.fn <- function(geoR.obj) { 
  temp.df <- data.frame(x = geoR.obj$coords[,1], y = geoR.obj$coords[,2], z = geoR.obj$data)
  temp.df <- dplyr::arrange(temp.df, x, y)
  mat <- matrix(temp.df$z, ncol = length(unique(temp.df$x)), nrow = length(unique(temp.df$y)), byrow = FALSE)
  #mat <- mat[nrow(mat):1, ]
  return(mat)
}


### Function to extract values from an environmental grid (raster) at point locations in a ppp object. 
# Inputs are the ppp and geoR.env, which is a geoR object holding the environmental data layer.
# Extract works with a buffer around each point (10cm in this case)
# geoR.env <- geo.elev

extract.env.fn <- function(ppp.dat, geoR.env, xmin, xmax, ymin, ymax, quad.size = quad.size) {
  
  ppp.df <- data.frame(x = ppp.dat$x, y = ppp.dat$y)   # create a dataframe from the ppp object
  
  X = geoR.env$coords[, 1] + quad.size / 2
  Y = geoR.env$coords[, 2] + quad.size / 2
  
  rast <- raster()  # create empty raster to add data to
  extent(rast) <- extent(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) 
  ncol(rast) <- xmax / quad.size 
  nrow(rast) <- ymax / quad.size 
  raster.env <- rasterize(cbind(X, Y), rast, geoR.env$data)
  
  env.dat <- extract(x = raster.env, y = ppp.df, df = TRUE) # extract environmental data at tree locations
  
  env.out <- data.frame(x = ppp.df$x, y = ppp.df$y, z = env.dat$layer)
  env.geo <- as.geodata(env.out, coords.col = 1:2, data.col = 3)
  return(env.geo)
}

#------------------------------------
#### Functions to run null modelling
#------------------------------------
  
#### Function to generate a list of ppp objects under three different null models
  # nrand = the number of randomisations you want
  # model = the chosen model under which to randomise the point pattern
  # marks = whether or not the point pattern contains marks, e.g., DBH values
  
ppp.null.model <- function(ppp.dat, nrand, model = c("RLM", "HomP", "Tor"), marks = TRUE)   {
  
  ppp.out <- vector("list", nrand) # create output list object  
  
  if(model == "RLM"){  # Random labelling model
    for(i in 1:nrand) { # start loop to generate randomisations
      ppp.out[[i]] <- rlabel(ppp.dat, labels = marks(ppp.dat), permute = TRUE) # randomise marks
    }  # end randomisations loop
  }  # end RLM loop  
  
  if(model == "HomP") {  # Homogeneous Poisson model
    for(i in 1:nrand) {  # start loop to generate randomisations
      ppp.HomP <- rpoint(ppp.dat$n, win = ppp.dat$win) # randomise the observed ppp
      ppp.HomP$marks <- sample(ppp.dat$marks, replace = F) # assign shuffled marks to new ppp
      ppp.out[[i]] <- ppp.HomP # add new marked ppp to output list
    } # end randomisations loop    
  } # end HomP loop  
  
  if(model == "Tor") {  # Toroidal shift null model
    for(i in 1:nrand){  # start loop to generate randomisations
      ppp.out[[i]] <- rshift(ppp.dat, edge = "torus", group = NULL) # toroidal shift randomisation 
    }  # end randomisations loop
  }  # end toroidal shift  
  
  return(ppp.out)
} # end function


# List to array function for Co_disp null model output objects
list2ary = function(input.list) {  #input a list of lists
  temp.ls <- vector("list", length(input.list))
  for(i in 1:length(input.list)) { temp.ls[[i]] <- input.list[[i]]  } # take the dataframes out of the list and put them in a new list
  rows.cols <- dim(temp.ls[[1]])
  sheets <- length(temp.ls)
  output.ary <- array(unlist(temp.ls), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(temp.ls[[1]])
  row.names(output.ary) <- row.names(temp.ls[[1]])
  return(output.ary)    # output as a 3-D array
}


# Function to return a data frame with the null model comparison results
  # null.input.ary = array of codispersion outputs from null model randomisations
# CoDisp_obs = observed codispersion result

codisp.compare <- function(CoDisp_obs, null.input.ary, round = FALSE) {
  out.df <- CoDisp_obs # observed Codispersion result df
  for(i in 1:length(null.input.ary[, 1, 1])) { # loop through each cell
    nsims <- length(null.input.ary[1, 1, ])
    obser <- out.df$Codispersion[i] # observed codispersion value
    expec <- null.input.ary[i, 3, ]
    prop.greater.than <- length(which(expec > obser)) / nsims
    prop.less.than <- length(which(expec < obser)) / nsims
    out.df$P.value[i] <- min(prop.greater.than, prop.less.than)
  } # end cell loop
  
  out.df$null_mean <- apply(null.input.ary[, 3, ], MARGIN = 1, mean)  # calculate mean codispersion value for each cell from the array of null model results
  out.df$Difference <- out.df$Codispersion-out.df$null_mean  # observed minus expected
  out.df$P.value.cat <- factor(ifelse(out.df$P.value < 0.025, "Sig.", "Non-sig.")) # significance at alpha = 0.05
  
  if(round == TRUE){  # for printing table of results
    out.df$X <- round(out.df$X, 1)
    out.df$Y <- round(out.df$Y, 1)
    out.df$Codispersion <- round(out.df$Codispersion, 3)
    out.df$P.value <- round(out.df$P.value, 3)
    out.df$null_mean <- round(out.df$null_mean, 3)
    out.df$Difference <- round(out.df$Difference, 3)
  }
  return(out.df)
}

### Function for simulating a bivariate normal distribution
bivariate <- function(x,y){
  mu1 <- 0     # expected value of x
  mu2 <- 0     # expected value of y
  sig1 <- 1    # variance of x
  sig2 <- 1    # variance of y
  rho <- 0.5   # corr(x, y)
  term1 <- 1 / (2 * pi * sig1 * sig2 * sqrt(1 - rho^2))
  term2 <- (x - mu1)^2 / sig1^2
  term3 <- -(2 * rho * (x - mu1)*(y - mu2))/(sig1 * sig2)
  term4 <- (y - mu2)^2 / sig2^2
  z <- term2 + term3 + term4
  term5 <- term1 * exp((-z / (2 *(1 - rho^2))))
  return (term5)
}


# Function to calculate species diversity indices from a list of ppp objects
# Output is a geoR object

ppp.diversity.to.geoR.fn <- function(ppp.ls = ppp.ls, spp.list = spp.list, xmin, xmax, ymin, ymax, quad.size, index = c("total.abundance", "richness", "shannon", "simpson", "pielou.evenness", "simpson.evenness", "jaccard.dissimilarity", "bray-curtis.dissimilarity") ){
  
  for(i in 1:length(ppp.ls)) {
    
    x <- ppp.ls[[i]]$x # extract x coordinate
    y <- ppp.ls[[i]]$y # extract y coordinate
    xt <- cut(x, seq(xmin, xmax, quad.size)) # cut x coordinates using quad.size spacing
    yt <- cut(y, seq(ymin, ymax, quad.size)) # cut y coordinates using quad.size spacing
    coords <- dimnames(table(yt, xt)) # extract quadrat coordinate lists
    qx <- rep(seq(xmin, xmax - quad.size, length = length(coords$xt)), each = length(coords$yt)) # vector of x coordinates for the bottom left corner of the quadrat
    qy <- rep(seq(ymin, ymax - quad.size, length = length(coords$yt)), length(coords$xt)) # vector of y coordinates for the bottom left corner of the quadrat
    
    out.grid <- table(yt, xt)  # count the trees in each quadrat
    out.grid[is.na(out.grid) == T] <- 0  # replace NAs in table with zeros for empty quadrats
    
    if(i == 1) { spp.df <- data.frame(spp1 = as.vector(out.grid)) }
    if(i > 1)  { spp.df <- cbind(spp.df, as.vector(out.grid)) }
  }
  
  empty.vec <- ifelse(rowSums(spp.df) == 0, 1, 0) # vector of empty quadrats
  
  names(spp.df) <- spp.list  # put species names in dataframe
  
  if(index == "total.abundance") { out.vec <- rowSums(spp.df) }   
  if(index == "richness") { out.vec <- specnumber(spp.df) }  
  if(index == "shannon") { out.vec <- exp(diversity(spp.df, index = "shannon")) }   
  if(index == "simpson") { out.vec <- diversity(spp.df, index = "simpson") 
  #out.vec <- ifelse(out.vec == Inf, 0, out.vec) 
  }  
  if(index == "pielou.evenness") { 
    out.vec <- exp(diversity(spp.df, index = "shannon"))
    out.vec <- ifelse(is.na(out.vec) == T, 0, out.vec) }
  if(index == "simpson.evenness") { 
    out.vec <- diversity(spp.df, index = "invsimpson")
    out.vec <- ifelse(out.vec == Inf, 0, out.vec) }
  if(index == "jaccard.dissimilarity") { 
    out.vec <- colMeans(as.data.frame(as.matrix(vegdist(cbind(spp.df, empty.vec), method = "jaccard", upper = T, diag = T, na.rm = T))))
    out.vec <- ifelse(is.na(out.vec) == T, 1, out.vec)
    print("Note that empty quadrats have the species 'empty' added")   }
  if(index == "bray-curtis.dissimilarity") { 
    out.vec <- colMeans(as.data.frame(as.matrix(vegdist(cbind(spp.df, empty.vec), method = "bray", upper = T, diag = T, na.rm = T))))
    #out.vec <- ifelse(is.na(out.vec) == T, 1, out.vec)
    print("Note that empty quadrats have the species 'empty' added for the calculation of bray-curtis.dissimilarity")  }
  
  out.df <- data.frame(qx, qy, out.vec)
  out.geo <- as.geodata(out.df, coords.col = 1:2, data.col = 3)
  
  return(out.geo)
  
}  # end function


### basal area function: calculates basal area from DBH values (must be in cm)
basal.area.fn <- function(x){ (pi*(x)^2)/40000 } # calculate basal area in m^2


### Function to calculate species importance values (IV) from a list of ppp objects
# IV = relative frequency + relative density + relative basal area (all calculated as percent)
# Output is a dataframe with species names, RF, RD, RBA and IV, sorted by IV

importance.value.fn <- function(ppp.ls = ppp.ls) {
  
  for(i in 1:length(ppp.ls)) {
    
    x <- ppp.ls[[i]]$x # extract x coordinate
    y <- ppp.ls[[i]]$y # extract y coordinate
    z <- ppp.ls[[i]]$marks # extract DBH values
    ba <- (pi * (z)^2) / 40000 # calculate basal area in m^2 
    xt <- cut(x, seq(xmin, xmax, quad.size)) # cut x coordinates using quad.size spacing
    yt <- cut(y, seq(ymin, ymax, quad.size)) # cut y coordinates using quad.size spacing
    coords <- dimnames(table(yt, xt)) # extract quadrat coordinate lists
    qx <- rep(seq(xmin, xmax - quad.size, length = length(coords$xt)), each = length(coords$yt)) # vector of x coordinates for the bottom left corner of the quadrat
    qy <- rep(seq(ymin, ymax - quad.size, length = length(coords$yt)), length(coords$xt)) # vector of y coordinates for the bottom left corner of the quadrat
    
    out.grid <- table(yt, xt)  # abundance
    out.grid[is.na(out.grid) == T] <- 0  # replace NAs in table with zeros for empty quadrats
    
    out.grid.ba <- tapply(ba, list(yt, xt), sum) # calculate total ba in each quadrat    
    out.grid.ba[is.na(out.grid.ba) == T] <- 0 
    
    if(i == 1) { spp.df <- data.frame(spp1 = as.vector(out.grid)) }
    if(i > 1) { spp.df <- cbind(spp.df, as.vector(out.grid)) }

    if(i == 1) { ba.df <- data.frame(spp1 = as.vector(out.grid.ba)) }
    if(i > 1) { ba.df <- cbind(ba.df, as.vector(out.grid.ba)) }
    
  }
  
  names(spp.df) <- all.spp.list  # put species names on dataframe
  names(ba.df) <- all.spp.list  # put species names on dataframe
  
  total.abundance <- sum(spp.df) # total number of all stems in the plot
  total.ba <- sum(ba.df) # total basal area of all stems in the plot
  
  relative.density <- colSums(spp.df) / total.abundance * 100 # percent of stems belonging to each sp
  relative.frequency <- colSums(spp.df > 0) / dim(spp.df)[1] * 100 # percent of occupied cells 
  relative.basal.area <- colSums(ba.df) / total.ba * 100 # percent of basal area belonging to each sp

  out.df <- data.frame(species = all.spp.list, relative.density, relative.frequency, relative.basal.area)
  
  out.df$importance.value <- out.df$relative.density + out.df$relative.frequency + out.df$relative.basal.area

  out.df <- arrange(out.df, desc(importance.value))

  return(out.df)

  } # end IV function


