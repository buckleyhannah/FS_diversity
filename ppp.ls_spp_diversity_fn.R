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