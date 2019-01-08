## Analyse the relationship between HF foundation species vs. abundance and diversity using codispersion analysis

###################################
#### Load required packages
###################################
library(tidyverse) 
library(geoR)
library(readr)
library(spatstat)
library(vegan) # for calculating diversity indices
library(akima)
library(fields)

###################################
#### Load required functions
###################################
source("Codispersion-analysis-functions.R")
#dyn.load("hCodismap.dll")
dyn.load("hCodismap.so") # on mac, need to set wd to 'mac update' to run this
source("Noises.R")
source("Input.R")
source("Codismap.R")
source("hCodismap_hlb.R")

########################
#### Set the quadrat size 
########################

quad.size = 5
quad.size = 10
quad.size = 20

########################
#### Set the number of randomisations 
########################

n.rand = 199

########################
#### Dataset loop
########################

dataset.names <- c("HF", "WR", "TY", "BCI", "LFDP", "AM")

for(dataset in 1:length(dataset.names)) {
  
  # Massachusetts-HarvardForest
  if(dataset == 1) { 
    HF_raw <- read.csv("hf253-03-trees-2014.csv", header = TRUE)
    lat = 42.5
    xmin = 0; xmax = 700; ymin = 0; ymax = 500  # xmax = 700
    dat <- HF_raw %>% 
      filter(status == "A") %>%
      dplyr::select(sp, gx, gy, dbh) %>%
      distinct(.)
    rm(HF_raw)
    dat$sp <- as.factor(dat$sp)      # make species code a factor
    dat$ba <- basal.area.fn(dat$dbh) # calculate basal area for each tree in the dataset
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size)) # generate vectors for 20 x 20 quadrat grid
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size))
    dat <- dat[ order(dat[,"gx"]), ]
    spp.names.list <- c("Acer rubrum", "Pinus strobus", "Quercus rubra", "Tsuga canadensis")
    spp.list <- c("acerru", "pinust", "querru", "tsugca")
  }
  
  # Washington-Wind River
  if(dataset == 2) {
    lat = 45.82
    xmin = 0; xmax = 800; ymin = 0; ymax = 320  
    windriver_raw <- read.csv("WFDP_Tree_20141108_for_Aaron.csv", header=TRUE)
    dat <- windriver_raw %>% 
      filter(dbh != "NULL", plot_x >= 0, plot_x <= xmax, plot_y <= ymax) %>%
      filter(as.character(tree_tag) == as.character(stem_tag)) %>%
      mutate(dbh = as.numeric(as.character(dbh)) ) %>%  # convert dbh to cm
      dplyr::select(sp = species, gx = plot_x, gy = plot_y, dbh) %>%
      distinct() %>%
      drop_na(.)
    rm(windriver_raw)
    dat$sp <- as.factor(dat$sp)      
    dat$ba <- basal.area.fn(dat$dbh) 
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size), include.lowest = T) 
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size), include.lowest = T)
    dat <- dat[ order(dat[,"gx"]), ]
    spp.names.list <- c("Pseudotsuga menzeisii", "Tsuga heterophylla", "Acer circinatum", "Abies amabilis", "Taxus brevifolia")
    spp.list <- c("PSME", "TSHE", "ACCI", "ABAM", "TABR")
  }
  
  # Missouri-TysonResearchCenter 
  if(dataset == 3) {
    lat = 38.52
    xmin = 0; xmax = 480; ymin = 0; ymax = 420  # xmax = 480
    tyson_raw <- read.csv("TRCP_Census4_20140805_ForAaronEllison.csv", header=TRUE)
    dat <- tyson_raw %>% 
      mutate(gx = gx - 19, gy = gy -19) %>% # plot is offset by 19 m
      filter(stem == "main", status == "alive", dbh != "NULL", gx <= xmax, gy <= ymax) %>%
      mutate(dbh = as.numeric(as.character(dbh)) / 10) %>%  # convert dbh to cm
      dplyr::select(sp = spcode, gx, gy, dbh) %>%
      distinct() %>%
      drop_na(.)
    rm(tyson_raw)
    dat$sp <- as.factor(dat$sp) 
    dat$ba <- basal.area.fn(dat$dbh) 
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size)) 
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size))
    dat <- dat[ order(dat[,"gx"]), ]
    spp.names.list <- c("Quercus alba", "Quercus rubra", "Quercus velutina", "Carya tomentosa", "Carya ovata", "Carya glabra", "Cornus florida", "Lindera benzoin", "Asimina triloba", "Frangula caroliniana")
    spp.list <- c("quealb", "querub", "quevel", "cartom", "carova", "cargla", "corflo", "linben", "asitri", "fracar")
  }
  
  # Panama-BCI
  if(dataset == 4) {
    lat = 9.15
    xmin = 0; xmax = 1000; ymin = 0; ymax = 500  # xmax = 1000
    bci_raw <- read.table("PlotDataReport04-22-2018_census-8.txt", header = T, sep = "\t")
    dat <- bci_raw %>% 
      mutate(dbh = as.numeric(as.character(DBH)) / 10) %>%  # convert dbh to cm
      dplyr::select(sp = Mnemonic, gx = PX, gy = PY, dbh) %>%
      distinct()
    rm(bci_raw)
    dat$sp <- as.factor(dat$sp) 
    dat$ba <- basal.area.fn(dat$dbh) 
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size), include.lowest = T) 
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size), include.lowest = T)
    dat <- dat[ order(dat[,"gx"]), ]
    spp.names.list <- c("Alseis blackiana", "Oenocarpus mapora", "Ficus popenoei")
    spp.list <- c("alsebl", "oenoma", "ficupo")
  }
  
  # PuertoRico-LFDP
  if(dataset == 5) {
    lat = 18.33
    xmin = 0; xmax = 320; ymin = 0; ymax = 500  # ymax = 500
    LFDP_raw1 <- read.csv("LFDP_Census3-Part1.txt", header = T)
    LFDP_raw2 <- read.csv("LFDP_Census3-Part2.txt", header = T)
    LFDP_raw <- rbind(LFDP_raw1, LFDP_raw2)
    dat <- LFDP_raw %>% 
      separate(codes, into = c("stem", "status"), sep = ";", extra = "drop") %>%
      filter(stem == "main", status == "A" | status == "AB", dbh > 0) %>%
      mutate(dbh = dbh / 10) %>%  # convert dbh to cm
      dplyr::select(sp = spcode, gx, gy, dbh) %>%
      drop_na(.)
    rm(LFDP_raw1, LFDP_raw2, LFDP_raw)
    dat$sp <- as.factor(dat$sp)      
    dat$ba <- basal.area.fn(dat$dbh) # calculate basal area for each tree in the dataset
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size)) 
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size))
    dat <- dat[ order(dat[,"gx"]), ]
    spp.names.list <- c("Pisonia subcordata", "Prestoea acuminata", "Cecropia schreberiana", "Dacryodes excelsa")
    spp.list <- c("PISSUB", "PREMON", "CECSCH", "DACEXC")
  }
  
  # Amacayacu 
  if(dataset == 6) {
    lat=-3.80917
    xmin = 0; xmax = 500; ymin = 0; ymax = 500  
    load("amacayacu.full1.RData")
    dat <- amacayacu.full1 %>%
      filter(DFstatus == "alive", sp != "Unidunid") %>% # exclude missing status stems and unidentified species
      mutate(dbh = dbh / 10) %>%  # convert dbh to cm
      dplyr::select(sp, gx, gy, dbh) %>%
      drop_na(.)
    rm(amacayacu.full1)
    dat$sp <- as.factor(dat$sp)    
    dat$ba <- basal.area.fn(dat$dbh) 
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size)) 
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size))
    dat <- dat[ order(dat[,"gx"]), ]
    spp.names.list <- c("Eschweilera coriacea", "Eschweilera itayensis", "Eschweilera rufifolia", "Guarea pubescens", "Otoba glycycarpa", "Rinorea lindeniana")
    spp.list <- c("Eschcori", "Eschitay", "Eschrufi", "Guarpube", "otobglyc", "Rinolind")
  }

########################
#### Set inputs for codispersion analysis and graphing 
########################
#l = xmax/quad.size/(2*4) # plot size / quadrat size / 2 * the number we use to get max lag
l = 4
r = 1
gcolours = colorRampPalette(c("#0000FF","#FFFFFF","#FF6666"))(128)
focal_species_variables <- rep(rep(c("Total\nabundance", "Total basal\narea", "Mean basal\narea"), 5), length(unique(spp.list)))
community_variables <- rep(c(rep("Total\nabundance",3), rep("Richness",3), rep("Shannon",3), rep("Simpson",3), rep("Mean Bray-Curtis\ndissimilarity",3)), length(unique(spp.list)))
  
###################################
## Species loop
###################################

for(species in 1:length(spp.list)) { # loop through each putative FS species
  
  print(paste("quadrat.size =", quad.size, "dataset =", dataset, "species =", spp.list[species]))
  
###################################
## Create input objects for putative foundation species
###################################

putative.FS <- spp.list[species]
  
ppp.sp1 <- ppp(dat$gx[dat$sp == putative.FS], dat$gy[dat$sp == putative.FS], xrange = c(xmin, xmax), yrange = c(ymin, ymax), marks = dat$dbh[dat$sp == putative.FS])

obs.tab.sp1 <- geoR.to.matrix.fn(ppp.to.geoR.fn(ppp.sp1, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "abundance"))
obs.tba.sp1 <- geoR.to.matrix.fn(ppp.to.geoR.fn(ppp.sp1, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "total.ba"))
obs.mba.sp1 <- geoR.to.matrix.fn(ppp.to.geoR.fn(ppp.sp1, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "mean.ba"))

###################################
## Create geodata objects for diversity
###################################

ppp.ls.other <- vector("list", (length(unique(dat$sp)) - 1))  # empty list to hold all OTHER species' point patterns

other.spp.list <- as.character(unique(dat$sp[!(dat$sp %in% putative.FS)]))

for (i in 1:length(other.spp.list)){  # generate ppp objects for all species except FS
  ppp.ls.other[[i]] <- ppp(dat$gx[dat$sp == other.spp.list[i]], dat$gy[dat$sp == other.spp.list[i]], xrange = c(xmin, xmax), yrange = c(ymin, ymax), marks = dat$dbh[dat$sp == other.spp.list[i]])
}

# generate geoR objects for total abundance of all stems and for each diversity measure
obs.tot_ab <- geoR.to.matrix.fn(ppp.diversity.to.geoR.fn(ppp.ls.other, spp.list = other.spp.list, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, index = "total.abundance", quad.size = quad.size))
obs.richness <- geoR.to.matrix.fn(ppp.diversity.to.geoR.fn(ppp.ls.other, spp.list = other.spp.list, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, index = "richness", quad.size = quad.size))
obs.shannon <- geoR.to.matrix.fn(ppp.diversity.to.geoR.fn(ppp.ls.other, spp.list = other.spp.list, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, index = "shannon", quad.size = quad.size))
obs.simpson <- geoR.to.matrix.fn(ppp.diversity.to.geoR.fn(ppp.ls.other, spp.list = other.spp.list, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, index = "simpson", quad.size = quad.size))
obs.braycurtis <- geoR.to.matrix.fn(ppp.diversity.to.geoR.fn(ppp.ls.other, spp.list = other.spp.list, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, index = "bray-curtis.dissimilarity", quad.size = quad.size))


###################################
## Run the codispersion analysis for species diversity data
###################################

cd <- as.list(numeric(15))

cd[[1]] <- hCodismap(obs.tot_ab, r = r, l = l, obs.tab.sp1)
cd[[2]] <- hCodismap(obs.tot_ab, r = r, l = l, obs.tba.sp1)
cd[[3]] <- hCodismap(obs.tot_ab, r = r, l = l, obs.mba.sp1)

cd[[4]] <- hCodismap(obs.richness, r = r, l = l, obs.tab.sp1)
cd[[5]] <- hCodismap(obs.richness, r = r, l = l, obs.tba.sp1)
cd[[6]] <- hCodismap(obs.richness, r = r, l = l, obs.mba.sp1)

cd[[7]] <- hCodismap(obs.shannon, r = r, l = l, obs.tab.sp1)
cd[[8]] <- hCodismap(obs.shannon, r = r, l = l, obs.tba.sp1)
cd[[9]] <- hCodismap(obs.shannon, r = r, l = l, obs.mba.sp1)

cd[[10]] <- hCodismap(obs.simpson, r = r, l = l, obs.tab.sp1)
cd[[11]] <- hCodismap(obs.simpson, r = r, l = l, obs.tba.sp1)
cd[[12]] <- hCodismap(obs.simpson, r = r, l = l, obs.mba.sp1)

cd[[13]] <- hCodismap(obs.braycurtis, r = r, l = l, obs.tab.sp1)
cd[[14]] <- hCodismap(obs.braycurtis, r = r, l = l, obs.tba.sp1)
cd[[15]] <- hCodismap(obs.braycurtis, r = r, l = l, obs.mba.sp1)


######################################################
## Create sets of null randomised matrices
######################################################

rnd.stats <- rep(list(as.list(numeric(n.rand))),5)

for(rand in 1:n.rand) {
  rnd.stats[[1]][[rand]] <- obs.tot_ab[ sample(nrow(obs.tot_ab)), sample(ncol(obs.tot_ab)) ]
  rnd.stats[[2]][[rand]] <- obs.richness[ sample(nrow(obs.richness)), sample(ncol(obs.richness)) ]
  rnd.stats[[3]][[rand]] <- obs.shannon[ sample(nrow(obs.shannon)), sample(ncol(obs.shannon)) ]
  rnd.stats[[4]][[rand]] <- obs.simpson[ sample(nrow(obs.simpson)), sample(ncol(obs.simpson)) ]
  rnd.stats[[5]][[rand]] <- obs.braycurtis[ sample(nrow(obs.braycurtis)), sample(ncol(obs.braycurtis)) ]
}

######################################################
## Compute codispersion for randomised matrices
######################################################

temp <- rep(list(rep(list(as.list(numeric(n.rand))),3)),5)

for(rand in 1:n.rand) {
  
print(paste("randomisation =", rand, "dataset =", dataset, "species =", spp.list[species]))

temp[[1]][[1]][[rand]] <- hCodismap(rnd.stats[[1]][[rand]], r = r, l = l, obs.tab.sp1)$z
temp[[1]][[2]][[rand]] <- hCodismap(rnd.stats[[1]][[rand]], r = r, l = l, obs.tba.sp1)$z
temp[[1]][[3]][[rand]] <- hCodismap(rnd.stats[[1]][[rand]], r = r, l = l, obs.mba.sp1)$z

temp[[2]][[1]][[rand]] <- hCodismap(rnd.stats[[2]][[rand]], r = r, l = l, obs.tab.sp1)$z
temp[[2]][[2]][[rand]] <- hCodismap(rnd.stats[[2]][[rand]], r = r, l = l, obs.tba.sp1)$z
temp[[2]][[3]][[rand]] <- hCodismap(rnd.stats[[2]][[rand]], r = r, l = l, obs.mba.sp1)$z

temp[[3]][[1]][[rand]] <- hCodismap(rnd.stats[[3]][[rand]], r = r, l = l, obs.tab.sp1)$z
temp[[3]][[2]][[rand]] <- hCodismap(rnd.stats[[3]][[rand]], r = r, l = l, obs.tba.sp1)$z
temp[[3]][[3]][[rand]] <- hCodismap(rnd.stats[[3]][[rand]], r = r, l = l, obs.mba.sp1)$z

temp[[4]][[1]][[rand]] <- hCodismap(rnd.stats[[4]][[rand]], r = r, l = l, obs.tab.sp1)$z
temp[[4]][[2]][[rand]] <- hCodismap(rnd.stats[[4]][[rand]], r = r, l = l, obs.tba.sp1)$z
temp[[4]][[3]][[rand]] <- hCodismap(rnd.stats[[4]][[rand]], r = r, l = l, obs.mba.sp1)$z

temp[[5]][[1]][[rand]] <- hCodismap(rnd.stats[[5]][[rand]], r = r, l = l, obs.tab.sp1)$z
temp[[5]][[2]][[rand]] <- hCodismap(rnd.stats[[5]][[rand]], r = r, l = l, obs.tba.sp1)$z
temp[[5]][[3]][[rand]] <- hCodismap(rnd.stats[[5]][[rand]], r = r, l = l, obs.mba.sp1)$z
}

# Extract null model results to a list of arrays
null.out <- as.list(numeric(15))
null.out[[1]] <- list2ary(temp[[1]][[1]])
null.out[[2]] <- list2ary(temp[[1]][[2]])
null.out[[3]] <- list2ary(temp[[1]][[3]])
null.out[[4]] <- list2ary(temp[[2]][[1]])
null.out[[5]] <- list2ary(temp[[2]][[2]])
null.out[[6]] <- list2ary(temp[[2]][[3]])
null.out[[7]] <- list2ary(temp[[3]][[1]])
null.out[[8]] <- list2ary(temp[[3]][[2]])
null.out[[9]] <- list2ary(temp[[3]][[3]])
null.out[[10]] <- list2ary(temp[[4]][[1]])
null.out[[11]] <- list2ary(temp[[4]][[2]])
null.out[[12]] <- list2ary(temp[[4]][[3]])
null.out[[13]] <- list2ary(temp[[5]][[1]])
null.out[[14]] <- list2ary(temp[[5]][[2]])
null.out[[15]] <- list2ary(temp[[5]][[3]])
  
rnd.results <- rep(list(as.list(numeric(3))),5)

for(k in 1:length(cd)){ # loop through observed codispersion results
  
  temp_res <- cd[[k]] # extract observed results
  null.input.ary <- null.out[[k]] # extract null model results
  
  Codispersion_output <- as.vector(temp_res$z)
  Codispersion_output[is.nan(Codispersion_output)] <- 999
  temp_out <- data.frame(cell = 1:length(Codispersion_output), Codispersion = Codispersion_output)
  
# null model comparisons
for(n.col in 1: ncol(null.input.ary)) { # loop through each cell
  for(n.row in 1: nrow(null.input.ary)) {
    null.vec <- null.input.ary[n.row, n.col, ]
    if(n.row == 1) { out <- data.frame(row = n.row, col = n.col, null.cd = null.vec) }
    if(n.row > 1) { out <- rbind(out, data.frame(row = n.row, col = n.col, null.cd = null.vec)) }
    } # end col loop    
    if(n.col == 1) { res.out <- out }
    if(n.col > 1) { res.out <- rbind(res.out, out) }
  } # end row loop
  
res.out$cell <- rep(1:length(unique(paste(res.out$row, res.out$col))), each = n.rand)
  
for(cell.no in 1:length(unique(res.out$cell))) {
    obs.cd.value <- temp_out$Codispersion[temp_out$cell == cell.no]
    null.cd.vec <- res.out$null.cd[res.out$cell == cell.no]
    prop.greater.than <- length(which(null.cd.vec > obs.cd.value)) / n.rand
    prop.less.than <- length(which(null.cd.vec < obs.cd.value)) / n.rand
    temp_out$P.value[cell.no] <- min(prop.greater.than, prop.less.than)
    temp_out$Mean.null.cd[cell.no] <- mean(null.cd.vec, na.rm = T) 
  } # end cell loop
  
  temp_out$Mean.null.cd <- ifelse(temp_out$Codispersion == 999, 999, temp_out$Mean.null.cd)
  temp_out$P.value <- ifelse(temp_out$Codispersion == 999, 999, temp_out$P.value)
  temp_out <- temp_out[complete.cases(temp_out), ]
  temp_out$dataset <- dataset
  temp_out$species = spp.list[species]
  temp_out$quadrat.size = quad.size
  temp_out$xcirc = temp_res$xcirc
  temp_out$ycirc = temp_res$ycirc
  temp_out$Species_name <- putative.FS
  temp_out$Focal_species_variable <- focal_species_variables[k]
  temp_out$Community_variable <- community_variables[k]

  if(k == 1) { out1 <- temp_out }
  if(k > 1) { out1 <- rbind(out1, temp_out) }
  
  } # end k loop
 
if(species == 1) { out2 <- out1 }
if(species > 1) { out2 <- rbind(out2, out1) }

} # end FS species loop

if(dataset == 1) { full.output <- out2 }
if(dataset > 1) { full.output <- rbind(full.output, out2) }

} # end dataset loop

output_5 <- full.output
output_10 <- full.output
output_20 <- full.output

save.image("Results_19-01-07.RData")
#load("Results_19-01-06.RData")

######################################################
## Draw graphs of null model outputs
######################################################

g.output <- output_5 ; quad.size = 5
g.output <- output_10 ; quad.size = 10
g.output <- output_20 ; quad.size = 20

for(d in 1:6) { # dataset loop
  sdat <- g.output[g.output$dataset == d, ]
  for(s in 1:length(unique(sdat$species))){ # species loop
    gdat <- sdat[sdat$species == unique(sdat$species)[s],]
# Observed codispersion plots
gdat %>% 
  mutate(Codisp = ifelse(Codispersion == 999, NA, Codispersion)) %>%
  filter(complete.cases(.)) %>%
  mutate(x = xcirc*quad.size, y = ycirc*quad.size) %>%
  mutate(new_metric = fct_recode(Community_variable, 'Inverse\nSimpson' = "Simpson")) %>%
  ggplot(aes(x = x, y = y, fill = Codisp)) +
  geom_raster(interpolate = F) +
  scale_fill_gradient2(high = "#FF6666", mid = "#FFFFFF", low = "#0000FF", midpoint = 0, limits = c(-1, 1), guide = FALSE) +
  coord_fixed(ratio = 1) + 
  facet_grid(Focal_species_variable ~ new_metric) +
  xlab("Spatial lag in X (m)") + ylab("Spatial lag in Y (m)") +
  theme_bw(base_size = 16)
ggsave(filename = paste("obs", dataset.names[d],unique(sdat$species)[s],quad.size,".jpg", sep="_"), width = 12, height = 6)

# Mean null codispersion plots
gdat %>% 
  mutate(Codisp = ifelse(Codispersion == 999, NA, Codispersion)) %>%
  filter(complete.cases(.)) %>%
  mutate(x = xcirc*quad.size, y = ycirc*quad.size) %>%
  mutate(new_metric = fct_recode(Community_variable, 'Inverse\nSimpson' = "Simpson")) %>%
  ggplot(aes(x = x, y = y, fill = Mean.null.cd)) +
  geom_raster(interpolate = F) +
  scale_fill_gradient2(high = "#FF6666", mid = "#FFFFFF", low = "#0000FF", midpoint = 0, limits = c(-1, 1), guide = FALSE) +
  coord_fixed(ratio = 1) + 
  facet_grid(Focal_species_variable ~ new_metric) +
  xlab("Spatial lag in X (m)") + ylab("Spatial lag in Y (m)") +
  theme_bw(base_size = 16)
ggsave(filename = paste("mean_null", dataset.names[d],unique(sdat$species)[s],quad.size,".jpg", sep="_"), width = 12, height = 6)

# Significance plots
my.cols <- c("steelblue3", "firebrick3")  # select colours for graph
gdat$P.value.cat <- as.factor(ifelse(gdat$P.value < 0.05, "Sig.", "NS"))
if(levels(gdat$P.value.cat)[1] == "Sig.") { my.cols <- c("firebrick3") }
gdat %>% 
  mutate(Codisp = ifelse(Codispersion == 999, NA, Codispersion)) %>%
  mutate(new_metric = fct_recode(Community_variable, 'Inverse\nSimpson' = "Simpson")) %>%
  filter(complete.cases(.)) %>%
  mutate(x = xcirc*quad.size, y = ycirc*quad.size) %>%
  ggplot(aes(x = x, y = y, fill = Mean.null.cd)) +
  geom_tile(aes(fill = P.value.cat)) + 
  scale_fill_manual(values = my.cols) +
  scale_colour_discrete(name = "P.value.cat", limits = c(0,1),guide = FALSE) +
  coord_fixed(ratio = 1) + 
  facet_grid(Focal_species_variable ~ new_metric) +
  xlab("Spatial lag in X (m)") + ylab("Spatial lag in Y (m)") +
  theme_bw(base_size = 16)
ggsave(filename = paste("sig", dataset.names[d],unique(sdat$species)[s],quad.size,".jpg", sep="_"), width = 14, height = 8)
  } # end species graphing null results loop s
} # end dataset graphing null results loop d

######################################################
## SITE graphs for paper (20 x 20 m grids)
######################################################

spp.layer = "overstory"

# Species selection for each site
# HF
dataset = 1; site = "HF"; spp.set <- c("acerru", "pinust", "querru", "tsugca")

# WR
dataset = 2; site = "WR"; spp.set <- c("PSME", "TSHE", "ABAM"); understory.spp <- c("ACCI", "TABR")

# TY
dataset = 3; site = "TY"; spp.set <- c("quealb", "querub", "quevel", "cartom", "carova", "cargla"); understory.spp <- c("asitri", "corflo", "fracar", "linben")

# BCI
dataset = 4; site = "BCI"; spp.set <- c("alsebl", "oenoma")

# LFDP
dataset = 5; site = "LFDP"; spp.set <- c("CECSCH", "DACEXC"); understory.spp <- c("PREMON")

# AM
dataset = 6; site = "AM"; spp.set <- c("Eschcori", "Eschitay", "Eschrufi", "Guarpube", "otobglyc"); understory.spp <- c("Rinolind")

### Observed Codispersion values
g.output %>%
  mutate(new_species = fct_recode(Species_name, `Acer rubrum` = "acerru", `Pinus strobus` = "pinust", `Quercus rubra` = "querru", `Tsuga canadensis` = "tsugca", `Pseudotsuga menzeisii` = "PSME", `Tsuga heterophylla` = "TSHE", `Abies amabilis` = "ABAM", `Acer circinatum` = "ACCI", `Taxus brevifolia` = "TABR", `Quercus alba` = "quealb", `Quercus rubra` = "querub", `Quercus velutina` = "quevel", `Carya tomentosa` = "cartom", `Carya ovata` = "carova", `Carya glabra` = "cargla", `Asimina triloba` = "asitri", `Cornus florida` = "corflo", `Frangula caroliniana` = "fracar", `Lindera benzoin` = "linben", `Alseis blackiana` = "alsebl", `Oenocarpus mapora` = "oenoma", `Pisonia subcordata` = "PISSUB", `Cecropia schreberiana` = "CECSCH", `Dacryodes excelsa` = "DACEXC", `Prestoea acuminata` = "PREMON", `Eschweilera coriacea` = "Eschcori", `Eschweilera itayensis` = "Eschitay", `Eschweilera rufifolia` = "Eschrufi", `Guarea pubescens` = "Guarpube", `Otoba glycycarpa` = "otobglyc", `Rinorea lindeniana` = "Rinolind")) %>%
  filter(dataset == dataset) %>%
  filter(Focal_species_variable == "Total basal\narea") %>%
  filter(Species_name %in% spp.set) %>%
  mutate(new_metric = fct_recode(Community_variable, 'Inverse\nSimpson' = "Simpson")) %>%
  mutate(Codisp = ifelse(Codispersion == 999, NA, Codispersion)) %>%
  filter(complete.cases(.)) %>%
  mutate(x = xcirc*quad.size, y = ycirc*quad.size) %>%
  ggplot(aes(x = x, y = y, fill = Codisp)) +
  geom_raster(interpolate = F) +
  scale_fill_gradient2(high = "#FF6666", mid = "#FFFFFF", low = "#0000FF", midpoint = 0, limits = c(-1, 1), guide = FALSE) +
  coord_fixed(ratio = 1) + 
  facet_grid(new_metric ~ new_species) +
  xlab("Spatial lag in X (m)") + ylab("Spatial lag in Y (m)") +
  theme_bw(base_size = 16) +
  theme(strip.text.x = element_text(face = "italic"))
ggsave(filename = paste(site, "Total_basal_area", spp.layer, quad.size, "codisp", ".jpg", sep="_"), width = 20, height = 10)

### Significance graphs (20 x 20 m grids)
my.cols <- c("steelblue3", "firebrick3")  # select colours for graph
g.output$P.value.cat <- as.factor(ifelse(g.output$P.value < 0.05, "Sig.", "NS"))
if(levels(g.output$P.value.cat)[1] == "Sig.") { my.cols <- c("firebrick3") }
g.output %>%
  mutate(new_species = fct_recode(Species_name, `Acer rubrum` = "acerru", `Pinus strobus` = "pinust", `Quercus rubra` = "querru", `Tsuga canadensis` = "tsugca", `Pseudotsuga menzeisii` = "PSME", `Tsuga heterophylla` = "TSHE", `Abies amabilis` = "ABAM", `Acer circinatum` = "ACCI", `Taxus brevifolia` = "TABR", `Quercus alba` = "quealb", `Quercus rubra` = "querub", `Quercus velutina` = "quevel", `Carya tomentosa` = "cartom", `Carya ovata` = "carova", `Carya glabra` = "cargla", `Asimina triloba` = "asitri", `Cornus florida` = "corflo", `Frangula caroliniana` = "fracar", `Lindera benzoin` = "linben", `Alseis blackiana` = "alsebl", `Oenocarpus mapora` = "oenoma", `Pisonia subcordata` = "PISSUB", `Cecropia schreberiana` = "CECSCH", `Dacryodes excelsa` = "DACEXC", `Prestoea acuminata` = "PREMON", `Eschweilera coriacea` = "Eschcori", `Eschweilera itayensis` = "Eschitay", `Eschweilera rufifolia` = "Eschrufi", `Guarea pubescens` = "Guarpube", `Otoba glycycarpa` = "otobglyc", `Rinorea lindeniana` = "Rinolind")) %>%
  filter(dataset == dataset) %>%
  filter(Focal_species_variable == "Total basal\narea") %>%
  filter(Species_name %in% spp.set) %>%
  mutate(new_metric = fct_recode(Community_variable, 'Inverse\nSimpson' = "Simpson")) %>%
  mutate(Codisp = ifelse(Codispersion == 999, NA, Codispersion)) %>%
  filter(complete.cases(.)) %>%
  mutate(x = xcirc*quad.size, y = ycirc*quad.size) %>%
  ggplot(aes(x = x, y = y, fill = Mean.null.cd)) +
  geom_tile(aes(fill = P.value.cat)) + 
  scale_fill_manual(values = my.cols) +
  scale_colour_discrete(name = "P.value.cat", limits = c(0,1),guide = FALSE) +
  coord_fixed(ratio = 1) + 
  facet_grid(new_metric ~ new_species) +
  xlab("Spatial lag in X (m)") + ylab("Spatial lag in Y (m)") +
  theme_bw(base_size = 16) +
  theme(strip.text.x = element_text(face = "italic"), legend.position = "none")
ggsave(filename = paste(site, "Total_basal_area", spp.layer, quad.size, "sig", ".jpg", sep="_"), width = 20, height = 10)

#### Understory graphs: Temperate vs. tropical

spp.set <- c("asitri", "corflo", "fracar", "linben", "ACCI", "TABR"); TT = "Temperate"
spp.set <- c("PREMON", "Rinolind"); TT = "Tropical"

## Codispersion graphs
g.output %>%
  mutate(new_species = fct_recode(Species_name, `Acer rubrum` = "acerru", `Pinus strobus` = "pinust", `Quercus rubra` = "querru", `Tsuga canadensis` = "tsugca", `Pseudotsuga menzeisii` = "PSME", `Tsuga heterophylla` = "TSHE", `Abies amabilis` = "ABAM", `Acer circinatum` = "ACCI", `Taxus brevifolia` = "TABR", `Quercus alba` = "quealb", `Quercus rubra` = "querub", `Quercus velutina` = "quevel", `Carya tomentosa` = "cartom", `Carya ovata` = "carova", `Carya glabra` = "cargla", `Asimina triloba` = "asitri", `Cornus florida` = "corflo", `Frangula caroliniana` = "fracar", `Lindera benzoin` = "linben", `Alseis blackiana` = "alsebl", `Oenocarpus mapora` = "oenoma", `Pisonia subcordata` = "PISSUB", `Cecropia schreberiana` = "CECSCH", `Dacryodes excelsa` = "DACEXC", `Prestoea acuminata` = "PREMON", `Eschweilera coriacea` = "Eschcori", `Eschweilera itayensis` = "Eschitay", `Eschweilera rufifolia` = "Eschrufi", `Guarea pubescens` = "Guarpube", `Otoba glycycarpa` = "otobglyc", `Rinorea lindeniana` = "Rinolind")) %>%
  #filter(dataset == 2 | dataset == 3) %>% # temperate
  filter(dataset == 5 | dataset == 6) %>% # tropical
  #mutate(newer_species = factor(new_species, levels=c("Acer circinatum", "Taxus brevifolia", "Asimina triloba", "Cornus florida", "Frangula caroliniana", "Lindera benzoin"))) %>%
  filter(Focal_species_variable == "Total basal\narea") %>%
  filter(Species_name %in% spp.set) %>%
  mutate(new_metric = fct_recode(Community_variable, 'Inverse\nSimpson' = "Simpson")) %>%
  mutate(Codisp = ifelse(Codispersion == 999, NA, Codispersion)) %>%
  filter(complete.cases(.)) %>%
  mutate(x = xcirc*quad.size, y = ycirc*quad.size) %>%
  ggplot(aes(x = x, y = y, fill = Codisp)) +
  geom_raster(interpolate = F) +
  scale_fill_gradient2(high = "#FF6666", mid = "#FFFFFF", low = "#0000FF", midpoint = 0, limits = c(-1, 1), guide = FALSE) +
  coord_fixed(ratio = 1) + 
  facet_grid(new_metric ~ new_species) +  # newer_species
  xlab("Spatial lag in X (m)") + ylab("Spatial lag in Y (m)") +
  theme_bw(base_size = 16) +
  theme(strip.text.x = element_text(face = "italic"))
ggsave(filename = paste(TT, "understory", "Total_basal_area", quad.size, "codisp", ".jpg", sep="_"), width = 20, height = 10)

### Significance graphs (20 x 20 m grids)
my.cols <- c("steelblue3", "firebrick3")  # select colours for graph
g.output$P.value.cat <- as.factor(ifelse(g.output$P.value < 0.05, "Sig.", "NS"))
if(levels(g.output$P.value.cat)[1] == "Sig.") { my.cols <- c("firebrick3") }
g.output %>%
  mutate(new_species = fct_recode(Species_name, `Acer rubrum` = "acerru", `Pinus strobus` = "pinust", `Quercus rubra` = "querru", `Tsuga canadensis` = "tsugca", `Pseudotsuga menzeisii` = "PSME", `Tsuga heterophylla` = "TSHE", `Abies amabilis` = "ABAM", `Acer circinatum` = "ACCI", `Taxus brevifolia` = "TABR", `Quercus alba` = "quealb", `Quercus rubra` = "querub", `Quercus velutina` = "quevel", `Carya tomentosa` = "cartom", `Carya ovata` = "carova", `Carya glabra` = "cargla", `Asimina triloba` = "asitri", `Cornus florida` = "corflo", `Frangula caroliniana` = "fracar", `Lindera benzoin` = "linben", `Alseis blackiana` = "alsebl", `Oenocarpus mapora` = "oenoma", `Pisonia subcordata` = "PISSUB", `Cecropia schreberiana` = "CECSCH", `Dacryodes excelsa` = "DACEXC", `Prestoea acuminata` = "PREMON", `Eschweilera coriacea` = "Eschcori", `Eschweilera itayensis` = "Eschitay", `Eschweilera rufifolia` = "Eschrufi", `Guarea pubescens` = "Guarpube", `Otoba glycycarpa` = "otobglyc", `Rinorea lindeniana` = "Rinolind")) %>%
  #filter(dataset == 2 | dataset == 3) %>% # temperate
  filter(dataset == 5 | dataset == 6) %>% # tropical
  #mutate(newer_species = factor(new_species, levels=c("Acer circinatum", "Taxus brevifolia", "Asimina triloba", "Cornus florida", "Frangula caroliniana", "Lindera benzoin"))) %>%
  filter(Focal_species_variable == "Total basal\narea") %>%
  filter(Species_name %in% spp.set) %>%
  mutate(new_metric = fct_recode(Community_variable, 'Inverse\nSimpson' = "Simpson")) %>%
  mutate(Codisp = ifelse(Codispersion == 999, NA, Codispersion)) %>%
  filter(complete.cases(.)) %>%
  mutate(x = xcirc*quad.size, y = ycirc*quad.size) %>%
  ggplot(aes(x = x, y = y, fill = Mean.null.cd)) +
  geom_tile(aes(fill = P.value.cat)) + 
  scale_fill_manual(values = my.cols) +
  scale_colour_discrete(name = "P.value.cat", limits = c(0,1),guide = FALSE) +
  coord_fixed(ratio = 1) + 
  facet_grid(new_metric ~ new_species) + # newer_species
  xlab("Spatial lag in X (m)") + ylab("Spatial lag in Y (m)") +
  theme_bw(base_size = 16) +
  theme(strip.text.x = element_text(face = "italic"), legend.position = "none")
ggsave(filename = paste(TT, "understory", "Total_basal_area", quad.size, "sig", ".jpg", sep="_"), width = 20, height = 10)


######################################################
## Codispersion results table
######################################################

# site, no.stems, mean DBH, total BA, mean (SD) codisp, Range in codisp (min, max)
tab1 <- g.output %>%
  filter(Species_name != "PISSUB", Species_name != "ficupo") %>%
  filter(Focal_species_variable == "Total basal\narea") %>%
  mutate(new_species = fct_recode(Species_name, `Eschweilera coriacea` = "Eschcori", `Eschweilera itayensis` = "Eschitay", `Eschweilera rufifolia` = "Eschrufi", `Guarea pubescens` = "Guarpube", `Otoba glycycarpa` = "otobglyc", `Rinorea lindeniana` = "Rinolind",  `Cecropia schreberiana` = "CECSCH", `Dacryodes excelsa` = "DACEXC", `Prestoea acuminata` = "PREMON", `Alseis blackiana` = "alsebl", `Oenocarpus mapora` = "oenoma", `Quercus alba` = "quealb", `Quercus rubra` = "querub", `Quercus velutina` = "quevel", `Carya tomentosa` = "cartom", `Carya ovata` = "carova", `Carya glabra` = "cargla", `Asimina triloba` = "asitri", `Cornus florida` = "corflo", `Frangula caroliniana` = "fracar", `Lindera benzoin` = "linben", `Pseudotsuga menzeisii` = "PSME", `Tsuga heterophylla` = "TSHE", `Abies amabilis` = "ABAM", `Acer circinatum` = "ACCI", `Taxus brevifolia` = "TABR", `Acer rubrum` = "acerru", `Pinus strobus` = "pinust", `Quercus rubra` = "querru", `Tsuga canadensis` = "tsugca")) %>%
  mutate(new_dataset = fct_recode(as.factor(dataset), "HF" = "1", "WR" = "2", "TY" = "3", "BCI" = "4", "LFDP" = "5", "AM" = "6")) %>%
  mutate(latitude = fct_recode(as.factor(dataset), "43.8" = "1", "42.5" = "2", "38.52" = "3", "9.15" = "4", "18.33" = "5", "-3.80917" = "6")) %>%
  mutate(Codisp = ifelse(Codispersion == 999, NA, Codispersion)) %>%
  group_by(new_dataset, new_species, Community_variable) %>%
  summarise(mean_codisp = round(mean(Codisp, na.rm = T), 2), sd_codisp = round(sd(Codisp, na.rm = T), 2), med_codisp = round(median(Codisp, na.rm = T), 2), min_codispersion = round(min(Codisp, na.rm = T), 2), max_codisp = round(max(Codisp, na.rm = T), 2), lat = mean(as.numeric(as.character(latitude))))
  
write.csv(tab1, "Codisp_results_table_new_all.csv", row.names = F)







