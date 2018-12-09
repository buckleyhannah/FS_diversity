## FS diversity paper figures for PIs on FS

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
library(ggrepel)

###################################
#### Load required functions
###################################
source("Codispersion-analysis-functions.R")
source("ppp.ls_spp_diversity_fn.R")

########################
#### Select the quadrat size 
########################

quad.size = 20

########################
#### Dataset loop
########################

dataset.names <- c("HF", "WR", "TY", "BCI", "LFDP", "AM")
all.FS.spp <- c("acerru", "pinust", "querru", "tsugca", "PSME", "TSHE", "quealb", "querub", "quevel", "cartom", "carova", "cargla", "alsebl", "oenoma", "ficupo", "PISSUB", "PREMON", "CECSCH", "DACEXC", "Eschcori", "Eschitay", "Eschrufi", "Guarpube", "otobglyc", "Rinolind")
all.FS.names <- c("Acer rubrum", "Pinus strobus", "Quercus rubra", "Tsuga canadensis", "Pseudotsuga menzeisii", "Tsuga heterophylla", "Quercus alba", "Quercus rubra", "Quercus velutina", "Carya tomentosa", "Carya ovata", "Carya glabra", "Alseis blackiana", "Oenocarpus mapora", "Ficus popenoei", "Pisonia subcordata", "Prestoea acuminata", "Cecropia schreberiana", "Dacryodes excelsa", "Eschweilera coriacea", "Eschweilera itayensis", "Eschweilera rufifolia", "Guarea pubescens", "Otoba glycarpa", "Rinorea lindeniana")

for(dataset in 1:length(dataset.names)) {
  
  if(dataset == 1) {
    # 1. Massachusetts-HarvardForest
    #HF_raw <- read.csv("http://harvardforest.fas.harvard.edu/data/p25/hf253/hf253-03-trees-2014.csv", header = TRUE)
    HF_raw <- read.csv("hf253-03-trees-2014.csv", header = TRUE)
    lat = 42.5
    xmin = 0; xmax = 700; ymin = 0; ymax = 500
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
    
    # subset data to create square plot
    #dat <- dat[dat$gx < 500, ]
    spp.names.list <- c("Acer rubrum", "Pinus strobus", "Quercus rubra", "Tsuga canadensis")
    spp.list <- c("acerru", "pinust", "querru", "tsugca")
  }
  
  # 2. California-Yosemite
  # yosemite_full <- read.csv("../data/California-Yosemite/YFDP_Tree_20141108_for_Aaron.csv", header=TRUE)
  # names(yosemite_full) <- c("stem_tag","tree_tag","quadrat","sp","dbh","gx","gy") # make names consistent with HF plot
  # yosemite <- subset(yosemite_full,gx>=0&gx<=800&gy>=0&gy<=320) # limit trees in the dataset to those rooted in the 800 x 300m plot (removes 42 stems)
  
  if(dataset == 2) {
    # 3. Washington-Wind River
    #spp.list <- c("ACCI","ABAM","PSME","TABR","THPL","TSHE");lat=43.8;long=121.9558;area=25.6;grid.loc=2;gg.loc=2;dataset_code="WR";focal.spp=c("PSME","TSHE");focal=c("P. menziesii","T. heterophylla")
    lat = 43.8
    xmin = 0; xmax = 320; ymin = 0; ymax = 320  
    windriver_raw <- read.csv("WFDP_Tree_20141108_for_Aaron.csv", header=TRUE)
    dat <- windriver_raw %>% 
      filter(dbh != "NULL", plot_x >= 0, plot_x <= xmax, plot_y <= xmax) %>%
      mutate(dbh = as.numeric(as.character(dbh)) / 10) %>%  # convert dbh to cm
      dplyr::select(sp = species, gx = plot_x, gy = plot_y, dbh) %>%
      distinct() %>%
      drop_na(.)
    rm(windriver_raw)
    dat$sp <- as.factor(dat$sp)      
    dat$ba <- basal.area.fn(dat$dbh) 
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size)) 
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size))
    dat <- dat[ order(dat[,"gx"]), ]
    spp.names.list <- c("Pseudotsuga menzeisii", "Tsuga heterophylla")
    spp.list <- c("PSME", "TSHE")
  }
  
  if(dataset == 3) {
    # 4. Missouri-TysonResearchCenter 
    #dat <- tyson; dataset_title="Tyson Research Center, Missouri";spp.list <- c("corflo","fracar","linben","quealb","querub","quesp","quevel");lat=38.52;long=-90.6;area=22.05;grid.loc=4;gg.loc=4;dataset_code="TY";focal.spp=c("quealb","querub","quevel");focal=c("Q. alba","Q. rubra","Q. velutina")
    lat = 38.52
    xmin = 0; xmax = 420; ymin = 0; ymax = 420; xmax = 480
    tyson_raw <- read.csv("TRCP_Census4_20140805_ForAaronEllison.csv", header=TRUE)
    dat <- tyson_raw %>% 
      mutate(gx = gx - 19, gy = gy -19) %>% # plot is offset by 19 m
      filter(status == "alive" | stem == "main", dbh != "NULL", gx <= xmax, gy <= ymax) %>%
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
    # subset data to create square plot
    #dat <- dat[dat$gx < 420, ]
    spp.names.list <- c("Quercus alba", "Quercus rubra", "Quercus velutina", "Carya tomentosa", "Carya ovata", "Carya glabra")
    spp.list <- c("quealb", "querub", "quevel", "cartom", "carova", "cargla")
  }
  
  # 5. Michigan-ESGeorgeReserve
  # esgeorge_full <- read.csv("../data/Michigan-ESGeorgeReserve/esgr2008.csv", header=TRUE)
  # esgeorge <- esgeorge_full[esgeorge_full$x<=300&esgeorge_full$x>=0,] # exclude extra plot sections
  # names(esgeorge) <- c("sp","gx","gy","dbh")
  # esgeorge$sp <- factor(esgeorge$sp)
  
  # 6. Virginia-SCBI
  # SCBI_full <- read.csv("../data/Virginia-SCBI/SCBI_initial_woody_stem_census_2012.csv", header=TRUE)
  # names(SCBI_full) <- c("Plot","sp","Latin","Quadrat","gx","gy","TreeID","Tag","StemID","Census","dbh","HOM","Date","Codes","Stem","Status")
  # SCBI_sub1 <- SCBI_full[SCBI_full$Status=="alive",] # exclude dead trees
  # SCBI <- SCBI_sub1[SCBI_sub1$dbh>0,] # exclude records for which dbh is "-999"
  # SCBI$dbh <- SCBI$dbh/10 # convert dbh from mm to cm
  # SCBI$sp <- factor(SCBI$sp)
  
  if(dataset == 4) {
    # 7. Panama-BCI
    #spp.list <- c("Alseis blackiana","Anacardium excelsum","Apeiba membranacea","Ceiba pentandra","Faramea occidentalis","Hura crepitans","Hybanthus prunifolius","Jacaranda copaia","Trichilia tuberculata")
    lat = 9.15
    xmin = 0; xmax = 500; ymin = 0; ymax = 500; xmax = 1000
    bci_raw <- read.table("PlotDataReport04-22-2018_census-8.txt", header = T, sep = "\t")
    dat <- bci_raw %>% 
      filter(!is.na(PX)) %>% 
      filter(Latin != "Unidentified species") %>%
      mutate(dbh = as.numeric(as.character(DBH)) / 10) %>%  # convert dbh to cm
      dplyr::select(sp = Mnemonic, gx = PX, gy = PY, dbh) %>%
      distinct()
    rm(bci_raw)
    dat$sp <- as.factor(dat$sp) 
    dat$ba <- basal.area.fn(dat$dbh) 
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size)) 
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size))
    dat <- dat[ order(dat[,"gx"]), ]
    # subset data to create square plot
    #dat <- dat[dat$gx < 500, ]
    spp.names.list <- c("Alseis blackiana", "Oenocarpus mapora", "Ficus popenoei")
    spp.list <- c("alsebl", "oenoma", "ficupo")
  }
  
  # 8. Panama-Cocoli
  # ## this has multiple stems....need to figure out this
  # cocoli_full <- read.csv("../data/Panama-Cocoli/cocoli.csv", header=TRUE)
  # cocoli_sub1 <- cocoli_full[cocoli_full$recr3=="A",] # exclude dead trees
  # cocoli_sub <- subset(cocoli_sub1,select=c(spcode,x,y,dbh3))
  # names(cocoli_sub) <- c("sp","gx","gy","dbh")
  # cocoli_sub1 <- cocoli_sub[cocoli_sub$dbh>0,]  # exclude records for which dbh is "-2" (plants not yet in the census)
  # cocoli <- cocoli_sub1[cocoli_sub1$gx<=100,] # exclude extra plot section
  # cocoli$dbh <- as.numeric(as.character(cocoli$dbh))/10 # convert to cm
  # cocoli$sp <- factor(cocoli$sp)
  
  # 9. Panama-Sherman
  # sherman_full <- read.table("../data/Panama-Sherman/sherman.txt", sep="\t", header=TRUE)
  # names(sherman_full) <- c("tag","sp","gx","gy","dbh1","dbh2","dbh3","recr1","recr2","recr3","pom1","pom2","pom3","code1","code2","code3","mult1","mult2","mult3","date1","date2","date3")
  # sherman_sub <- sherman_full[sherman_full$sp!="*",]   # remove species coded as "*"
  # sherman_sub$dbh <- sherman_sub$dbh3/10
  # sherman_sub1 <- sherman_sub[sherman_sub$recr3=="A",] # exclude dead trees
  # sherman <- subset(sherman_sub1,gx>=140)   # remove extra plot data from X = 0-140, Y = 0-140. This leaves a 100 x 400 m plot
  # sherman$gx=sherman$gx-140  # make coordinates start at (0,0) 
  # sherman$gy=sherman$gy-40
  # sherman$sp <- factor(sherman$sp)
  
  if(dataset == 5) {
    # 10. PuertoRico-LFDP
    lat = 18.33
    xmin = 0; xmax = 320; ymin = 0; ymax = 320; ymax = 500
    LFDP_raw1 <- read.csv("LFDP_Census3-Part1.txt", header = T)
    LFDP_raw2 <- read.csv("LFDP_Census3-Part2.txt", header = T)
    LFDP_raw <- rbind(LFDP_raw1, LFDP_raw2)
    dat <- LFDP_raw %>% 
      separate(codes, into = c("temp", "status"), sep = ";", extra = "drop") %>%
      filter(status == "A" | status == "AB", dbh > 0) %>%
      mutate(dbh = dbh / 10) %>%  # convert dbh to cm
      dplyr::select(sp = spcode, gx, gy, dbh) %>%
      drop_na(.)
    rm(LFDP_raw1, LFDP_raw2, LFDP_raw)
    dat$sp <- as.factor(dat$sp)      
    dat$ba <- basal.area.fn(dat$dbh) # calculate basal area for each tree in the dataset
    dat$xt <- cut(dat$gx, seq(0, round(max(dat$gx), 0), quad.size)) 
    dat$yt <- cut(dat$gy, seq(0, round(max(dat$gy), 0), quad.size))
    dat <- dat[ order(dat[,"gx"]), ]
    # subset data to create square plot
    #dat <- dat[dat$gy < 320, ]
    spp.names.list <- c("Pisonia subcordata", "Prestoea acuminata", "Cecropia schreberiana", "Dacryodes excelsa")
    spp.list <- c("PISSUB", "PREMON", "CECSCH", "DACEXC")
  }
  
  if(dataset == 6) {
    #11. Amacayacu 
    #dataset_title="Amacayacu, Colombia"; spp.list <- c("Eschcori","Eschitay","Eschrufi","Guarpube","otobglyc","Rinolind"); lat=-3.80917; area=25; grid.loc=11; gg.loc=11; dataset_code="AM"; focal.spp=c("Eschcori","Eschitay","Eschrufi","Guarpube","otobglyc","Rinolind"); focal=c("E. coriacea","E. itayensis","E. rufifolia","G. pubescens","O. glycycarpa","R. lindeniana")  
    #AM.spp = c("Eschcori","Eschitay","Eschrufi","Guarpube","otobglyc","Rinolind")
    # Eschweilera coriacea, Rinorea lindeniana, Guarea pubescens, Otoba glycycarpa, Eschweilera itayensis, Eschweilera rufifolia
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
    spp.names.list <- c("Eschweilera coriacea", "Eschweilera itayensis", "Eschweilera rufifolia", "Guarea pubescens", "Otoba glycarpa", "Rinorea lindeniana")
    spp.list <- c("Eschcori", "Eschitay", "Eschrufi", "Guarpube", "otobglyc", "Rinolind")
  }
  
  ########################
  #### Set inputs for graphing 
  ########################
  gcolours = colorRampPalette(c("#0000FF","#FFFFFF","#FF6666"))(128)
  focal_species_variables <- rep(rep(c("Total\nabundance", "Total basal\narea", "Mean basal\narea"), 5), length(unique(spp.list)))
  community_variables <- rep(c(rep("Total\nabundance",3), rep("Richness",3), rep("Shannon",3), rep("Simpson",3), rep("Mean Bray-Curtis\ndissimilarity",3)), length(unique(spp.list)))
  
  print(paste("quadrat.size =", quad.size, "dataset =", dataset))
    
###################################
## Create input objects for the set of putative foundation species
###################################
    
ppp.FS <- ppp(dat$gx[dat$sp %in% spp.list], dat$gy[dat$sp %in% spp.list], xrange = c(xmin, xmax), yrange = c(ymin, ymax), marks = dat$dbh[dat$sp %in% spp.list])

obs.tab.sp1 <- geoR.to.matrix.fn(ppp.to.geoR.fn(ppp.FS, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "abundance"))
obs.tba.sp1 <- geoR.to.matrix.fn(ppp.to.geoR.fn(ppp.FS, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "total.ba"))
obs.mba.sp1 <- geoR.to.matrix.fn(ppp.to.geoR.fn(ppp.FS, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "mean.ba"))

site <- dataset.names[dataset]

for(i in 1:length(spp.list)) {
  
  species.name <- spp.list[i]
  
  ppp.FS <- ppp(dat$gx[dat$sp == spp.list[i]], dat$gy[dat$sp == spp.list[i]], xrange = c(xmin, xmax), yrange = c(ymin, ymax), marks = dat$dbh[dat$sp == spp.list[i]])
  obs.tab.sp1 <- ppp.to.geoR.fn(ppp.FS, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "abundance")
  obs.tba.sp1 <- ppp.to.geoR.fn(ppp.FS, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "total.ba")
  obs.mba.sp1 <- ppp.to.geoR.fn(ppp.FS, quad.size = quad.size, xmin, xmax, ymin, ymax, method = "mean.ba")
  
  temp <- data.frame(Dataset = site, Latitude = lat, Species = species.name, x = obs.tab.sp1$coords[,1], y = obs.tab.sp1$coords[,2], Abundance = obs.tab.sp1$data, Basal.area = obs.tba.sp1$data, Mean.basal.area = obs.mba.sp1$data)
  
  if(i == 1) { out <- temp }
  if(i > 1) { out <- rbind(temp, out) }
} # end of FS loop

all.spp.list <- unique(dat$sp)

out1 <- dat %>% 
  dplyr::select(sp, dbh, ba) %>%
  group_by(sp) %>%
  summarise(Abundance = n(), Mean.DBH = mean(dbh), Total.BA = sum(ba), Mean.ba = mean(ba)) %>%
  mutate(Dataset = site, Latitude = lat, FS = ifelse(sp %in% all.FS.spp, 1, 0))

out1$FS <- ifelse(out1$sp %in% all.FS.spp, 1, 0)

dat$Dataset <- site
dat$Latitude <- lat

if(dataset == 1) { FS.gdat <- out }
if(dataset > 1) { FS.gdat <- rbind(FS.gdat, out) }

if(dataset == 1) { AB.gdat <- out1 }
if(dataset > 1) { AB.gdat <- rbind(AB.gdat, out1) }

if(dataset == 1) { ALL.dat <- dat }
if(dataset > 1) { ALL.dat <- rbind(ALL.dat, dat) }

} # end of dataset loop

table(FS.gdat$Dataset)
table(FS.gdat$Species)
table(FS.gdat$Dataset, FS.gdat$Species)

table(AB.gdat$Dataset)
table(AB.gdat$sp)
table(AB.gdat$Dataset, AB.gdat$FS)

table(ALL.dat$sp)

######################################################
## Graph of total abundance and total basal area per 20 x 20 m
######################################################

# HF, WR, TY, BCI, LFDP, AM
cols <- c("firebrick3","seagreen4","steelblue3","purple",  
          "firebrick3","seagreen4",
          "firebrick3","seagreen4","steelblue3","purple","violetred3","black",
          "firebrick3","seagreen4","steelblue3", 
          "firebrick3","seagreen4","steelblue3","purple",  
          "firebrick3","seagreen4","steelblue3","purple","violetred3","black" )

cols1 <- c("firebrick3","seagreen4","steelblue3","purple")

FS.gdat %>%
  mutate(dataset_f = factor(Dataset, levels=c("WR","HF","TY","LFDP","BCI","AM"))) %>%
  ggplot(aes(x = Abundance, y = Basal.area, col = Species)) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = cols) +
  facet_wrap(. ~ dataset_f, scales = "free", ncol=2) +
  xlab("Number of individuals per 20 x 20-m plot") + ylab("Total basal area (m^2) per 20 x 20-m plot") +
  theme_bw(base_size = 16) #+ 
  #theme(legend.position = "none")
ggsave("AB-BA_FS_20m.png", width = 8, height = 11, units="in", dpi = 300)


HF <- FS.gdat[FS.gdat$Dataset=="HF",]
 ggplot(HF, aes(x = Abundance, y = Basal.area, col = Species)) +
  geom_point(size = 1) +
  scale_colour_manual(values = cols1) +
#  facet_wrap(. ~ dataset_f, scales = "free", ncol=2) +
  xlab("Number of individuals per 20 x 20-m plot") + ylab("Total basal area (m^2) per 20 x 20-m plot") +
  theme_bw(base_size = 16) #+ 
  #theme(legend.position = "none")
ggsave("Fig2-HF-AB-BA_20m.png", width = 5, height = 5, units="in", dpi = 300)

######################################################
## Graph of mean DBH and abundance for all species
######################################################
FS.cols <- c("black", "red")
  
AB.ggdat <- AB.gdat %>%
  #filter(FS == 1) %>%
  mutate(FS = as.factor(FS)) %>%
  mutate(dataset_f = factor(Dataset, levels=c("WR","HF","TY","LFDP","BCI","AM")))

ggplot(AB.ggdat, aes(x = Mean.DBH, y = Abundance, col = FS)) +
  geom_point(size = 0.5) +
  geom_text_repel(data = subset(AB.ggdat, FS == "1"), aes(label = sp), show.legend = FALSE) +
  scale_colour_manual(values = FS.cols) +
  facet_wrap(. ~ dataset_f, scales = "free", ncol=2) +
  xlab("Mean DBH per tree (cm)") + ylab("Total number of individuals") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none")
ggsave("AB-dbh_AS.png", width = 8, height = 11, units="in", dpi = 300)

ggplot(HF.gdat.f, aes(x = Mean.DBH, y = Abundance, col = as.factor(FS))) +
  geom_point(size = 1) +
  geom_text_repel(data = subset(HF.gdat.f, FS == "1"), aes(label = sp), show.legend = FALSE) +
  scale_colour_manual(values = FS.cols) +
#  facet_wrap(. ~ dataset_f, scales = "free", ncol=2) +
  xlab("Mean DBH per tree (cm)") + ylab("Total number of individuals") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none")
ggsave("Fig1-HF-dbh_abundance.png", width = 5, height = 5, units="in", dpi = 300)




######################################################
## Map of trees where size is proportional to the BA
######################################################
# graphs for FS
ALL.dat %>%
  filter(sp %in% all.FS.spp) %>%
  mutate(dataset_f = factor(dataset, levels=c("WR","HF","TY","LFDP","BCI","AM"))) %>%
  dplyr::select(dataset_f, Latitude, gx, gy, sp, ba) %>%
  ggplot(aes(x = gx, y = gy, fill = sp, size = ba)) +
  geom_point(stat = "identity") +
  #scale_size_area()
  coord_fixed(ratio = 1) +
  facet_grid(dataset_f ~ .) +
  xlab("X (m)") + ylab("Y (m)") +
  theme_bw(base_size = 16)
ggsave("All_datasets_FS_points_BA.png", width = 8, height = 8, dpi = 600)








### Calculates abundance and occupancy for species at each location
FS_80_03 <- FS_80_02 <- FS_80_01 <- vector("list") # output objects
for(i in 1:ndatasets){
  location <- datasets.out$dataset_code[i]
  quad.sp.df <- geo.df.ras[geo.df.ras$dataset_code==location,]
  quad.sp.df$count <- ifelse(quad.sp.df$ab>0,1,0)
  quad.sp.df$sp <- factor(quad.sp.df$sp)
  occ <- tapply(quad.sp.df$count,quad.sp.df$sp,sum) # species' quadrat occupancy
  t.ab <- tapply(quad.sp.df$ab,quad.sp.df$sp,sum) # species' total abundance
  mn.qab <- tapply(quad.sp.df$ab,quad.sp.df$sp,mean) # species' mean abundance
  mn.qabp <- tapply(quad.sp.df$ab[quad.sp.df$ab>0],quad.sp.df$sp[quad.sp.df$ab>0],mean) # species mean abundance where present
  sum.ba <- tapply(quad.sp.df$ba[quad.sp.df$ab>0],quad.sp.df$sp[quad.sp.df$ab>0],sum) # total basal area in plot across quadrats
  mn.ba <- tapply(quad.sp.df$ba[quad.sp.df$ab>0],quad.sp.df$sp[quad.sp.df$ab>0],mean) # mean total basal area per quadrat where present
  
  aor <- as.data.frame(cbind(occ,t.ab,mn.qab,mn.qabp,sum.ba,mn.ba))
  aor$sp <- names(occ)
  aor$location <- location




######################################################
## Generate dataset for mapping observed values and diversity stats
######################################################
    
rgdat <- as.data.frame(obs.tot_ab) %>%
      mutate(ygrid = row.names(.)) %>%
      gather(key = coln, value = obs.tot_ab, -ygrid) %>%
      separate(coln, into = c("temp", "xgrid"), sep = "V") %>%
      dplyr::select(xgrid, ygrid, obs.tot_ab)
rgdat$dataset = dataset.names[dataset]
rgdat$quad.size = quad.size
rgdat$Latitude = lat
rgdat$obs.richness <- as.data.frame(obs.richness) %>%
      mutate(ygrid = row.names(.)) %>%
      gather(key = coln, value = obs.richness, -ygrid) %>%
      pull(obs.richness)
rgdat$obs.shannon <- as.data.frame(obs.shannon) %>%
      mutate(ygrid = row.names(.)) %>%
      gather(key = coln, value = obs.shannon, -ygrid) %>%
      pull(obs.shannon)
rgdat$obs.simpson <- as.data.frame(obs.simpson) %>%
      mutate(ygrid = row.names(.)) %>%
      gather(key = coln, value = obs.simpson, -ygrid) %>%
      pull(obs.simpson)
rgdat$obs.braycurtis <- as.data.frame(obs.braycurtis) %>%
      mutate(ygrid = row.names(.)) %>%
      gather(key = coln, value = obs.braycurtis, -ygrid) %>%
      pull(obs.braycurtis)
rgdat$obs.tab.sp1 <- as.data.frame(obs.tab.sp1) %>%
      mutate(ygrid = row.names(.)) %>%
      gather(key = coln, value = obs.tab.sp1, -ygrid) %>%
      pull(obs.tab.sp1)
rgdat$obs.tba.sp1 <- as.data.frame(obs.tba.sp1) %>%
      mutate(ygrid = row.names(.)) %>%
      gather(key = coln, value = obs.tba.sp1, -ygrid) %>%
      pull(obs.tba.sp1)
rgdat$obs.mba.sp1 <- as.data.frame(obs.mba.sp1) %>%
      mutate(ygrid = row.names(.)) %>%
      gather(key = coln, value = obs.mba.sp1, -ygrid) %>%
      pull(obs.mba.sp1)

    if(dataset == 1) { raw.gdat <- rgdat }

    if(dataset > 1) { raw.gdat <- rbind(raw.gdat, rgdat) }

    } # end dataset loop

gdat20 <- raw.gdat
gdat10 <- raw.gdat
gdat5 <- raw.gdat

quad.size = 20

# graphs for FS
gdat20 %>%
  mutate(x = as.numeric(xgrid)*quad.size, y = as.numeric(ygrid)*quad.size) %>%
  mutate(dataset_f = factor(dataset, levels=c("AM","BCI","LFDP","TY","HF","WR"))) %>%
  dplyr::select(dataset_f, Latitude, x, y, obs.tab.sp1, obs.tba.sp1, obs.mba.sp1) %>%
  group_by(dataset_f, x, y) %>%
  gather(key = diversity_stat, value = measurement, -x, -y, -dataset_f, -Latitude) %>%
  group_by(diversity_stat) %>%
  mutate(Diversity = (measurement - mean(measurement))/sd(measurement)) %>%
  mutate(diversity_f = factor(diversity_stat, levels=c("obs.tab.sp1","obs.tba.sp1","obs.mba.sp1"))) %>%
  mutate(diversity_f = fct_recode(diversity_f, Abundance = "obs.tab.sp1", `Basal area` = "obs.tba.sp1", `Mean basal area` = "obs.mba.sp1")) %>%
  ggplot(aes(x = x, y = y, fill = Diversity)) +
    geom_raster(stat = "identity") +
    scale_fill_gradientn(colours = terrain.colors(10)) +
    coord_fixed(ratio = 1) +
    facet_grid(dataset_f ~ diversity_f) +
    xlab("X (m)") + ylab("Y (m)") +
    theme_bw(base_size = 16)
ggsave("All_datasets_FS_stats_20m.png", width = 8, height = 8, dpi = 600)


# graphs for community
gdat20 %>%
  mutate(x = as.numeric(xgrid)*quad.size, y = as.numeric(ygrid)*quad.size) %>%
  mutate(dataset_f = factor(dataset, levels=c("AM","BCI","LFDP","TY","HF","WR"))) %>%
  dplyr::select(dataset_f, Latitude, x, y, obs.tot_ab, obs.richness, obs.shannon, obs.simpson, obs.braycurtis) %>%
  group_by(dataset_f, x, y) %>%
  gather(key = diversity_stat, value = measurement, -x, -y, -dataset_f, -Latitude) %>%
  group_by(diversity_stat) %>%
  mutate(Diversity = (measurement - mean(measurement))/sd(measurement)) %>%
  mutate(diversity_f = factor(diversity_stat, levels=c("obs.tot_ab","obs.richness","obs.shannon","obs.simpson","obs.braycurtis"))) %>%
  mutate(diversity_f = fct_recode(diversity_f, Abundance = "obs.tot_ab", Richness = "obs.richness", Shannon = "obs.shannon", Simpson = "obs.simpson", `Bray-Curtis` = "obs.braycurtis")) %>%
  ggplot(aes(x = x, y = y, fill = Diversity)) +
    geom_raster(stat = "identity") +
    scale_fill_gradientn(colours = terrain.colors(10)) +
    coord_fixed(ratio = 1) +
    facet_grid(dataset_f ~ diversity_f) +
    xlab("X (m)") + ylab("Y (m)") +
    theme_bw(base_size = 16)
ggsave("All_datasets_community_stats_20m.png", width = 12, height = 8, dpi = 600)
  





  