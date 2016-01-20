# Report_Code_SDM
# Maxent code for results in report III



## This code is for Mammals. To chaange code to other taxa, substitute "Mammals" with taxa name

## Mammal: spp of interest

#setwd("//Users/kpradhan/Desktop/02_masters stuff/winter break work/sdm/") #For working on mac. Comment when working on lab pc
# Load the following packages

library(dismo)
library(maptools)
library(rgeos)
library(grid)
library(ClustOfVar)
library(rgdal)
library(sm)
library(foreach)

# STEP: Read species of interest data ####
sppIntOrg <- read.csv("sppOfInterest.csv")
sppIntOrg2 <- read.csv("sppOfInterest2.csv")
sppIntOrg3 <- read.csv("sppOfInterest3.csv")
sppIntOrg4 <- read.csv("panthera_pardus.csv")
# STEP: Load info file ####
sppInfo <- read.csv("sppInfo.csv")
sppInt <- sppIntOrg[,c(10,16,17,22)]
colnames(sppInt) <- c("species", "lat", "long", "recordtype")
colnames(sppIntOrg2) <- c("species", "lat", "long", "recordtype")
colnames(sppIntOrg3) <- c("species", "lat", "long", "recordtype")
colnames(sppIntOrg4) <- c("species", "lat", "long", "recordtype")
sppInt <- rbind(sppInt, sppIntOrg2, sppIntOrg3, sppIntOrg4)
mammal <- sppInfo$Species[sppInfo$Type == "mammal"]
idx <- sppInt$species %in% mammal
sppInt <- sppInt[idx,]
# Get all mammal species names
usp <- unique(sppInt$species)
usp

# STEP: Get climate rasters ####
climAll <- getData("worldclim", var= "bio", res=10)

# STEP: Get region of interest and buffer around it ####
regionBuff10 <- readOGR(dsn="Data", "regionFinal_buffered_10")
plot(regionBuff10)
ext10 <- extent(regionBuff10)

# Crop climate data to extent
climExt10 <- crop(climAll, ext10)
rm(climAll)

# STEP: Remove global species ####
# 1. Remove Felis silvestris- wide rangeing species
sppInt <- sppInt[sppInt$species != "Felis silvestris",]

# STEP: Now check for number of non-duplicated observations in each species in each cell 
# to remove species with too few observations
# 1. remove duplicated observations
dups <- duplicated(sppInt[,c("lat", "long", "species")])
# remove duplicates
sppInt <- sppInt[!dups,]
# 2. assign site id to remove duplicates within sites
sppInt$site <- cellFromXY(climExt10[[1]],sppInt[,c("long", "lat")])
dups <- duplicated(sppInt[,c("site", "species")])
sppIntSite <- sppInt[!dups,]
usp <- unique(sppIntSite$species)
usp

# 3. crop occurences to region
# On a species by species basis, retain only one observation per grid cell
# select occurences in our desired extent
occAllSp <- sppIntSite
coordinates(occAllSp) <- occAllSp[, c("long", "lat")]
proj4string(occAllSp) <- proj4string(regionBuff10)
region2 <- as(regionBuff10, "SpatialPolygons")
selectOcc <- over(occAllSp, region2)
sppIntSite$select <- selectOcc
pres10 <- subset(sppIntSite, sppIntSite$select != "NA")
length(unique(pres10$species))
pres10 <- pres10[,-4]
usp <- unique(pres10$species)

# 4. Check number of occurences per species

for(i in usp){
  print(i)
  df <- subset(pres10, pres10$species==i)
  print(length(df$species))
}

# Remove species with < 10 occurences
# Oryx leucoryx  only has 4 observations
pres10 <- pres10[pres10$species != "Oryx leucoryx",]
# Gazella subgutturosa only has 4 observations
pres10 <- pres10[pres10$species != "Gazella subgutturosa",]
# Vulpes rueppellii only has 5 observations
pres10 <- pres10[pres10$species != "Vulpes rueppellii",]
# Get final list of species to work with
usp <- unique(pres10$species)

# check predictor variables for selection
## I checked the predictor variables and looking at cluster dendrograms for each species selected variables. Here I am making  a list of the variables for each species to work with in the for loop
CI <- c(1,4,13,19)
GG <- c(1,4,11,12,18)
CC <- c(1,4,11,12,18)  
MC <- c(1,2,4,10,12,15,19) 
VC <- c(1,4,10,12)
PP <- c(1,2,4,12,18,19)

modList <- list()
varList <- list(CI,GG,CC,MC,VC,PP)
j <- 1
for(i in usp){
  # subset data
  spPres10 <- subset(pres10, pres10$species==i) 
  spPres10 <- gridSample(spPres10[,c(2,3)], climExt10, n=1)
  spPres10 <- as.data.frame(spPres10)
  l <- length(spPres10$lat)
  
  #######TrendSurfaceAnalysis
  tsaData <- surf.ls(np = 3, x = spPres10$lat, )
  
  
  # create kernel density estimation bias raster
  # create raster mask
  maskRast <- climExt10[[1]]>-132
  
  # get coordinates of presences
  pres10XY <- spPres10
  coordinates(pres10XY) <- pres10XY[,c("long", "lat")]
  bias <- cellFromXY(maskRast, pres10XY)
  cells <- unique(sort(bias))
  kernelXY <- xyFromCell(maskRast, cells)
  colnames(kernelXY) <- c("x", "y")
  samps <- as.numeric(table(bias))
  
  # code to make KDE raster
  #########figure out what ngrid to use and what it means########
  KDEsur <- sm.density(kernelXY, weights=samps, display="none", ylim=c(-5.833333,59.33333), xlim=c(-25.33333,84.83333), nbins=0, ngrids = 661)
  KDErast=SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
  KDErast = SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, length(KDEsur$estimate))))
  KDErast <- raster(KDErast)
  KDErast <- resample(KDErast, maskRast)
  KDErast <- KDErast*maskRast
  KDEpts <- rasterToPoints(KDErast)
  
  # make sets of training and testing data
  # split data for training and testing
  # set 30% testing and 70% training
  splitPercent <- 0.30 
  
  # function to partition data into t/t
  fold <- function(splitPercent){ 
    k = round(1/splitPercent, 0)
    fold <- kfold(spPres10, k=k)
  }
  
  sets = 10 # number of sets
  folds <- replicate(sets, fold(splitPercent))
  spPres10 <- spPres10[,c(2,1)] 
  #looping to make training and testing datasets
  spTrain <- list()
  spTest <- list()
  for(h in 1:sets){
    spTrain[[h]] <- spPres10[folds[,h]!=1,]
    spTest[[h]] <- spPres10[folds[,h]==1,]
  }
  
  # Run Maxent with KDE bias
  system.file("java", package="dismo")
  for (k in 1:sets){
    spMaxMods_allF_bias <- foreach(k=1:sets, .verbose=T, .packages=c("dismo", "rJava", "rgdal")) %dopar% maxent(climExt10[[varList[[j]]]], spTrain[[k]], a=KDEpts[sample(seq(1:nrow(KDEpts)), size=10000, replace=T, prob=KDEpts[,"layer"]),1:2], path="maxent", args=c("responsecurves=true"))
    
    spMaxStack_allF_bias <- predict(spMaxMods_allF_bias[[1]], climExt10[[varList[[j]]]])
    spMaxStackRAW_allF_bias <- predict(spMaxMods_allF_bias[[1]], climExt10[[varList[[j]]]], args='outputformat=raw')
    
    for(l in 2:sets){
      mod_allF_bias <- predict(spMaxMods_allF_bias[[l]], climExt10[[varList[[j]]]])
      modRAW_allF_bias <- predict(spMaxMods_allF_bias[[l]], climExt10[[varList[[j]]]], args='outputformat=raw')
      spMaxStack_allF_bias <- stack(spMaxStack_allF_bias, mod_allF_bias)
      spMaxStackRAW_allF_bias <- stack(spMaxStackRAW_allF_bias, modRAW_allF_bias)
    }
  }
  save.image(file=paste(i, "_withBiasCorrection.RData", sep=""))
  j <- j+ 1
}  


## To plot current predictions saved in the above RData files ####

MammalFiles <- list.files(pattern="_withBiasCorrection.RData", recursive=T)
for(m in 1:length(MammalFiles)){
  load(MammalFiles[m])
  fileName <- strsplit(MammalFiles[m], split="_")[[1]][1]
  writeRaster(sum(spMaxStack_allF_bias)/dim(spMaxStack_allF_bias)[3], paste(fileName, "_maxent_biasCorrected.tiff", sep=""), type="GTiff", overwrite=T)
  writeRaster(calc(spMaxStack_allF_bias, sd), paste(fileName, "_maxent_sd_biasCorrected.tiff", sep=""), type="GTiff", overwrite=T)
}

## To evaluate model prediction ####
MammalFiles <- list.files(pattern="_withBiasCorrection.RData")
maxentEval_bias <- list()  
for(n in 1:length(MammalFiles)){
  load(MammalFiles[n])
  AUCmod <- corrMod <- AVI <- thresh95 <- predArea95 <- threshMPA <- MPA <- meanProb <- meanBG <- BoyceI <- NULL
  
  # predicted probability at random background points
  probBG <- extract(spMaxStack_allF_bias, randomPoints(climExt10, 10000))
  
  for(nn in 1:dim(spMaxStack_allF_bias)[3]){    
    # predicted probability at test points
    probTest <- as.numeric(na.omit(extract(spMaxStack_allF_bias[[nn]], spTest[[nn]])))
    evalDismo <- evaluate(p=probTest, a=probBG[,nn])
    AUCmod[[nn]] <- evalDismo@auc
    corrMod[[nn]] <- evalDismo@cor
    meanProb[[nn]] <- mean(probTest)
    meanBG[[nn]] <- mean(probBG[,nn])
    AVI[[nn]] <- sum(probTest>0.5)/length(probTest)
    thresh95[nn] <- sort(probTest, decreasing=T)[round(length(probTest)*0.95,0)]
    #predArea95[[nn]] <- count(spMaxStack_allF_bias[[nn]]>=thresh95[nn],1)/area
  }
  
  maxentEval_bias[[n]] <- rbind(thresh95, threshMPA, MPA, AVI, AUCmod,
                                meanProb, meanBG, corrMod)
  
  cat(round(n/length(MammalFiles)*100,2), "%", sep="", fill=T)
}
dput(maxentEval_bias, "MammalEval_bias.txt")

## To plot future predictions ####

# Load future climate Rasters

#indir <- "//mf-drobo/projects/ClimateData/future/2_5min/"
year <- c("2030", "2050", "2070")
scenario <- "/IPCC_AR5/ncar_ccsm4/"
rcp <- c("rcp4_5", "rcp6_0", "rcp8_5")
bioList <- c("/bio_1.asc","/bio_2.asc","/bio_3.asc","/bio_4.asc","/bio_5.asc","/bio_6.asc","/bio_7.asc","/bio_8.asc","/bio_9.asc","/bio_10.asc","/bio_11.asc","/bio_12.asc","/bio_13.asc","/bio_14.asc","/bio_15.asc","/bio_16.asc","/bio_17.asc","/bio_18.asc","/bio_19.asc")

MammalFiles <- list.files(pattern="_withBiasCorrection.RData", recursive=T)

# for mac. comment on lab pc
indir <- "//Users/kpradhan/Desktop/02_masters stuff/winter break work/"
year <- "2030"
rcp <- "rcp4_5"
bioList <- c("/bio_1.asc","/bio_2.asc","/bio_3.asc","/bio_4.asc","/bio_5.asc","/bio_6.asc","/bio_7.asc","/bio_8.asc","/bio_9.asc","/bio_10.asc","/bio_11.asc","/bio_12.asc","/bio_13.asc","/bio_14.asc","/bio_15.asc","/bio_16.asc","/bio_17.asc","/bio_18.asc","/bio_19.asc")


for(y in year){
  for(r in rcp){
    path <- paste(indir, r, bioList,".txt", sep="")
    glo1 <- raster(path[1])
    glo2 <- raster(path[2])
    glo3 <- raster(path[3])
    glo4 <- raster(path[4])
    glo5 <- raster(path[5])
    glo6 <- raster(path[6])
    glo7 <- raster(path[7])
    glo8 <- raster(path[8])
    glo9 <- raster(path[9])
    glo10 <- raster(path[10])
    glo11 <- raster(path[11])
    glo12 <- raster(path[12])
    glo13 <- raster(path[13])
    glo14 <- raster(path[14])
    glo15 <- raster(path[15])
    glo16 <- raster(path[16])
    glo17 <- raster(path[17])
    glo18 <- raster(path[18])
    glo19 <- raster(path[19])
    climAllFut <- stack(glo1, glo2, glo3, glo4, glo5,glo6, glo7, glo8, glo9, glo10,glo11, glo12, glo13, glo14, glo15,glo16, glo17, glo18, glo19)
    rm(glo1)
    rm(glo2)
    rm(glo3)
    rm(glo4)
    rm(glo5)
    rm(glo6)
    rm(glo7)
    rm(glo8)
    rm(glo9)
    rm(glo10)
    rm(glo11)
    rm(glo12)
    rm(glo13)
    rm(glo14)
    rm(glo15)
    rm(glo16)
    rm(glo17)
    rm(glo18)
    rm(glo19)
    
    climExtAllFut <- crop(climAllFut, ext10)
    names(climExtAllFut) <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18", "bio19")
    climFut10 <- climExtAllFut
    j <- 1
    for(m in 1:length(MammalFiles)){
      load(MammalFiles[m])
      fileName <- strsplit(MammalFiles[m], split="_")[[1]][1]
      
      spMaxStack_bias_Fut <- predict(spMaxMods_allF_bias[[1]], climFut10[[varList[[j]]]])
      spMaxStackRAW_bias_Fut <- predict(spMaxMods_allF_bias[[1]], climFut10[[varList[[j]]]], args='outputformat=raw')
      
      for(jj in 2:sets){
        mod_bias_Fut <- predict(spMaxMods_allF_bias[[jj]], climFut10[[varList[[j]]]])
        modRAW_bias_Fut <- predict(spMaxMods_allF_bias[[jj]], climFut10[[varList[[j]]]], args='outputformat=raw')
        spMaxStack_bias_Fut <- stack(spMaxStack_bias_Fut, mod_bias_Fut)
        spMaxStackRAW_bias_Fut <- stack(spMaxStackRAW_bias_Fut, modRAW_bias_Fut)
      }
      writeRaster(sum(spMaxStack_bias_Fut)/dim(spMaxStack_bias_Fut)[3], paste(fileName, "_", y, "_", r, "_maxent_biasCorrected.tiff", sep=""), type="GTiff", overwrite=T)
      writeRaster(calc(spMaxStack_bias_Fut, sd), paste(fileName, "_", y, "_", r, "_maxent_sd_biasCorrected.tiff", sep=""), type="GTiff", overwrite=T)
    }
  }
}


##################### here is where you stopped ########### 









# future predictions
# load future climate
indir <- "//mf-drobo/projects/ClimateData/future/2_5min/"
year <- c("2030", "2050", "2070")
scenario <- "/IPCC_AR5/ncar_ccsm4/"
rcp <- c("rcp4_5", "rcp6_0", "rcp8_5")
bioList <- c("/bio_1.asc","/bio_2.asc","/bio_3.asc","/bio_4.asc","/bio_5.asc","/bio_6.asc","/bio_7.asc","/bio_8.asc","/bio_9.asc","/bio_10.asc","/bio_11.asc","/bio_12.asc","/bio_13.asc","/bio_14.asc","/bio_15.asc","/bio_16.asc","/bio_17.asc","/bio_18.asc","/bio_19.asc")

for(y in year){
  for(r in rcp){
    path <- paste(indir, y, scenario, r, bioList, sep="")
    glo1 <- raster(path[1])
    glo2 <- raster(path[2])
    glo3 <- raster(path[3])
    glo4 <- raster(path[4])
    glo5 <- raster(path[5])
    glo6 <- raster(path[6])
    glo7 <- raster(path[7])
    glo8 <- raster(path[8])
    glo9 <- raster(path[9])
    glo10 <- raster(path[10])
    glo11 <- raster(path[11])
    glo12 <- raster(path[12])
    glo13 <- raster(path[13])
    glo14 <- raster(path[14])
    glo15 <- raster(path[15])
    glo16 <- raster(path[16])
    glo17 <- raster(path[17])
    glo18 <- raster(path[18])
    glo19 <- raster(path[19])
    climAllFut <- stack(glo1, glo2, glo3, glo4, glo5,glo6, glo7, glo8, glo9, glo10,glo11, glo12, glo13, glo14, glo15,glo16, glo17, glo18, glo19)
    rm(glo1)
    rm(glo2)
    rm(glo3)
    rm(glo4)
    rm(glo5)
    rm(glo6)
    rm(glo7)
    rm(glo8)
    rm(glo9)
    rm(glo10)
    rm(glo11)
    rm(glo12)
    rm(glo13)
    rm(glo14)
    rm(glo15)
    rm(glo16)
    rm(glo17)
    rm(glo18)
    rm(glo19)
    
    climExtAllFut <- crop(climAllFut, ext10)
    names(climExtAllFut) <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18", "bio19")
    climExtAllFut[[1]]
    mask <- climExtAllFut[[1]]>-186
    climExtFut <- mask*climExtAllFut
    j <- 1
    for (i in usp){ 
      climExtFuture <- climExtAllFut[[varList[[j]]]]
      modelSP <- modList[[j]]
      predFut <- predict(climExtFuture, modelSP, args=c("outputformat=logistic"))
      plot(predFut)
      writeRaster(predFut, filename=paste(i,"_",y, "_", r,"_Maxent.tif", sep=""))
      j <- j+1
    }
    
  }
}



predFut <- predict(climExtFuture, modMX, args=c("outputformat=logistic"))
plot(predFut)
writeRaster(predFut, filename=paste(i,"_",y, "_", r,"_Maxent.tif", sep=""))

