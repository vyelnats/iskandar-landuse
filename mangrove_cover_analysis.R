
#############################################################################################################
                                        ## PRE-PROCESSING ##
############################################################################################################

library(raster)
library(rgdal)
library(lubridate)
library(ggplot2)
library(sf)
library(reshape2)
library(randomForest)
library(caret)
library(dplyr)

# Change raster options to store large rasters in temp files on disk
rasterOptions(maxmemory = 1e6)
data.path <- '/Users/vic/EarthObs/Mangroves/Level2Data' 


# 1) Unpack the archives.
files.1989 <- list.files('/Users/vic/EarthObs/Mangroves/Level2Data/1989', pattern="[1-9].tif$", full.names=T, recursive =F)
files.2018 <- list.files('/Users/vic/EarthObs/Mangroves/Level2Data/2018', pattern="[1-9].tif$", full.names=T, recursive =F)# 2) Stack the image bands of interest.
#files2014 <- list.files('/Users/vic/EarthObs/Mangroves/Level2Data/2014x', pattern="[1-9].tif$", full.names=T, recursive =F)

#stack
stack.1989 <- stack(files.1989) #stack files from Landsat 5 TM 
stack.2018 <- stack(files.2018)

# L5: Band 1 - Blue, Band 2 - Green, Band 3 - Red, Band 4 - NIR, Band 5 - SWIR1, Band 7 - SWIR2
names(stack.1989) <- c('blue', 'green', 'red', 'nIR', 'swIR1', 'swIR2') 
# L8: same but includes ultrablue
names(stack.2018) <- c('ultra blue', 'blue', 'green', 'red', 'nIR', 'swIR1', 'swIR2') 

# 3) Crop the image(s) to the study region extent.
extent.1989 <- extent(stack.1989)
extent.2018 <- extent(stack.2018)
  x1 <- xmin(extent.2018) 
  x2 <- xmax(extent.2018)
  y1 <- ymin(extent.1989)
  y2 <- ymax(extent.1989)

#crop by common extent 
common.extent <- extent(x1, x2, y1, y2) 
#write cropped image to disk 
image.1989 <- crop(stack.1989, common.extent) 
image.2018 <-crop(stack.2018, common.extent)
r89 <- writeRaster(image.1989, filename="rawimage1989", datatype ="FLT4S", format ="GTiff", overwrite=TRUE)
r18 <- writeRaster(image.2018, filename="rawimage1918", datatype ="FLT4S", format ="GTiff", overwrite=TRUE)

#############################################################################################################
                                       ## CLOUD MASKING ##
############################################################################################################
# 4) Mask clouds & cloud shadows in the image(s). 
# From landsat dataproduct level 2, pixel_qa bit index 
fill.index <- 1 # 0 = 1 from R's perspective  
cloud.index <- 6 #5 -> 6 in R 
shadow.index <- 4#3 - > 4 in R 

#Creats a mask that is true for any pixel that is marked as cloud or cloud shadow or fill#
cloud.finder <- function (x) {
 bits <- intToBits(x)
 return (as.logical((bits[cloud.index] | bits[fill.index] | bits[shadow.index]))) }

#load BQA Band and make into raster object apply mask to raster/stack USING CALC
pixelqa.1989 <- raster(list.files('/Users/vic/EarthObs/Mangroves/Level2Data/1989', pattern="pixel", full.names=T, recursive =F))
pixelqa.2018 <- raster(list.files('/Users/vic/EarthObs/Mangroves/Level2Data/2018', pattern="pixel", full.names=T, recursive =F))

#calculate mask values from pixel qa then try cropping mask to Iskendar extent in QGIS 
msk89 <- calc(pixelqa.1989, cloud.finder)
m89 <- writeRaster(msk89, filename="cloudmask89", datatype ="FLT4S", format ="GTiff", overwrite=TRUE)
msk18 <- calc(pixelqa.2018, cloud.finder) 
m18 <- writeRaster(msk18, filename="cloudmask18", datatype ="FLT4S", format ="GTiff", overwrite=TRUE)

#load cropped cloud mask and cropped file
croppedmask98 <- stack('/Users/vic/EarthObs/Mangroves/CroppedCloudLayer.tif')
croppedMap98 <- stack('/Users/vic/EarthObs/Mangroves/CroppedImage98.tif') 
#2018
croppedmask18 <- stack('/Users/vic/EarthObs/Mangroves/CroppedCloudLayer18.tif')
croppedMap18 <- stack('/Users/vic/EarthObs/Mangroves/CroppedImage18.tif') 

#mask BOA data
stk89 <- mask(croppedMap98, croppedmask98, maskvalue=1, updatevalue=NA) 
names(stk89) <- c('blue', 'green', 'red', 'nIR', 'swIR1', 'swIR2') 
s89x <- writeRaster(stk89, filename="image89x.tiff", datatype ="FLT4S", format ="GTiff", overwrite=TRUE)

####REDO
stk18 <- mask(croppedMap18, croppedmask18, maskvalue=1, updatevalue=NA) 
names(stk18) <- c('ultra blue', 'blue', 'green', 'red', 'nIR', 'swIR1', 'swIR2') 
s18x <- writeRaster(stk18, filename="image18x.tiff", datatype ="FLT4S", format ="GTiff", overwrite=TRUE)


###########################################################################################################
          ###   Advanced Pre-Processing, ie calculating spectral/temporal matricies  ###
#################################fi##########################################################################
#read masked raster back in bc had to restart R 
#image89 <- stack('/Users/vic/image89.tif')
#normalize 
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
s18 <- range01(stk18)
s89 <- range01(stk89)

# 5) Calculate NDVI and EVI for Comparitive Overview on Vegetation Health
red.1989 <- stk89[[3]] 
nIR.1989 <- stk89[[4]]
blue.1989 <- stk89[[1]]
green.1989 <- stk89[[2]]

red.2018 <- s18[[4]]  
nIR.2018 <- s18[[5]]
blue.2018 <- s18[[2]]
green.2018 <- s18[[3]]

ndvi <- function(x, y) {
  (x - y) / (x + y)
}

ndvi.1989 <- overlay(nIR.1989, red.1989, fun = ndvi) 
ndvi18 <- (nIR.2018 - red.2018) / (nIR.2018 + red.2018)
ndvi.2018 <- overlay(nIR.2018, red.2018, fun = ndvi) 
plot(ndvi.1989, main = "Normalized Difference Vegetation Index (NDVI) L5 1989") 
plot(ndvi18, main = "Normalized Difference Vegetation Index (NDVI) L8 2018")

# RUN AT THE END BC THESE TAKE A LONG TIME
#evi.1989 <-  2.5 * ((nIR.1989 - red.1989) / (nIR.1989 + 6 * red.1989 - 7.5 * blue.1989 + 10000))
#evi.2018 <-  2.5 * ((nIR.2018 - red.2018) / (nIR.2018 + 6 * red.2018 - 7.5 * blue.2018 + 10000))

##save evi and NVDI
ndvi1989 <- writeRaster(ndvi.1989, filename="1989.ndvi.tiff", datatype ="INT2S", format ="GTiff", overwrite=TRUE)
ndvi2018 <- writeRaster(ndvi.2018, filename="2018.ndvi.tiff", datatype ="INT2S", format ="GTiff", overwrite=TRUE)

##############################################################################################################
          ################ 6) collect training data ##############
#############################################################################################################

#read in ESRI shapefile of Iskandar region
shp <- readOGR(dsn = "/Users/vic/EarthObs/Mangroves/malaysia/", layer = "malaysia")
xp <- stack('/Users/vic/TESTcrop.tif') 

shp <- spTransform(shp, proj4string(xp)) 
plot(xp)
plot(shp, add=T)
xi <- crop(xp, extent(shp))
xp <- mask(xi, shp)  

#crop stacks via shapefile 
# crop the lidar raster using the vector extent
s1989_crop <- crop(image89, extent(shape))

plot(s1989_crop , main = "Cropped Stack 1989") 
# add shapefile on top of the existing raster
plot(image89)
plot(crop_extent, add = TRUE)


##############################################################################################################
################ 7) Random Forest ##############
#############################################################################################################

# load files and create data stack
st89 <- stack('/Users/vic/image89x.tif')
st18 <- stack('/Users/vic/image18x.tif')

#stack the stacks
stack1 <- stack(c(s89, s18))
names(stack1) <- c('blue-89', 'green-89', 'red-89', 'nIR-89', 'swIR1-89', 'swIR2-89', 'ultra blue-18', 'blue-18', 'green-18', 'red-18', 'nIR-18', 'swIR1-18', 'swIR2-18') 
stack1x <- writeRaster(stack1, filename="MangroveStackTest.tif", datatype ="INT2S", format ="GTiff", overwrite=TRUE)

#read my shapefile in with readOGR
mangrove.shapefile <- readOGR(dsn = "/Users/vic/EarthObs/Mangroves", layer = "mangroveclassshapefile")

#extract 
mangrove.points <-as.data.frame(extract(stack1, mangrove.shapefile))
names(mangrove.points) <- c('blue-89', 'green-89', 'red-89', 'nIR-89', 'swIR1-89', 'swIR2-89', 'ultra blue-18', 'blue-18', 'green-18', 'red-18', 'nIR-18', 'swIR1-18', 'swIR2-18') 
mangrove.shapefile@data = data.frame(mangrove.shapefile@data, mangrove.points[match(rownames(mangrove.shapefile@data), rownames(mangrove.points)),])

data1 <- data.frame(class = as.factor(mangrove.shapefile@data$class), mangrove.points)
row.has.na <- apply(data1, 1, function(x){any(is.na(x))})
sum(row.has.na)
data1 <- data1[!row.has.na,]

#############################################################################
# 3) Model Training and classification 
#############################################################################
set.seed(998)
# Split into Train and Validation sets
# Training Set : Validation Set = 70 : 30 (random)
train <- sample(nrow(data1), 0.7*nrow(data1), replace = FALSE)
TrainSet <- data1[train,]
ValidSet <- data1[-train,]
summary(TrainSet)
summary(ValidSet)

# 5) Paramaeterize the model   #na rough fix 
model1 <- randomForest(formula = TrainSet$class~ ., data= TrainSet, mtry=2, ntree = 500, importance = TRUE, na.action=na.roughfix)
model1# Predicting on train set
predTrain <- predict(model1, TrainSet, type = "class", na.action=na.roughfix)
# Checking classification accuracy
table(predTrain, TrainSet$class)  
# Predicting on Validation set
predValid <- predict(model1, ValidSet, type = "class", na.action=na.roughfix)
# Checking classification accuracy
#mean(predValid == ValidSet$class)                    
table(predValid,ValidSet$class)

#to check important variables
importance(model1)
varImpPlot(model1)
partialPlot(model1, data1, class, main= "")

importance = importance(model1)
varImportance = data.frame(Variables = row.names(importance),
                           Importance =round(importance[, "MeanDecreaseGini"],2))

rankImportance=varImportance%>%mutate(Rank=paste('#',dense_rank(desc(Importance))))
  
ggplot(rankImportance,aes(x=reorder(Variables,Importance),
                          y=Importance,fill=Importance))+ 
  geom_bar(stat='identity') + 
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'white') +
  labs(x = 'Variables') +
  coord_flip() + 
  theme_classic() + ggtitle("Importance of variables in the model") +
  theme(plot.title = element_text(hjust = 0.5))   



   
ggplot(rankImportance,aes(x=reorder(Variables,Importance), y=Importance,fill=Importance))+ 
                                                    geom_bar(stat='identity') + 
                                                    geom_text(aes(x = Variables, y = 0.5, label = Rank),
                                                              hjust=0, vjust=0.55, size = 4, colour = 'white') +
                                                    
                                                    labs(x = 'Variables') 
                                                    coord_flip() + 
                                                    theme_classic() + ggtitle("Importance of variables in the model") +
  theme(plot.title = element_text(hjust = 0.5))   
                                                  
                                                  
                                                  

#Write the resulting map to disk in GTiff
prediction <- predict(model1, data1, type = 'class', filename= "prediction", na.rm = TRUE, format = 'GTiff', progress = 'text', index = 1, na.rm=TRUE, overwrite=TRUE)
prediction.valid <- predict(model1, TrainSet, type = 'class', filename= "prediction.trainset", na.rm = TRUE, OOB = TRUE, format = 'GTiff', progress = 'text', index = 1, na.rm=TRUE, overwrite=TRUE)
prediction <- predict(stack1, model1)

### K fold Validation ### 

# fix the parameters of the algorithm
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",   #repeated cross validation
  number = 10, # number of folds in k-fold 
  ## repeated ten times
  repeats = 10)

set.seed(825)
gbmFit1 <- train(class ~ ., data = data1, 
                 method = "gbm", #gradient boosting machine
                 trControl = fitControl,
                 verbose = FALSE)
gbmFit1

gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9), 
                        n.trees = (1:30)*50, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20)

nrow(gbmGrid)

set.seed(825)
gbmFit2 <- train(class ~ ., data = data1, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneGrid = gbmGrid)
gbmFit2



