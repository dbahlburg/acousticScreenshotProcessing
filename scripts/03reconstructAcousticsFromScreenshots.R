# 19.6.2021
# Dominik Bahlburg
# The aim of this script is to extract relative backscattering intensities
# from screenshots of the acoustic data recorded onboard a commercial krill fishing vessel.
# The script loads a screenshot/photo of a record of acoustic data, matches the 
# colour value of each pixel with the colourscale used for data visualization
# and assigns a relative biomass value (signal strength) to each pixel.
# In the end some post-processing is done in order to remove non-data elements
# of the screenshot and in order to properly scale the x- and y-axes.
#----------------------------------------------------------------------------------#
# Load packages used for the image- and data-processing
library(jpeg)
library(tidyverse)
library(data.table)
library(foreach)
library(parallel)

# list files of acoustic data, attach depth info and extract date and time stamp
acousticDataList <- tibble(filePath = list.files('data/screenshotsSample',pattern = '.jpg', full.names = T),
                           fileName = list.files('data/screenshotsSample',pattern = '.jpg', full.names = F)) 

acousticDataList <- acousticDataList %>% 
  left_join(.,read.csv('data/acousticsDepthInfo.csv', sep = ';')) %>% 
  mutate(timeInfo = regmatches(fileName, gregexpr("[[:digit:]]+", fileName))) %>% 
  filter(!is.na(maxDepth)) %>% 
  rowwise() %>% 
  mutate(
    DateTimeStart = as.POSIXct(paste(paste(unlist(timeInfo)[1:3], collapse = '-'),
                                     paste(unlist(timeInfo)[4:6], collapse = ':'), sep = ' ')),
    DateTimeStop = as.POSIXct(paste(paste(unlist(timeInfo)[7:9], collapse = '-'),
                                    paste(unlist(timeInfo)[10:12], collapse = ':'), sep = ' '))) %>% 
  dplyr::select(filePath, maxDepth, DateTimeStart, DateTimeStop)

#----------------------------------------------------------------------------------#
# create the colour palette used for data visualization by Echoview.
# The colours were extracted from the Echoview-Software-Manual.
colourPalette <- c('#9C8AA8', '#4E4848','#006CFF','#240CAE',
                   '#2AD27E','#078460','#FFFF2A','#FC7831',
                   '#FC5AA8','#FF1836','#B43C30','#962B3C')

# create palette generator which interpolates between the specified colours above
colourPaletteGenerator <- colorRampPalette(colourPalette)

# colour tibble containing the colourscale and dummy-variables (x and y) to 
# plot colourscale. It also contains the rgb-components of each colour
colourTib <- tibble(x = seq(0,1,length.out = 600),
                    y = 1,
                    colour = colourPaletteGenerator(600))  %>% 
  rowwise() %>% 
  mutate(colRGB = list(col2rgb(colour)),
         red = colRGB[[1]],
         green = colRGB[[2]],
         blue = colRGB[[3]])

# plot the colourscale:
ggplot(data = colourTib,
       aes(x = x, y = y, fill = colour)) +
  geom_raster() +
  scale_fill_identity() +
  theme(panel.background = element_rect(fill = NA, colour = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

# start parallelized processing
# temporal resolution
tempResolution <- '1 min'

# vertical resolution in meters
depthResolution <- 0.5

# detect number of cores to paralellize processing (subtract 1 for safety)
nCores <- parallel::detectCores() - 1

#create the cluster
myCluster <- parallel::makeCluster(
  nCores, 
  type = "PSOCK"
)

#register cluster to be used by %dopar%
doParallel::registerDoParallel(cl = myCluster)

# start the parallelized loop (it takes 12 minutes on my machine to process the 70 files)
processedAcoustics <- foreach(i = 1:nrow(acousticDataList), .combine='rbind') %dopar% {
  
  acoustics <- jpeg::readJPEG(acousticDataList$filePath[i])
  
  #convert image into long data format, extract rgb-values for each pixel
  acousticsColours <- tidyr::expand_grid(col = 1:ncol(acoustics),
                                         row = 1:nrow(acoustics)) |> 
    dplyr::mutate(colour = rgb(acoustics[,,1],acoustics[,,2],acoustics[,,3],maxColorValue = 1),
                  red = c(acoustics[,,1]) * 255,
                  green = c(acoustics[,,2] * 255),
                  blue = c(acoustics[,,3]) * 255)
  
  #----------------------------------------------------------------------------------#
  #In this section the colours of the image are matched with the colours from
  #the colourscale. In the following lines, each pixel is matched with all colours
  #from the colourscale. The colours are then matched based on the smallest
  #Euklidean distance in the rgb-colour-space.
  acousticsColoursMatched <- acousticsColours |> 
    dplyr::rowwise() |> 
    dplyr::mutate(closestMatch = which.min(((red - colourTib$red)^2 +
                                              (green - colourTib$green)^2 +
                                              (blue - colourTib$blue)^2)^(1/2)),
                  relativeIntensity = closestMatch/nrow(colourTib))
  
  #Add the exact colour from the legend which has been matched in the previous lines
  #Add red, green and blue components of the colourscale-colours.
  acousticsColoursMatched$legendColour <- colourTib$colour[acousticsColoursMatched$closestMatch]
  acousticsColoursMatched$legendColourRed <- colourTib$red[acousticsColoursMatched$closestMatch]
  acousticsColoursMatched$legendColourGreen <- colourTib$green[acousticsColoursMatched$closestMatch]
  acousticsColoursMatched$legendColourBlue <- colourTib$blue[acousticsColoursMatched$closestMatch]
  
  #----------------------------------------------------------------------------------#
  #Sea floor detection. 
  #Sea floor detection is needed in order to separate the sea floor from the acoustic data.
  #Sea floor detection is also needed to properly scale the y-axis and assign correct depth values.
  
  #add the info about the known sea floor depth at this point (extracted from image info)
  imageDepth <- acousticDataList$maxDepth[i]
  
  #add depth info
  acousticsColoursMatched <- acousticsColoursMatched |> 
    dplyr::ungroup() |> 
    dplyr::mutate(depth = row/max(row) * imageDepth) |> 
    dplyr::mutate(dateTime = acousticDataList$DateTimeStart[i] + col/max(col) * 
                    as.numeric(difftime(acousticDataList$DateTimeStop[i],acousticDataList$DateTimeStart[i], 
                                        units = 'secs')))
  
  #reduce temporal resolution of the data and scale pixels along y-scale accordingly
  acousticsColoursBinned <- acousticsColoursMatched |> 
    dplyr::ungroup() |> 
    dplyr::filter(col < max(col)) |> 
    dplyr::mutate(dateTime = as.POSIXct(cut(dateTime, breaks = tempResolution)),
                  depthBin = cut(depth, breaks = seq(0, 300, by = depthResolution), labels = F),
                  depth = depthBin/max(depthBin) * imageDepth) |> 
    dplyr::group_by(depth, dateTime) |> 
    dplyr::summarise(relativeIntensity = mean(relativeIntensity, na.rm = T)) 
  
  return(acousticsColoursBinned)
}

# close the cluster
parallel::stopCluster(cl = myCluster)

# save reconstructed data as RDS file with contained time period as file name
startTime <- str_replace(as.character(min(acousticDataList$DateTimeStart)), pattern = ' ', replacement = '_')
stopTime <- str_replace(as.character(max(acousticDataList$DateTimeStop)), pattern = ' ', replacement = '_')

write_rds(processedAcoustics, paste('data/reconstructedData/signalReconstructed_',startTime,'_to_',stopTime,'.RDS', sep = ''))

# visualize results
processedAcoustics %>% 
  ggplot(.,aes(x = dateTime, y = -depth, fill = relativeIntensity)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme(legend.position = 'bottom')
