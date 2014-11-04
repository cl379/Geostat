# welcome
cat('\n\n *** auto-Kriging loaded *** \n\n')
# =======================================================================================
library(sp)
library(rgdal)
library(raster)
library(gstat)
# =======================================================================================
# enter dataset
dataset <- read.table('dataset.csv', sep = ';', dec = '.', header = T)
# subsample
example <- data.frame(PHI = dataset$PHI, X = dataset$X, Y = dataset$Y)
# coerces to SpatialPointsDataFrame class
coordinates(example) <- ~ X + Y
# =======================================================================================
# function for generating a SpatialPixels class with a SpatialPointsDataFrame object
# parameter nCells is the number of cells in the X axis 
makeGrid <- function(spatialPoints, nCells) {
	# establish the x and y ranges and modify with margins parameter
	xRange <- range(spatialPoints@coords[, 1])
	yRange <- range(spatialPoints@coords[, 2])
	# establish cell size
	cellSize <- (xRange[2] - xRange[1]) / nCells
	# create the x and y points
	xPoints <- seq(from = xRange[1], to = xRange[2], by = cellSize)
	yPoints <- seq(from = yRange[1], to = yRange[2], by = cellSize)
	# create the 'SpatialPixels' class
	grd <- expand.grid(X = xPoints, Y = yPoints)
	coordinates(grd) <- ~ X + Y
	gridded(grd) <- TRUE
	# return the grid ('SpatialPixels' class)
	return(grd)
}
# =======================================================================================
# function for generating an Inverse Distance Weighting interpolation with a SpatialPointsDataFrame
# (only with a 1 variable)
invDistWei <- function(spatialPoints, grd, makeOuputs=TRUE) {
	# generate the 'SpatialPixelsDataFrame' class
	idw <- idw(formula = spatialPoints@data[, 1] ~ 1, 
			locations = spatialPoints, 
			newdata = grd)
	# create the 'raster' class
	rasterIdw <- raster(idw)
	# optional plot and export to GTiff (using raster and rgdal packages)
	if (makeOuputs) {
		png('rasterIDW.png', width = 600, height = 600, units = "px")
		plot(rasterIdw, 
				ylab = 'Y coord', 
				xlab = 'X coord', 
				main = 'Inverse Distance Weighting')
		dev.off()
		writeRaster(rasterIdw, 
				filename = 'rasterIdw.tif', 
				format = 'GTiff', 
				overwrite = TRUE)
	}
	# return the interpolation ('SpatialPixelsDataFrame' class)
	return(idw)
}
# =======================================================================================
# function for generating an Ordinary Kriging interpolation with a SpatialPointsDataFrame class
# (only with a 1 variable)
ordinaryKriging <- function(spatialPoints, grd, makeOuputs=TRUE){
	# construction of variogram
	semivariog <- variogram(spatialPoints@data[, 1] ~ 1, 
			locations = spatialPoints, 
			data = spatialPoints)
	# estimate the parameters (don't round)
	rang <- semivariog$dist[length(semivariog$dist)]
	sill <- semivariog$gamma[length(semivariog$gamma)]
	nugg <- semivariog$gamma[1]
	# optional ouputs
	if (makeOuputs) {
		# plot semivariogram (through a function because a gstat's code problem)
		png('variogram.png', width = 600, height = 645, units = "px")
		plot(plotVariog(semivariog, 
					paste('SEMIVARIOGRAM\n\nrang:',rang,'\nsill:',sill,'\nnugg:',nugg)))
		dev.off()
		# generate the models of semivariogram, plot it and export to GTiff
		for (i in c('Exp', 'Sph', 'Gau', 'Lin')){
			makeModel(semivariog, sill, nugg, rang, i, spatialPoints, grd)
		}
	}
	# Return the semivariogram
	return(semivariog)
}
# =======================================================================================
# auxiliar function for generating the Ordinary Kriging models and ouputs
makeModel <- function(semivariog, sill, nugg, rang, model, spatialPoints, grd){
	# make model and fit it
	modelVariog <- vgm(psill = sill, model = model, nugget = nugg, range = rang)
	fitVariog <- fit.variogram(semivariog, modelVariog)
	# plot semivariogram (through a function because a gstat's code problem)
	png(paste(model, '_fit.png', sep = ''), width = 600, height = 600, units = "px")
	plot(plotVariogFit(semivariog, fitVariog, model))
	dev.off()
	# generate the kriging surface prediction
	OK <- krige(formula = spatialPoints@data[, 1] ~ 1, 
			locations = spatialPoints, 
			newdata = grd, 
			model = modelVariog)
	# convert to raster class
	rasterOK <- raster(OK)
	# plot the image (using raster package)
	png(paste(model, '.png', sep = ''), width = 600, height = 600, units = "px")
	plot(rasterOK,  main = model)
	dev.off()
	# export to GTiff
	writeRaster(rasterOK, 
			filename = paste(model, '.tif', sep = ''), 
			format = 'GTiff', 
			overwrite = TRUE)
}
# =======================================================================================
# functions for plotting the variograms and fit semivariograms (see the 76 and 93 lines)
plotVariog <- function(semivariog, main){
	return(plot(semivariog, main = main))
}
plotVariogFit <- function(semivariog, fitVariog, model){
	return(plot(semivariog, fitVariog, main = model))
}
# =======================================================================================
