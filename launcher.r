source("<REAL_PATH>\\R-function.r")

#SETTINGS
#Change path 
path = "<REAL_PATH>"
outpath = paste(path, "outputs\\", sep="")	#output path [default = "outputs"]

level = 3			#level of the Corine Land Cover (europe only) (1,2,3) to be considered in the analysis [default = 3]
threshold = 0.75	#threshold below which the combinations of experts are excluded based on the concordance between ahp weights [default = 0.75]
npoints = 70		#number of sampling sites (traps)
mindist = 8000		#minimum distance (meters) between traps [default = 4000]
radius = 2000		#radius (meters) of the buffer, for the Land Cover analisys	[default = 2000]
imax = 1000			#number of iteration

#studyarea <- readOGR(paste(path,"data\\shp\\study_area.shp",sep=""), "study_area", verbose=TRUE)	#study area boundaries
studyarea <- readOGR(paste(path,"<REAL_PATH><STUDY_AREA_NAME>.shp",sep=""), "<STUDY_AREA_NAME>", verbose=TRUE)	#study area boundaries
csvweights <- read.csv(file=paste(outpath,"<CLC_WEIGHTS>.csv",sep=""), head=TRUE, sep=";")				#csv weights 
#rstclc <- raster(paste(path, "data\\rst\\g100_06_as1.tif", sep=""))
rstclc <- raster(paste(path, "<REAL_PATH><CLC_RASTER>.tif", sep=""))
csvclass <- read.csv(file=paste(path, "<REAL_PATH><CLC_CLASSIFICATION>.csv",sep=""), head=TRUE, sep=";")

#GET DATA
dataframe <- ahp.get.data(path=path, level=level, view=FALSE)

#GEOMETRIC MEAN CALCULATION WITH BEST EXPERTS COMBINATION (if export = TRUE select the export file location)
ahp.val.gm(path=outpath, dataframe=dataframe, threshold=threshold, view=FALSE)

#GEO ANALYSIS
outputs <- ahp.geo.analysis(csvWeights=csvweights, rstClc=rstclc, csvClass=csvclass, studyArea=studyarea,  npoints=npoints, mindist=mindist, radius=radius, imax=imax, showtime=TRUE, view=TRUE)

#VIEW DATA
#prjcrs="+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +units=m +no_defs"
#imgpoints <- SpatialPoints(df.points[c("x","y")], proj4string=CRS(prjcrs), bbox=shpareastudio@bbox)
#plot(rstclc)
#plot(shpareastudio, add=TRUE)
#plot(imgpoints, add=TRUE, col="black", pch=1)
plot(outputs[[1]])
plot(studyarea, add=TRUE)
plot(outputs[[2]],add=TRUE, col="black", pch=1)

#CHECK THE MODEL
modeltest <- raster(paste(path, "<REAL_PATH><CLC_RECLASS>.tif", sep=""))
check.outputs <- ahp.geo.check(rstModelTest=modeltest, list.data=outputs[[3]], studyArea=studyarea)



#change NA to 0
modeltest <- setValues(modeltest, ifelse(is.na(values(modeltest)),0,values(modeltest)))
spdataall <- data.frame(id = numeric(length=0), val = numeric(length=0), checkval = numeric(length=0))

#ciclo su list
for(i in 1:length(outputs[[3]])){
	spdata = SpatialPoints(outputs[[3]][[i]][2:3], proj4string=CRS(prjcrs), bbox = studyarea@bbox)
	spdatatmp <- data.frame(id = i, val = sum(outputs[[3]][[i]]$value), checkval = sum(extract(modeltest, spdata)))
	spdataall <- rbind(spdataall, spdatatmp)
}


#shp <- readOGR(paste(path,"out_20170315\\sites.shp",sep=""), "sites", verbose=TRUE)
#shp@proj4string@projargs = "+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +units=m +no_defs"
#p4s <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
#shpwgs84 <- spTransform(shp, CRS= p4s)
#writeOGR(obj=shpwgs84, dsn=paste(outpath,"kml\\sites.kml",sep=""), layer="sites", driver="KML")