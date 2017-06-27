#####################################################################################################
#																									#
#	FUNCTIONS PER LA VERIFICA DELLA CONCORDANZA PESI ESPERTI AHP									#
#	PROGETTO	= 	AHP VETTORI																		#
#	AUTHORS		=	PAOLO MULATTI, MATTEO MAZZUCATO													#
#	VERSIONE	=	1.0																				#
#																									#	
#####################################################################################################

# INPUT and VARs
# path = expert elicitation path files
# level = level (1,2,3) to be considered in the analysis 
# threshold = threshold below which the combinations of experts are excluded based on the concordance between ahp weights. Default value = 0.75(75%).

# npoints =	number of random points generated
# prjcrs = prject reference system. Defeault value(UTM32) = "+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +units=m +no_defs"
# imax = maximum number of iterations
# rstclc = raster("data\\rst\\clcl06clip.tif")

require(gtools)
require(irr)
require(Kendall)

#Get data from a specific path
ahp.get.data = function(path, level, view=FALSE){
	path <- paste(path,"data\\ahp\\",sep="")
	csv.list <- list.files(path = path, pattern="*.csv")
	csv.names <- strsplit(csv.list,"_")
	csv.names <- lapply(csv.names,function(x) x[1])
	ncsv <- do.call("c",csv.names)
	data.list <- list()
	d <- length(csv.list)
	for (i in 1:d) {
		data.list[[i]] <- read.csv(paste(path,csv.list[i],sep=""), header=TRUE, sep=";")
	}
	data.list = lapply(data.list,function(x) x[order(x$elemcod),])
	dlist = lapply(data.list,function(x) subset(x,x$level==level))
	dmat = data.frame(CLASS=dlist[[1]]$elemcod, do.call("cbind",lapply(dlist,function(x) x["ahp_res"])))
	names(dmat)[-1]=ncsv
	if(view == TRUE){	
		ahp.plot.data(dmat)
	}
	return(dmat)
}

#Plot ahp data with histograms
ahp.plot.data = function(dataframe){
	windows(height=10,width=20)
	par(mfrow=c((length(names(dataframe))-1)/5,5))
	for(i in 2:ncol(dataframe))
	{
		barplot(dataframe[,i], main=names(dataframe)[i], ylim=c(0,round(max(apply(dataframe[2:length(dataframe)], 2, function(x) max(x))),1)+0.1))
	}
}



#Best combination from AHP experts
ahp.combo.best = function(dataframe, threshold=0.75){
	n.df = names(dataframe)[-1]
	combo.w = list()
	w.list = list()
	
	for(i in 2:length(n.df))
	{
	 w.list[[i]] = list()
	 combo.nm = combinations(length(n.df),i,n.df)
 	 w = integer(length=0)
	 id = vector(length=0)
	 l = integer(length=0)
		for(j in 1:nrow(combo.nm))
		{
		 id[j] = paste(combo.nm[j,],collapse="-")
		 l = length(combo.nm[j,])
		 sub.data = dataframe[,names(dataframe)%in%combo.nm[j,]]
		 w.list[[i]][[j]] = kendall(sub.data,correct=T)	
		 w[j] = w.list[[i]][[j]]$value
		}
	 combo.w[[i]] = data.frame(id=id,w=w,l=l)
	}
	df.w = do.call("rbind",combo.w)
	
	u1=df.w
	u1$ok=u1$w*ifelse(u1$w>=threshold,1,NA)
	u2=split(u1,u1$l)
	u3=lapply(u2,function(x) x[which.max(x$ok),])
	u4=do.call("rbind",u3)
	u5=u4[which.max(u4$l),]
	u.id=unlist(strsplit(as.character(u5$id),"-"))
	df.w.s=dataframe[,u.id]
	
	return(df.w.s)
}

#Geometric Mean
gm_mean = function(x,na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#Export csv file for analysis
ahp.val.gm = function(path, dataframe, threshold=0.75, view=FALSE, type="gm"){
	df.w.s = ahp.combo.best(dataframe, threshold)
	dm <- as.matrix(df.w.s)
	
	gmean = switch(type,
              mean = apply(dm,1,mean),
              median = apply(dm,1,median),
              gm = apply(dm,1,gm_mean) 
			)
			
	gm.df=data.frame(CLASSLVL3=substr(dataframe$CLASS,4,length(dataframe$CLASS)),PESO=gmean/sum(gmean))
	if(view == TRUE){	
		gm.df	
	}
	write.table(gm.df,paste(path,"clcweights.csv",sep=""),sep=";",row.names=FALSE,col.names=TRUE)
	return(gm.df)
}

#SAW
require("raster")
require("rgeos")
require("spatstat")
require("rgdal")


ahp.geo.analysis = function(csvWeights, rstClc, csvClass, studyArea, npoints, mindist, radius=2000, prjcrs="+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +units=m +no_defs", imax, showtime=FALSE, view=FALSE){
	prt=proc.time()
	
	print("STARTING PROCESS...")
	
	i = 1
	nodataval = 128
	res.list <- list()
	dfClassLVL3 <- merge(csvClass, csvWeights, by.x = "CLC_CODE", by.y="CLASSLVL3",all.x=T)
	mtxPesiLVL3 <- cbind(csvClass$GRID_CODE, csvClass$GRID_CODE+1, csvClass$CLC_CODE)
	mtxPesiLVL3 <- rbind(mtxPesiLVL3,c(nodataval,nodataval+1,999))
	mtxPesiLVL3 <- cbind(dfClassLVL3$GRID_CODE, dfClassLVL3$GRID_CODE+1, dfClassLVL3$PESO)
	mtxPesiLVL3 <- rbind(mtxPesiLVL3,c(nodataval,nodataval+1,NA))
	rstclclvl3 <- reclassify(rstClc, mtxPesiLVL3, include.lowest=TRUE, right=FALSE)
	
	polyas=lapply(studyArea@polygons[[1]]@Polygons,function(x) x@coords[nrow(x@coords):2,])
	shpas <- owin(xrange=studyArea@bbox[1,], yrange=studyArea@bbox[2,], poly=polyas)
	
	rstbuf <- raster(nrows=rstclclvl3@nrows, ncols=rstclclvl3@ncols, xmn=rstclclvl3@extent@xmin, xmx=rstclclvl3@extent@xmax, ymn=rstclclvl3@extent@ymin, ymx=rstclclvl3@extent@ymax, crs=prjcrs)
	while(i<=imax){
		rndpoints <- rSSI(r=mindist, n=npoints, win=shpas)	
		if(rndpoints$n == npoints){
			dfrandompoints <- data.frame(zone=1:npoints,x=rndpoints$x,y=rndpoints$y)
			shppunti <- SpatialPoints(dfrandompoints[c("x","y")], proj4string=CRS(prjcrs), bbox = studyArea@bbox)
			buf2k = gBuffer(shppunti,byid=T,width=radius)
			rstbuf2k <- rasterize(buf2k, rstbuf)
			rstzonres <-zonal(rstclclvl3,rstbuf2k, sum)
			dftmp <- merge(dfrandompoints,rstzonres,by="zone")
			res.list[[i]] <- dftmp
			if(showtime==TRUE){
				print(paste("round=",i," - time(s)=", round((proc.time()-prt)[[3]],2)))
			}
			i <- i+1
		}
	}
	ris.fin=lapply(res.list,function(x) data.frame(PESO=sum(x$value),mindist=min(dist(x[c("x","y")]))))
	#ris.df=do.call("rbind",ris.fin)
	#(proc.time()-prt)[[3]]
	output <- do.call("rbind", ris.fin)					#output pesi
	#output.ord <- output[order(-output$PESO),] 		#pesi ordinati desc
	output.val = ris.fin[[which.max(output$PESO)]]		#lancio con MAX valore
	output.elem = res.list[[which.max(output$PESO)]]	#lista punti con MAX valore
	df.points <- data.frame(res.list[[which.max(output$PESO)]])
	
	print("END OF THE PROCESS!")
	print("")
	print(paste("BEST SIMULATION OF POINTS = ", which.max(output$PESO), ", WITH TOTAL VALUE = ", round(output$PESO[which.max(output$PESO)],2), sep= ""))
	print("")
	
	#exports
	print("STARTING EXPORT...")
	
	outSites <- SpatialPoints(df.points[c("x","y")], proj4string=CRS(prjcrs), bbox=studyArea@bbox)
	shpobj = SpatialPointsDataFrame(outSites@coords, data.frame(df.points))
	writeOGR(obj=shpobj, dsn=paste(outpath,"sites.shp",sep=""), layer="sites.shp", driver="ESRI Shapefile")
	print("sites.shp exported")
	
	outBuffers <- gBuffer(outSites,byid=T,width=radius)
	shpobj = SpatialPolygonsDataFrame(outBuffers, data.frame(dummy = rep(1,npoints)))
	writeOGR(obj=shpobj, dsn=paste(outpath,"buffers.shp",sep=""), layer="buffers.shp", driver="ESRI Shapefile")
	print("buffers.shp exported")
	print("EXPORT COMPLETE!")
	
	if(view == TRUE){
		windows(height=10,width=20)
		plot(rstclclvl3)
		plot(studyArea,add=TRUE)
		plot(outBuffers,add=TRUE,col="#6bc6ff",pch=1 )
		plot(outSites,add=TRUE,col="black",cex=0.5, pch=21)
		
	}
	outputs = list(rstclc = rstclclvl3, outSites = outSites, res.list = res.list, bestSimulation = which.max(output$PESO), bestValue = output$PESO[which.max(output$PESO)])
	writeRaster(outputs[[1]], paste(outpath, "clc_class.tif"), "GTiff", overwrite=TRUE)
	return(outputs)
}


#Check results
ahp.geo.check  = function(rstModelTest, list.data, studyArea, prjcrs="+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +units=m +no_defs"){
	#change NA to 0 in the check raster layer
	rstModelTest <- setValues(rstModelTest, ifelse(is.na(values(rstModelTest)),0,values(rstModelTest)))
	df.check <- data.frame(id = numeric(length=0), val = numeric(length=0), checkval = numeric(length=0), checkperc = numeric(length=0))

	for(i in 1:length(list.data)){
		spdata = SpatialPoints(list.data[[i]][2:3], proj4string=CRS(prjcrs), bbox = studyArea@bbox)
		spdatatmp <- data.frame(id = i, val = sum(list.data[[i]]$value), checkval = sum(extract(rstModelTest, spdata)), checkperc = sum(extract(rstModelTest, spdata))/nrow(list.data[[i]]))
		df.check <- rbind(df.check, spdatatmp)
	}
	return(df.check)
}
