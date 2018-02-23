#===============#
# Tools for GCM #
# analysis      #
#===============#

source("~/Dropbox/Rtools/my_image_plot.R")

global_mean = function(lon,lat,field){
		# assumes lat/lon in degrees!
		lat_rad  = pi/180*lat
		lon_rad  = pi/180*lon
		nlon     = length(lon)	
		nlat     = length(lat)
		dlon_rad = lon_rad[2]-lon_rad[1]
		dlat_rad = lat_rad[2]-lat_rad[1]
		dlonvec  = rep.int(dlon_rad,nlon)
		dlatvec  = rep.int(dlat_rad,nlat)
		if (length(dim(field)) ==3 ){
		   n3d = dim(field)[3]
		   d3dvec = rep.int(1,n3d)
		   domega = dlonvec%o%(cos(lat_rad)*dlatvec)%o%d3dvec
		   } else if ( length(dim(field)) ==2 ){	
		   domega = dlonvec%o%(cos(lat_rad)*dlatvec) 
		   } else {
		   print("dim(field) neq 2 or 3. Stopping")
		   return()
		   }	
		domega[is.na(field)] = NA
		mean = sum(field*domega,na.rm=TRUE)/sum(domega,na.rm=TRUE)
		mean
		}

plot_map <- function(lon,lat,field,main,zlim=range(field),
			do_mean = TRUE, units = NULL,cex = 1.15,cex_leg = 1.15){
            my.image.plot(lon,lat,field,zlim=zlim,
                        xlab = "Longitude (deg)",
			ylab = "Latitude (deg)", 
                        main = main,
			cex.lab = cex,
			cex.main = cex,
			cex.axis = cex,
			cex.legend  = cex_leg )
            map("world2",add=TRUE,interior=FALSE)
            title("")
            mean = global_mean(lon,lat,field)
            mtext(bquote("Global mean = "~.(round(mean,1))~.(units)),at=360)
            }


