#===============#
# Tools for GCM #
# analysis      #
#===============#

source("~/Dropbox/Rtools/my_image_plot.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
latlab = "Latitude (deg)"

calc_col_pars = function(p,dp,q,tabs){
	       # assume all vars on s levs
	       col_pars = list() 	 

	       # ensure increasing p, GCM style
	       if (diff(p)[1] < 0){
	          for (var in c("p","q","tabs","dp")){
		      assign(var,rev(eval(as.name(var))))
		  }
	       }
	       qstar   = qsat(tabs,p)

	       # set mid-trop 
	       np         = length(p)  
	       k_mid_trop = floor(np*3/4)

 	       # Determine tp index k_tp as lowest inversion above mid-trop
               if (any(diff(tabs[1:k_mid_trop])<=0)) { 
               	  k_tp <- max(which(diff(tabs[1:k_mid_trop])<=0)) + 1
               } else {
               	  k_tp <- 1  # TOA
               } 

	       # Determine Tbot index k_bot as highest inversion below k_tp
	       if (all(diff(tabs[k_tp:np])>0)) {
               	  k_bot <- np
               } else {
               	  k_bot <- min(which(diff(tabs[k_tp:np])<=0)) - 1 + k_tp
               } 

               ptp    = p[k_tp]
	       	   pbot   = p[k_bot]
	           Ttp    = tabs[k_tp]
               Tbot   = tabs[k_bot]
               if (k_tp >= k_bot ){	
		  return(NA)
	       }  else {
	          kvec   = k_tp:np  # integrate through inversion at k_bot for CRH!
      	      col_pars$CRH  = (dp[kvec]%*%q[kvec])/(dp[kvec]%*%qstar[kvec])
	       	  #col_pars$Gamma = g/Rd*log(Ttp/tabs[np])/log(ptp/p[np]) # down to np, not kbot 
	       	  col_pars$Gamma = g/Rd*log(Ttp/Tbot)/log(ptp/pbot) # down to kbot 
	       	  col_pars$ptp   = ptp
		  	  col_pars$Ttp   = Ttp
	       	  return(col_pars)
	       }
}	       


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
		return(mean)
		}

plot_map <- function(lon,lat,field,main,zlim=range(field,na.rm=TRUE),
					do_mean = TRUE, units = NULL, horizontal = FALSE,
					cex.lab = 1.15, cex.axis=1.15,cex.main=1.15,
					cex_mean=1.15, cex_leg = 1.15,add.legend=TRUE){
           	my.image.plot(lon,lat,field,zlim=zlim,
                        xlab = "Longitude (deg)",
						ylab = "Latitude (deg)", 
                        main = main,
						cex.lab = cex.lab,
						cex.main = cex.main,
						cex.axis = cex.axis,
						cex.legend  = cex_leg,
						horizontal  = horizontal,
						add.legend  = add.legend)
            map("world2",add=TRUE,interior=FALSE)
            title("")
	    if (do_mean){
	        mean = global_mean(lon,lat,field)
            	mtext(bquote("Global mean = "~.(round(mean,1))),
            	     at=265,cex=cex_mean)
	    }	     
            mtext(units,3,at=350,cex=1.15)
            }


