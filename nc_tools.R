#============#
# Tools for  #
# processing #
# netcdf     #
# files      #
#============#

source("~/Dropbox/Rtools/calculus_tools.R")

get_profile = function(nc,period,varname){
	      		# period assumed in days
		        time = get.var.ncdf(nc,"time")
				time_start = max(time) - period
				if (time_start < 0 ){
				   print("error: period longer than record")
				   return()
				}
				nt_start = which.min(abs(time-time_start))
		        var2D = get.var.ncdf(nc,start=c(1,nt_start),varname)
				var = apply(var2D,1,mean)
				return(var)      
				}

get_slice = function(nc,nt_avg,varname){
	        time = get.var.ncdf(nc,"time")
		nt = length(time)
	        var3D = get.var.ncdf(nc,start=c(1,1,nt-nt_avg),varname)
		var = apply(var3D,1:2,mean)
		return(var)      
		}

get_streamfunction = function(nc,nt_avg){
	        z    = get.var.ncdf(nc,"z")		   
		nz   = length(z)
		zint = zinterp(z)
		dz   = c(diff(zint),z[nz]-z[nz-1])
		rho    = get_slice(nc,nt_avg,"rho") 
		v    = get_slice(nc,nt_avg,"v") 
		w    = get_slice(nc,nt_avg,"w") 
		ny   = dim(v)[1]
		psi  = array(dim=c(ny,nz))
		psi[ ,1]<-(rho*v)[ ,1]*dz[1]
		for (j in 2:nz){
			psi[,j]<- ((rho*v)[ ,1:j])%*%dz[1:j]
			}
		return(psi)
		}


get_avg     = function(nc,period,varname){    
	      		# period assumed in days
		        time = get.var.ncdf(nc,"time")
				time_start = max(time) - period
				if (time_start < 0 ){
				   print("error: period longer than record")
				   return()
				}
				nt_start = which.min(abs(time-time_start))
                var1D = get.var.ncdf(nc,start=c(nt_start),varname)
                var = mean(var1D)
                return(var)
                }                            

