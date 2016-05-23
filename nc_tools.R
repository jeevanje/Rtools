#============#
# Tools for  #
# processing #
# netcdf     #
# files      #
#============#

get_profile = function(nc,nt_avg,varname){
	        time = get.var.ncdf(nc,"time")
		nt = length(time)
	        var2D = get.var.ncdf(nc,start=c(1,nt-nt_avg),varname)
		var = apply(var2D,1,mean)
		return(var)      
		}