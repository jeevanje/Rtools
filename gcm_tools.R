#===============#
# Tools for GCM #
# analysis      #
#===============#

global_mean = function(lon,lat,field){
		# assumes lat in degrees!
		R    = 6.36e6  # m
		nlon = length(lon)
		nlat = length(lat)
		dlon = lon[2]-lon[1]
		dlat = lat[2]-lat[1]
		dlonvec = rep.int(dlon,nlon)
		dlatvec= rep.int(dlat,nlat)
		area = (pi/180)^2*R^2*dlonvec%o%(cos(pi/180*lat)*dlatvec)
		mean = sum(field*area)/sum(area)
		mean
		}

