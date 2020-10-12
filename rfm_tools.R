# Post-processing tools for RFM

#================================#
# Function rfm_i2s               #
#                                #
# Produces scalar values from    #
# from interface values,         #
# assuming interface z           # 
#================================#

rfm_i2s = function(f){   # assume native grid is interfacial,
	  		 # output has vertical length nz-1
      if (length(dim(f))==1){
         nz = dim(f)[1]
         fs = (f[-1] + f[-nz])/2
      } else if (length(dim(f))==2){
      	 nz = dim(f)[2]
         fs = (f[ ,-1] + f[ ,-nz])/2
      }
      return(fs)
      }

#================================#
# Function rfm_lni2s             #
#                                #
# Produces scalar values from    #
# interface values using         # 
# logarithmic interpolation      #
#================================#

rfm_lni2s = function(f){   # assume native grid is interfacial,
	   		   # output has vertical length nz-1
	  f = as.array(f)
      if (length(dim(f))==1){
         nz   = dim(f)[1]
	 lnf  = log(f)
         lnfs = (lnf[-1] + lnf[-nz])/2
      } else if (length(dim(f))==2){
      	 nz = dim(f)[2]
    	 lnfs = array(dim=c(dim(f)[1],nz-1))	 
         if (sum(f[ ,nz]==0) >0){  # in case zeros at top
	    lnf = log(f[ ,1:(nz-1)]) 
	    lnfs[ ,1:(nz-2)] = (lnf[ ,-1]+lnf[ ,-(nz-1)])/2
 	    lnfs[ ,nz-1]     = 2*lnfs[ ,nz-2] - lnfs[ ,nz-3]
	  } else {   # no zeros
	    lnf  = log(f) 
	    lnfs = (lnf[ ,-1]+lnf[ ,-nz])/2
	  }  # if zeros
      }	 # if length(dim) 
      fs = exp(lnfs)
      return(fs)
}


#========================#
# Function running_mean  #
#
# Computes running mean  #
# over interval nrun     #
#========================#

running_mean = function(x,nrun){
	if (length(dim(x))==1){
		nx    = length(x)
		xmean = 0
		for (i in 1:(nrun)){
		    xmean = xmean + x[i:(i+(nx-nrun))]
		    }
		xmean = xmean/nrun
	} else if (length(dim(x))==2){
		nx    = dim(x)[1]
		xmean = 0
		for (i in 1:(nrun)){
		    xmean = xmean + x[i:(i+(nx-nrun)),]
		}
		xmean = xmean/nrun
	}
	return(as.array(xmean))   # length nx-nrun+1
}

#========================#
# Function coarse_grain  #
#                        #
# Coarse-grains input by #
# a factor n             #
#========================#

coarse_grain = function(x,n){
				   x       = as.array(x)
				   nx      = dim(x)[1]
				   xmean   = running_mean(x,n)
				   ivec    = n*(0:(nx/n -1)) + 1	
		
				   if (length(dim(xmean))==1){
					   xcoarse = xmean[ivec]
					} else if (length(dim(xmean))==2){
					   xcoarse = xmean[ivec,]
					}
					return(xcoarse)
}


#========================#
# Function calc_pem      #
#                        #
# Calculates emission    #
# pressures for given    #
# tau_em                 #
#========================#

calc_pem = function(tau,p,tau_em,na.rm=F){
				   nk    = dim(tau)[1]
				   pem   = numeric(nk)	
				   for (m in 1:nk){
				   	  	pem[m] = p[which.min(abs(tau[m, ]-tau_em))]
				   }
				   if (na.rm){
				   	pem[pem==max(p)] <- NA
				   }
				   return(pem)
}

#========================#
# Function calc_Tem      #
#                        #
# Calculates emission    #
# temps for given tau_em #
#========================#

calc_Tem = function(tau,tabs,tau_em,na.rm=F){
				   nk    = dim(tau)[1]
				   Tem   = numeric(nk)	
				   for (m in 1:nk){
				   	  	Tem[m] = tabs[which.min(abs(tau[m, ]-tau_em))]
				   }
				   if (na.rm){
				   	Tem[Tem==max(p)] <- NA
				   }
				   return(Tem)
}

#==========================#
# Function get_olr         #
#                          #
# Calculates OLR from      # 
# RFM output               #
#==========================#

get_olr = function(ncpath,klim="NULL"){
		nc	   = nc_open(ncpath)
		k	   = ncvar_get(nc,"k")
		dk	   = diff(k)[1]
		nk	   = length(k)
		np	   = nc$dim$p$len
		start  = c(1,np)
		count  = c(nk,1)
		OLRk   = ncvar_get(nc,"flx",start=start,count=count)/1e2  # W/m^2/m^-1
		if (klim=="NULL"){
			OLR    = sum(dk*OLRk)
		} else {
			ivec = (which((k>=klim[1])&(k<=klim[2])))
			OLR    = sum(dk*OLRk[ivec])
		}
		nc_close(nc)
		return(OLR)
}

#===========================#
# Function get_forcing      #
#                           #
# Calculates OLR difference #
# from RFM output file vec  #
#===========================#

get_forcing = function(ncpaths,klim="NULL"){
		# ncpaths should be c(ncpath1,ncpath2)
		OLR1    = get_olr(ncpaths[1],klim)
		OLR2    = get_olr(ncpaths[2],klim)
		return(OLR1-OLR2)
}


#==========================#
# Functions calc_cts       #
#                          #
# Calculates 1D and 2D CTS # 
# fields from RFM output   #
#==========================#

calc_cts2d = function(k,p_i,tabs_s,tau_i){
		# NOTE: tau_i should already be scaled by D!
        nk        = dim(tau_i)[1]
        np        = dim(tau_i)[2]
        trans     = exp(-tau_i)
        dtrans_dp = (trans[ ,2:np]-trans[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
        cts2d     = -pi*outer(k,tabs_s,planck_k)*dtrans_dp  #s, W/m^2/Pa/m^-1, pppf units
        return(cts2d)
}

calc_cts1d = function(k,p_i,tabs_s,tau_i){
		# NOTE: tau_i should already be scaled by D!
        cts2d     = calc_cts2d(k,p_i,tabs_s,tau_i)
        dk        = diff(k)[1]
        cts1d     = apply(cts2d,2,sum)*dk
        return(cts1d)
}
