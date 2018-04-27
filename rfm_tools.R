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
      if (length(dim(f))==1){
         nz   = dim(f)[1]
	 lnf  = log(f)
         lnfs = (lnf[-1] + lnf[-nz])/2
	 fs   = exp(lnfs)
      } else if (length(dim(f))==2){
      	 nz = dim(f)[2]
	 if (sum(f[ ,nz]==0) >0){
	    lnfs = array(dim=c(dim(f)[1],dim(f)[2]-1))
	    lnf = log(f[ ,1:(nz-1)]) 
	    lnfs[ ,1:(nz-2)] = (lnf[ ,-1]+lnf[ ,-(nz-1)])/2
 	    lnfs[ ,nz-1]     = 2*lnfs[ ,nz-2] - lnfs[ ,nz-3]
	    fs	= exp(lnfs)
	  }
      }	  
      return(fs)
}


#========================#
# Function running_mean  #
#
# Computes running mean  #
# over interval nrun     #
#========================#

running_mean = function(x,nrun){
	nx    = length(x)
	xmean = 0
	for (i in 1:(nrun)){
	    xmean = xmean + x[i:(i+(nx-nrun))]
	    }
	xmean = xmean/nrun
	return(xmean)   # length nx-nrun+1
}

