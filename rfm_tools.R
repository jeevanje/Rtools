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

#=======================#
# Function read_atm     #
#                       #
# Reads .atm RFM input  #
# and outputs single    #
# atmospheric field     #
#=======================#

read_atm = function(atmpath,skip,nlev){
	var = scan(atmpath,skip=skip,sep=",",nlines=ceiling(nlev/5),
			strip.white=TRUE,skipNul=TRUE) #km
        var = var[!is.na(var)]
        dim(var) <- nlev						       
        return(var)
}

#=======================#
# Function read_asc     #
#                       #
# Reads .asc RFM output #
# and outputs list with #
# thermodynamic fields  #
# as well as RFM output #
#=======================#

read_asc = function(case,type){
    data       = list()
    varlist    = c("k","z","p","tabs","field")

    rfmdir     = "~/rad_cooling2/rfm"
    drvfile    = paste(rfmdir,"/",case,"/rfm_",type,".drv",sep="")
    datadir    = paste(rfmdir,case,type,sep="/")
    atmfile    = scan(drvfile,skip=9,sep="/",nmax=3,what="raw")[3]
    atmpath    = paste(rfmdir,"/atm/",atmfile,sep="")

    nlev       = as.numeric(scan(atmpath,skip=2,nmax=1)[1])
    nlines     = ceiling(nlev/5)
    z  	       = 1e3*read_atm(atmpath,skip=4,nlev=nlev)   #m
    p	       = 1e2*read_atm(atmpath,skip=4+(1+nlines),nlev=nlev)  #Pa
    tabs       = read_atm(atmpath,skip=4+2*(1+nlines),nlev=nlev)    

    files      = list.files(datadir,pattern=type)
    file_ref   = paste(datadir,files[1],sep="/")
    metadata   = scan(file_ref,skip=3,nmax=5,what="raw")
    nk	       = as.numeric(metadata[1])
    k1	       = 1e2*as.numeric(metadata[2])  # m^-1
    dk	       = 1e2*as.numeric(metadata[3])  # m^-1
    k_nk       = 1e2*as.numeric(metadata[4])  # m^-1
    k	       = seq(from = k1, to = k_nk,by=dk)
    
    field      = array(dim=c(nk,nlev))

    for (m in 1:nlev){
       zval      = z[m]    # km
       if (zval < 1e5){
	  zstring = formatC(zval,format="d",width=5,flag="0")
       } else {
	  zstring = formatC(zval/1000,format="d",width=5,flag="0")
       }
       file   = paste(datadir,"/",type,"_",zstring,".asc",sep="")
       field[ ,m] = scan(file,skip=4)
    }
    for (var in varlist){
       
       data[[var]] = eval(as.name(var))
    }
    return(data)
}

#=====================#
# Function make_atm   #
#                     #
# Make atm file from  #
# given profiles      #
#=====================#

make_atm = function(file,description,z,p,tabs,h2o,co2){
	 # Preliminaries
	 labs  = c("*HGT [km]","*PRE [mb]","*TEM [K]",
	 	    "*H2O [ppmv]","*CO2 [ppmv]")
	 vars  = c("z","p","tabs","h2o","co2")  # SI units
	 facs  = c(1e-3,1e-2,1,1e6,1e6)		# For conversion
	 nvars = length(vars)
	 write_file = function(record){
	 	   write(record,file=file,append=TRUE,sep=", ")
         }
	 nlev = length(z)

	 # Write
	 write_file(description)
	 write_file("! Produced by make_atm in rfm_tools.R")
	 write_file(nlev)
	 for (i in 1:nvars){ 
	    lab  = labs[i]
	    var  = vars[i]
	    fac  = facs[i]
	    write_file(lab)
	    write_file(fac*eval(as.name(var)))
	 }
	 write_file("*END")
}