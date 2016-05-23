#================#
# Functions for  #
# computing      #
# gray radiative #
# fluxes         #
#================#

Rtoolsdir = "~/Dropbox/Rtools"
source(paste(Rtoolsdir,"/calculus_tools.R",sep=""))
source(paste(Rtoolsdir,"/thermo_tools.R",sep=""))

compute_wvp   = function(z,rhov){
	      	#ssi !
	      	nz = length(z)
		zint = zinterp(z)
		dz   = diff(zint)  # length = nz-1
		wvp = numeric(nz)
		wvp[nz] = rhov[nz]*dz[nz-1]  # extrapolate dzvec
		for (k in (nz-1):1){
                    wvp[k] = wvp[k+1] + dz[k]*rhov[k]
		    } 
		return(wvp)
		}

compute_Ugray = function(kappa,z,rhov,tabs,sst){
	        nz = length(z)
		zint = zinterp(z)
		dz   = diff(zint)  # length = nz-1
	      	Ugray = numeric(nz)
		Ugray[1] = B(sst)
		for (i in 1:(nz-1)){
		    # see rad_cooling notes 5/6/16
		    Ugray[i+1] = ( (1-dz[i]*kappa*rhov[i]/2)*Ugray[i] + 
		    	          kappa*rhov[i]*B(tabs[i])*dz[i] ) /
				  (1 + kappa*rhov[i]*dz[i]/2)
		    }
		 return(Ugray)
	}

compute_Dgray = function(kappa,z,rhov,tabs){
	        nz = length(z)
		zint = zinterp(z)
		dz   = diff(zint)  # length = nz-1
	      	Dgray = numeric(nz)
		Dgray[nz] = 0  # neglect emission from topmost grid cell
		for (i in (nz-1):1){
		    # see rad_cooling notes 5/6/16
		    Dgray[i] = ( (1-dz[i]*kappa*rhov[i]/2)*Dgray[i+1] + 
		    	          kappa*rhov[i]*B(tabs[i])*dz[i] ) /
				  (1 + kappa*rhov[i]*dz[i]/2)
		    }
		 return(Dgray)
	}

compute_pptf = function(kappa,z,U,D,tabs,rhov,lapse) {
	     U_sss = i2s(3,z,U)
	     D_sss = i2s(3,z,D)
	     pptf = kappa/lapse*rhov*(U_sss + D_sss - 2*B(tabs))
	     return(pptf)
	     }

compute_trans = function(tau){
	        # All quantities ssi!
	        nz = length(tau)
	      	trans = array(dim=c(nz,nz))
		trans = outer(exp(-tau),exp(tau),"*")
		}

tune_kappa = function(z,rhov,tabs,sst,olr){
	     nz = length(z)
	     cost = function(kappa){
	     	     U = compute_Ugray(kappa,z,rhov,tabs,sst)
		     return(U[nz]-olr)
		     }
	     kappa = uniroot(cost,interval=c(1e-2,1e1) )[[1]]
	     return(kappa)
	     }
