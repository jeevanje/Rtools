#================#
# Functions for  #
# computing      #
# gray radiative #
# fluxes         #
#================#

Rtoolsdir = "~/Dropbox/Rtools"
source(paste(Rtoolsdir,"/calculus_tools.R",sep=""))
source(paste(Rtoolsdir,"/thermo_tools.R",sep=""))

#===============#
# Optical depth #
# routines      #
#===============#

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

compute_trans = function(tau){
	        # All quantities ssi!
	        nz = length(tau)
	      	trans = array(dim=c(nz,nz))
		trans = outer(exp(-tau),exp(tau),"*")
		}

#====================#
# Radiative transfer #
#====================#

compute_Ugray = function(tau,tabs,Ts=tabs[1]){
	    # tau ssi, tabs sss
	    ntau  = length(tau)
	    dtau  = diff(c(tau,0))
	    Ugray = numeric(ntau)
		Ugray[1] = B(Ts)
		for (i in 1:(ntau-1)){
		    # see rad_cooling notes 10/17
		    Ugray[i+1] = ( (1+dtau[i]/2)*Ugray[i] -dtau[i]*B(tabs[i]) ) /(1-dtau[i]/2)
		    }
		 return(Ugray)
	}

compute_Dgray = function(tau,tabs){
	    # tau ssi, tabs sss
	    ntau  = length(tau)
	    dtau  = diff(c(tau,0))
	    Dgray = numeric(ntau)
		Dgray[ntau] = 0  # neglect emission from topmost grid cell
		for (i in (ntau-1):1){
		    # see rad_cooling notes 10/17
		    Dgray[i] = ( (1+dtau[i]/2)*Dgray[i+1] -dtau[i]*B(tabs[i]) ) /(1 - dtau[i]/2)
		    }
		 return(Dgray)
	}

#===============#
# Heating rates #
#===============#

compute_ex   = function(Bvals,tau){
	       # tau ssi, Bvals sss	       
		   epsilon=1e-10  # fudge for weird R bug in k_2tau
	       N     = length(tau)
		   ex    = list()
		   ex[['ax']] <- numeric(N)
		   ex[['sx']] <- numeric(N)		   
	       taus  = tau[1]
	       tau_s = c(tau[-N] + diff(tau)/2,tau[N]+diff(tau)[N-1]/2)
	       dtau  = -c(diff(tau),0) # slight fudge here
	       for (k in 1:N){   # k is s levels
 	       	   if (tau_s[k] < taus/2){
				  # ax		
	       	      k_2tau = max(which(tau>2*tau_s[k]))  # i level
				  kvec   = 1:(k_2tau-1)  #s levs
	       	      ax =  -sum(dtau[kvec]*(Bvals[kvec]-Bvals[k])*exp(-(tau_s[kvec]-tau_s[k])))
				  ex[['ax']][k] <- ax
				   
				  # sx	
				  kvec1  = (k+1):N
	       	      sx1    = -sum(dtau[kvec1]*(Bvals[kvec1]-Bvals[k])*exp(-(tau_s[k]-tau_s[kvec1])))
				  kvec2  = k_2tau:(k-1)	
				  sx2	 = -sum(dtau[kvec2]*(Bvals[kvec2]-Bvals[k])*exp(-(tau_s[kvec2]-tau_s[k])))		       	      
				  ex[['sx']][k]  <- sx1 + sx2
				  	
				} else if (tau_s[k]>=taus/2){
				  # ax	
	       	      k_2tau = min(which(tau <= (2*tau_s[k]-taus + epsilon)))  # i level
				  kvec   = k_2tau:N  #s levs
	       	      ex[['ax']][k]  <- -sum(dtau[kvec]*(Bvals[kvec]-Bvals[k])*exp(-(tau_s[k]-tau_s[kvec])))
	       	      
	       	      # sx
				  kvec1  = 1:(k-1)
	       	      sx1    = -sum(dtau[kvec1]*(Bvals[kvec1]-Bvals[k])*exp(-(tau_s[kvec1]-tau_s[k])))
				  kvec2  = (k+1):(k_2tau-1)	
				  sx2	 = -sum(dtau[kvec2]*(Bvals[kvec2]-Bvals[k])*exp(-(tau_s[k]-tau_s[kvec2])))		       	      
				  ex[['sx']][k]  <- sx1 + sx2
				}		      
			}   # k
			return(ex)
}

 
compute_pptf = function(kappa,z,U,D,tabs,rhov,lapse) {
	     U_sss = i2s(3,z,U)
	     D_sss = i2s(3,z,D)
	     pptf = kappa/lapse*rhov*(U_sss + D_sss - 2*B(tabs))
	     return(pptf)
	     }

#=======#
# other #
#=======# 

tune_kappa = function(z,rhov,tabs,sst,olr){
	     nz = length(z)
	     cost = function(kappa){
	     	     U = compute_Ugray(kappa,z,rhov,tabs,sst)
		     return(U[nz]-olr)
		     }
	     kappa = uniroot(cost,interval=c(1e-2,1e1) )[[1]]
	     return(kappa)
	     }
