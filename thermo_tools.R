####################
# Thermo_tools     #
#                  #
# R tools for      #
# thermodynamic    #
# computations     # 
####################

Rtools_dir = "~/Dropbox/Rtools/"
source(paste(Rtools_dir,"calculus_tools.R",sep=""))

L     = 2.5e6      # J/kg
L0    = 2.555e6    # J/kg, see Bryan (2008)
Rv    = 461.5      # J/kg/K
Rd    = 287        # J/kg/K
epsilon = Rd/Rv
Cp    = 1004       # J/kg/K
einf  = 2.53e11 # Pa
g     = 9.81       # m/s^2, consistent with DAM
ps    = 1e5        # Pa
sigmaSB=5.67e-8
h     = 6.626e-34  # Js
c     = 3e8        # m/s
k_b   = 1.38e-23   # J/K
N_avo = 6.022e23 
m_air = 28.971e-3  # molar mass, kg/mol
m_co2 = 44e-3
m_h2o = 18e-3

#==================#
# Basic quantities #
#==================#
esat<-function(tabs){
	einf*exp(-L/(Rv*tabs)) # in Pa
	}

qsat<-function(tabs,p){
 epsilon*esat(tabs)/p #p in Pa
	} # in kg/kg
	
gamma_m<-function(tabs,p){
		g*(1+qsat(tabs,p)*L/(Rd*tabs))/( Cp+qsat(tabs,p)*L^2/(Rv*tabs^2) )
		}
gamma_qv<-function(tabs,p){
		L*gamma_m(tabs,p)/(Rv*tabs^2)-g/(Rd*tabs) 		
		}
theta_func <-function(tabs,p){
      tabs*(ps/p)^(Rd/Cp)
      }

k2f <- function(tabs){
		return((tabs - 273.15)*9/5 +32)
}

k2c <- function(tabs){
		return(tabs - 273.15)
}

#================#
# moist adiabat  #
#================#
make_adiabat = function(z,Ts,gamma,Ttp=200){  
	     # Note: i,s levs not consistent
	     nz = length(z)
	     if (gamma != "m"){       # gamma = "m" or value in K/m
		tabs  = Ts - gamma*z
		p     = ps*(tabs/Ts)^(g/gamma/Rd)
	     } else if (gamma == "m"){
	        tabs    = numeric(nz)
		p       = numeric(nz)
		tabs[1] = Ts
		p[1]    = ps
		for (k in 2:nz){
		    tabs[k] = tabs[k-1] - gamma_m(tabs[k-1],p[k-1])*(diff(z)[k-1])
		    if (tabs[k] < Ttp){ tabs[k]=Ttp }
		    p[k]    = p[k-1]    - g*p[k-1]/Rd/tabs[k-1]*(diff(z)[k-1])
   	      	}
	     }   	
	     return(list(tabs,p))
}

#======================#
# Theta_e computation  #
#======================#
mix_ratio<-function(q){
		q/(1-q)
		}
rho_func<-function(q,p,tabs){
	p/( (1-q)*(Rd*tabs+q/(1-q)*Rv*tabs) )
	}
e<-function(q,p){
	q*p/epsilon
	}
RH_func<-function(q,p,tabs){
	e(q,p)/esat(tabs)
	}
thetae_func<-function(q,p,tabs,p0=1e5){
	pd = p-e(q,p)
	tabs*(p0/pd)^(Rd/Cp)*RH_func(q,p,tabs)^(-Rv*q/(1-q)/Cp)*
		exp(L0*q/(1-q)/(Cp*tabs)) 
	}

#=====#
# WVP #
#=====#

compute_wvp   = function(z,rhov){  #z, rhov on s levs
            nz = length(z)
            zint = zinterp(z)
            dz   = diff(zint)  # length = nz-1
			if (length(dim(rhov))==3){
	            wvp  = array(dim=dim(rhov)[1:2])
    	        wvp   = rhov[ , ,nz]*dz[nz-1]  # extrapolate dzvec
        	    for (k in (nz-1):1){
            	     wvp = wvp + dz[k]*rhov[ , ,k]
            	}
            } else if (length(dim(rhov))==1){
            	wvp  = rhov[nz]*dz[nz-1]
            	wvp  = wvp +rhov[-nz]%*%dz
            }
            return(wvp)
}



#====================#
# Buoyancy functions #
#====================#

buoyancy3d <- function(rho3d,far_field=TRUE){
	 nx     = (dim(rho3d))[1]
	 ny     = (dim(rho3d))[2]
	 nz     = (dim(rho3d))[3]
	 B3d   = array(dim=c(nx,ny,nz))
	 # compute reference density rhobar 
	 if (far_field){
	    rhobar = rho3d[nx,ny,]       # far-field values
	    } else {
	    rhobar = apply(rho3d,3,mean) # domain-mean
	    }
	 # compute B
	 for (k in 1:nz){
	     B3d[ , ,k]<- g*(rhobar[k]-rho3d[ , ,k])/rhobar[k]
	     }
	  return(B3d)
	  }

buoyancy4d <- function(rho4d,far_field=TRUE){
	   nx = (dim(rho4d))[1]
	   ny = (dim(rho4d))[2]
	   nz = (dim(rho4d))[3]
	   nt = (dim(rho4d))[4]
     	   B4d   = array(dim=c(nx,ny,nz,nt))	   
	   for (l in 1:nt){
	       B4d[ , , ,l] = buoyancy3d(rho4d[ , , ,l],far_field)
	       }
	    return(B4d)
	    }
	       
buoyancy_func <- function(rho,far_field=TRUE){
	      	 if (length(dim(rho)) == 4){
		    return(buoyancy4d(rho,far_field))
		    } else if (length(dim(rho)) == 3) {
		    return(buoyancy3d(rho,far_field))
		    } else {
		    print("dim(rho) not 3 or 4. Exit.")
		    stop
		    }
		 }

#==================#
# Planck functions #
#==================#

planck_lambda = function(lambda,T){
		return(2*h*c^2/lambda^5/(exp(h*c/(lambda*k_b*T))-1) )
	  	}   # W/m^2/sr/m

planck_nu     = function(nu,T){
	      	return(2*h*nu^3/c^2/(exp(h*nu/(k_b*T))-1) )
	      	}   # W/m^2/sr/Hz 

# wavenumber k = 1/lambda
planck_k      = function(k,T){
	      						return(2*h*c^2*k^3/(exp(h*c*k/(k_b*T))-1) )
								# W/m^2/sr/m^-1 
				}   

B	      	  = function(T){
	      	return(sigmaSB*T^4)
		}	 
		
planck_k_linear = function(T,k){
	      			return(2*k_b*c*k^2*T)
								# W/m^2/sr/m^-1 
				}   

planck_k_interval = function(T,k1,k2){
					kintvals = 1e2*((1e-2*k1):(1e-2*k2))  # dk = 1 cm^-1, int vals
#					kintvals = seq(k1,k2,length.out=1e5)
					n_kint 	 = length(kintvals)
		     		kvals	 = (kintvals[1:(n_kint-1)] + kintvals[2:n_kint])/2
					dkvec	 = diff(kintvals)
					return(pi*t(planck_k(kvals,T))%*%dkvec)  # W/m^2
}		     
