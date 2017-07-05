####################
# Thermo_tools     #
#                  #
# R tools for      #
# thermodynamic    #
# computations     # 
####################

L = 2.5e6      # J/kg
L0= 2.555e6    # J/kg, see Bryan (2008)
Rv = 461.5     # J/kg/K
Rd = 287       # J/kg/K
epsilon = Rd/Rv
Cp = 1004      # J/kg/K
einf = 2.53e11 # Pa
g = 9.81       # m/s^2, consistent with DAM
ps= 1e5        # Pa
sigmaSB=5.67e-8
h = 6.626e-34    # Js
c = 3e8        # m/s
k_b = 1.38e-23 # J/K


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

planck_lambda = function(T,lambda){
		return(2*h*c^2/lambda^5/(exp(h*c/(lambda*k_b*T))-1) )
	  	}   # W/m^2/sr/m

planck_nu     = function(T,nu){
	      	return(2*h*nu^3/c^2/(exp(h*nu/(k_b*T))-1) )
	      	}   # W/m^2/sr/Hz 

# wavenumber k = 1/lambda
planck_k      = function(T,k){
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
					return(pi*planck_k(T,kvals)%*%dkvec)  # W/m^2
}		     
