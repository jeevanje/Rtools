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
theta_calc <-function(tabs,p){
      tabs*(ps/p)^(Rd/Cp)
      }


#======================#
# Theta_e computation  #
#======================#
mix_ratio<-function(q){
		q/(1-q)
		}
rho<-function(q,p,tabs){
	p/( (1-q)*(Rd*tabs+q/(1-q)*Rv*tabs) )
	}
e<-function(q,p){
	q*p/epsilon
	}
RH<-function(q,p,tabs){
	e(q,p)/esat(tabs)
	}
thetae<-function(q,p,tabs,p0=1e5){
	pd = p-e(q,p)
	tabs*(p0/pd)^(Rd/Cp)*RH(q,p,tabs)^(-Rv*q/(1-q)/Cp)*
		exp(L0*q/(1-q)/(Cp*tabs)) 
	}


#===================#
# Buoyancy function #
#===================#
# Construct buoyancy, refer to far-field values (not domain mean)
# Only ingest 3D fields for now

buoyancy <- function(rho){
	 nx     = (dim(rho))[1]
	 ny     = (dim(rho))[2]
	 nz     = (dim(rho))[3]
	 rhobar = rho[nx,ny,]   # far-field values
	 buoy   = array(dim=c(nx,ny,nz))
	 for (k in 1:nz){
	     buoy[ , ,k]<- g*(rhobar[k]-rho[ , ,k])/rhobar[k]
	     }
	  buoy
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
		     
