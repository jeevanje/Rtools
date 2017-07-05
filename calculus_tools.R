  #================================#
  # Function zinterp               #
  #                                #
  # Produces interface values of   #
  # z from scalar values           # 
  #================================#

  zinterp<-function(z){
      nz=length(z)
      zint<-numeric()
      zint[1]=0
	  zint[2:nz] = 0.5*(z[-1]+z[-nz])
      return(zint)
      }

  #================================#
  # Function pinterp               #
  #                                #
  # Produces interface values of   #
  # p from scalar values           # 
  #================================#

  pinterp<-function(p){
	  # Assumes p increasing with k
      np=length(p)
      lnp = log(p)
      lnp_int<-numeric(np)
	  lnp_int[-np] = 0.5*(lnp[-1]+lnp[-np])
	  lnp_int[np] = 1.5*lnp[np] - 0.5*lnp[np-1]
      return(exp(lnp_int))
      }

  #================================#
  # Function i2s                   #
  #                                #
  # Moves field from               #
  # interface positions to scalar  #
  # with Dirichlet BCs             #
  #================================#

  i2s<-function(p,coord,f){

    # p=Coordinate rank, coord=coordinate vector, f=3d field

    # Quickly do 1D calculation
    if (length(f) == length(coord)) {
       nz = length(coord)
       i2s<-numeric(length(coord))
       # Initialize  zhalo
       zhalo = c( -coord[1],coord, 2*coord[nz]-coord[nz-1] )

       # Construct fhalo in z-direction with f = 0 on bottom boundary
       fhalo<-numeric(nz+2)
       fhalo[2:(nz+1)] = f[1:nz]
       fhalo[nz+2] = 0
    
       # Compute i2s
       for (k in 1:nz) {
 			i2s[k] = 1/(zhalo[k+2]-zhalo[k])*
			   ( fhalo[k+1]*(zhalo[k+2]-zhalo[k+1]) + 
			     fhalo[k+2]*(zhalo[k+1]-zhalo[k] ) )      
          }   # k  
       return(i2s)
       }  # if 1D


    # Get dimension sizes
    nx = dim(f)[1]
    ny = dim(f)[2]
    nz = dim(f)[3]
    i2s<-array(dim=c(nx,ny,nz))

    if (p == 1) {        # ppix

       # Construct f-halo on x-direction
       fhalo<-array(dim=c(nx+2,ny,nz))
       fhalo[2:(nx+1),1:ny,1:nz] = f
       fhalo[1, , ] = f[nx, , ]
       fhalo[nx+2, , ] = f[1, , ]

       # Compute i2s
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                i2s[i,j,k] = 1/2*( fhalo[i+2,j,k] + 
				      	fhalo[i+1,j,k] )
		# Klugey fudge in indexing between LHS and RHS
             }
           }
	}	  

    } else if (p == 2) {                  # ppiy

       # Construct f-halo on y-direction
       fhalo<-array(dim=c(nx,ny+2,nz))
       fhalo[1:nx,2:(ny+1),1:nz] = f
       fhalo[ ,1, ] = f[ ,ny, ]
       fhalo[ ,ny+2 , ] = f[ ,1 , ]

       # Compute i2s
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                i2s[i,j,k] = 1/2*( fhalo[i,j+2,k] + fhalo[i,j+1,k] )
                 		# Klugey fudge in indexing between LHS and RHS
             }
          }
       }


    } else if (p == 3) {        # ppiz

       # Initialize  zhalo
       zhalo<-numeric()
       zhalo[2:(nz+1)] = coord
       zhalo[1] = -coord[1]
       zhalo[nz+2] = 2.*coord[nz]-coord[nz-1]

       # Construct fhalo in z-direction with f = 0 on upper boundary
       fhalo<-array(dim=c(nx,ny,nz+2))
       fhalo[ , ,2:(nz+1)] = f
       fhalo[ , ,nz+2] = 0
    
       # Compute i2s
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                i2s[i,j,k] = 1/(zhalo[k+2]-zhalo[k])*
			   ( fhalo[i,j,k+1]*(zhalo[k+2]-zhalo[k+1]) + 
			     fhalo[i,j,k+2]*(zhalo[k+1]-zhalo[k] ) )
             }
          }     
       }
    }

    i2s
}

  #================================#
  # Function s2i                   #
  #                                #
  # Moves field from               #
  # scalar positions to interface  #
  # with Dirichlet BCs             #
  #================================#

  s2i<-function(p,coord,f){

    # p=Coordinate rank, coord=coordinate vector, f=3d field

 # Quickly do 1D calculation
    if (length(f) == length(coord)) {
       nz = length(coord)
       s2i<-numeric(length(coord))

       zhalo = c( -coord[1],coord, 2*coord[nz]-coord[nz-1] )

       # Construct fhalo in z-direction with f = 0 on bottom boundary
       fhalo<-numeric(nz+2)
       fhalo[2:(nz+1)] = f[1:nz]
       fhalo[1] = -f[1]
    
       # Compute s2i
       for (k in 1:nz) {
                s2i[k] = 1/2*(fhalo[k+1]+fhalo[k])
		}
       return(s2i)
       }  # if 1D

    # Get dimension sizes
    nx = dim(f)[1]
    ny = dim(f)[2]
    nz = dim(f)[3]
    s2i<-array(dim=c(nx,ny,nz))

    if (p == 1) {        # ppix

       # Construct f-halo on x-direction
       fhalo<-array(dim=c(nx+2,ny,nz))
       fhalo[2:(nx+1),1:ny,1:nz] = f
       fhalo[1, , ] = f[nx, , ]
       fhalo[nx+2, , ] = f[1, , ]

       # Compute s2i
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                s2i[i,j,k] = 1/2*( fhalo[i+1,j,k] + 
				      	fhalo[i,j,k] )
		# Klugey fudge in indexing between LHS and RHS
             }
           }
	}	  

    } else if (p == 2) {                  # ppiy

       # Construct f-halo on y-direction
       fhalo<-array(dim=c(nx,ny+2,nz))
       fhalo[1:nx,2:(ny+1),1:nz] = f
       fhalo[ ,1, ] = f[ ,ny, ]
       fhalo[ ,ny+2 , ] = f[ ,1 , ]

       # Compute s2i
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                s2i[i,j,k] = 1/2*( fhalo[i,j+1,k] + fhalo[i,j,k] )
                 		# Klugey fudge in indexing between LHS and RHS
             }
          }
       }


    } else if (p == 3) {        # ppiz

       # Construct fhalo in z-direction with f = 0 on bottom boundary
       fhalo<-array(dim=c(nx,ny,nz+2))
       fhalo[ , ,2:(nz+1)] = f
       fhalo[ , ,1] = -fhalo[ , ,2]
    
       # Compute s2i
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                s2i[i,j,k] = 1/2*(fhalo[i,j,k+1]+fhalo[i,j,k])
             }
          }     
       }
    }
    s2i
}

  #================================#
  # Function i2s4d                 #
  #                                #
  # Moves 4d field from            #
  # interface positions to scalar  #
  # with Dirichlet BCs             #
  #================================#

  i2s4d<-function(p,coord,f){

    # Get dimension sizes
    nx = dim(f)[1]
    ny = dim(f)[2]
    nz = dim(f)[3]
    nt = dim(f)[4]
    i2s4d<-array(dim=c(nx,ny,nz,nt))
    for (l in 1:nt){
    	i2s4d[ , , ,l]=i2s(p,coord,f[ , , ,l])
    }
    i2s4d
  }

  #================================#
  # Function s2i4d                 #
  #                                #
  # Moves 4d field from            #
  # scalar positions to interface  #
  # with Dirichlet BCs             #
  #================================#


  s2i4d<-function(p,coord,f){

    # Get dimension sizes
    nx = dim(f)[1]
    ny = dim(f)[2]
    nz = dim(f)[3]
    nt = dim(f)[4]
    s2i4d<-array(dim=c(nx,ny,nz,nt))
    for (l in 1:nt){
    	s2i4d[ , , ,l]=s2i(p,coord,f[ , , ,l])
    }
    s2i4d
  }

  #================================#
  # Function partialder_s2i        #
  #                                #
  # Computes partial derivative    #
  # of field with respect to ith   #
  # coordinate. Moves field from   #
  # scalar positions to inerfaces. #
  # Assumes Dirichlet BC at bottom #
  #================================#

  partialder_s2i<-function(p,coord,f){

    # p=Coordinate rank, coord=coordinate vector, f=3d field

    # Quickly do 1D calculation
    if (length(f) == length(coord) ) {
       nz = length(f)
       partialder_s2i<-numeric(length(coord))
       # Initialize  zhalo
       zhalo = c( -coord[1],coord, 2*coord[nz]-coord[nz-1] )

       # Construct fhalo in z-direction with f = 0 on bottom boundary
       fhalo<-numeric(nz+2)
       fhalo[2:(nz+1)] = f[1:nz]
       fhalo[1] = -f[1]
    
       # Compute partialder
       for (k in 1:nz) {
           partialder_s2i[k] = 1/(zhalo[k+1] - zhalo[k])*
			 ( fhalo[k+1]-fhalo[k] )             
          }   # k  
       return(partialder_s2i)
       }  # if 1D
       

    # Get dimension sizes
    nx = dim(f)[1]
    ny = dim(f)[2]
    nz = dim(f)[3]
    partialder_s2i<-array(dim=c(nx,ny,nz))

    if (p == 1) {        # ppix

       delta = coord[2] - coord[1]

       # Construct f-halo on x-direction
       fhalo<-array(dim=c(nx+2,ny,nz))
       fhalo[2:(nx+1),1:ny,1:nz] = f
       fhalo[1, , ] = f[nx, , ]
       fhalo[nx+2, , ] = f[1, , ]

       # Compute partialder
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                partialder_s2i[i,j,k] = 1/delta*( fhalo[i+1,j,k] - 
				      	fhalo[i,j,k] )
		# Klugey fudge in indexing between LHS and RHS
             }
           }
	}	  

    } else if (p == 2) {   # ppiy

       delta = coord[2] - coord[1]

       # Construct f-halo on y-direction
       fhalo<-array(dim=c(nx,ny+2,nz))
       fhalo[1:nx,2:(ny+1),1:nz] = f
       fhalo[ ,1, ] = f[ ,ny, ]
       fhalo[ ,ny+2 , ] = f[ ,1 , ]

       # Compute partialder
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                partialder_s2i[i,j,k] = 1/delta*( fhalo[i,j+1,k] - 
				      	fhalo[i,j,k] )
		# Klugey fudge in indexing between LHS and RHS
             }
          }
       }


    } else if (p == 3) {        # ppiz

       # Initialize  zhalo
       zhalo = c( -coord[1],coord, 2*coord[nz]-coord[nz-1] )

       # Construct fhalo in z-direction with f = 0 on bottom boundary
       fhalo<-array(dim=c(nx,ny,nz+2))
       fhalo[ , ,2:(nz+1)] = f[ , ,1:nz ]
       fhalo[ , ,1] = -f[ , ,1]
    
       # Compute partialder
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                partialder_s2i[i,j,k] = 1/(zhalo[k+1] - zhalo[k])*
			 ( fhalo[i,j,k+1]-fhalo[i,j,k] )
             }
          }     
       }
    }

    partialder_s2i
}

 

  #================================#
  # Function partialder_i2s        #
  #                                #
  # Computes partial derivative    #
  # of field with respect to ith   #
  # coordinate. Moves field from   #
  # interface positions to scalar. #
  # Assumes Dirichlet BC at top    #
  #================================#

  partialder_i2s<-function(p,coord,f){

    # p=Coordinate rank, coord=coordinate vector, f=3d field

    # Quickly do 1D calculation
    if (length(f) == length(coord)) {
       nz = length(coord)
       partialder_i2s<-numeric(length(coord))
       # Initialize  zhalo
       zhalo = c( -coord[1],coord, 2*coord[nz]-coord[nz-1] )

       # Construct fhalo in z-direction with f = 0 on bottom boundary
       fhalo<-numeric(nz+2)
       fhalo[2:(nz+1)] = f[1:nz]
       fhalo[nz+2] = 0
    
       # Compute partialder
       for (k in 1:nz) {
           partialder_i2s[k] = 2/(zhalo[k+2] - zhalo[k])*
			 ( fhalo[k+2]-fhalo[k+1] )             
          }   # k  
       return(partialder_i2s)
       }  # if 1D


    # Get dimension sizes
    nx = dim(f)[1]
    ny = dim(f)[2]
    nz = dim(f)[3]
    partialder_i2s<-array(dim=c(nx,ny,nz))

    if (p == 1) {        # ppix

       delta = coord[2] - coord[1]

       # Construct f-halo on x-direction
       fhalo<-array(dim=c(nx+2,ny,nz))
       fhalo[2:(nx+1),1:ny,1:nz] = f
       fhalo[1, , ] = f[nx, , ]
       fhalo[nx+2, , ] = f[1, , ]

       # Compute partialder
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                partialder_i2s[i,j,k] = 1/delta*( fhalo[i+2,j,k] - 
				      	fhalo[i+1,j,k] )
		# Klugey fudge in indexing between LHS and RHS
             }
           }
	}	  

    } else if (p == 2) {   # ppiy

       delta = coord[2] - coord[1]

       # Construct f-halo on y-direction
       fhalo<-array(dim=c(nx,ny+2,nz))
       fhalo[1:nx,2:(ny+1),1:nz] = f
       fhalo[ ,1, ] = f[ ,ny, ]
       fhalo[ ,ny+2 , ] = f[ ,1 , ]

       # Compute partialder
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                partialder_i2s[i,j,k] = 1/delta*( fhalo[i,j+2,k] - 
				      	fhalo[i,j+1,k] )
		# Klugey fudge in indexing between LHS and RHS
             }
          }
       }


    } else if (p == 3) {        # ppiz

       # Initialize  zhalo
       zhalo = c( -coord[1],coord, 2*coord[nz]-coord[nz-1] )

       # Construct fhalo in z-direction with f = 0 on upper boundary
       fhalo<-array(dim=c(nx,ny,nz+2))
       fhalo[ , ,2:(nz+1)] = f[ , ,1:nz]
       fhalo[ , ,nz+2] = 0
    
       # Compute partialder_i2s
       for (i in 1:nx) {
          for (j in 1:ny) {
             for (k in 1:nz) {
                partialder_i2s[i,j,k] = 2/(zhalo[k+2] - zhalo[k])*
			 ( fhalo[i,j,k+2]-fhalo[i,j,k+1] )
             }
          }     
       }
    }

    partialder_i2s
}

  #=================================#
  # Function partialder_i2s4d       #
  #                                 #
  # Computes partial derivative     #
  # of 4d field with respect to ith #
  # coordinate. Moves field from    #
  # interface positions to scalar.  #
  # Assumes Dirichlet BC at top     #
  #=================================#

  partialder_i2s4d<-function(p,coord,f){

    # Get dimension sizes
    nx = dim(f)[1]
    ny = dim(f)[2]
    nz = dim(f)[3]
    nt = dim(f)[4]
    partialder_i2s4d<-array(dim=c(nx,ny,nz,nt))
    for (l in 1:nt){
    	partialder_i2s4d[ , , ,l]=partialder_i2s(p,coord,f[ , , ,l])
    }
    partialder_i2s4d
  }

  #==========================#
  # function ddz2matrix      #
  #                          #
  # For given vertical grid, #
  # computes matrix that     #
  # implements (d/dz)^2 on   #
  # quantity f. Switch bc    #
  # specifies Dirichlet BCs  #
  # with f(0)=f(H)=0 or      # 
  # Neumann with df/dz(0)=   #
  # df/dz(H)=0               #
  #                          #
  #==========================#

  ddz2matrix<-function(z,s_or_i,bc){

      # Check boundary conditions
      if ( (bc != "d") & (bc != "n") & (bc != 'fs')) {
         stop("Boundary conditions must be either 'd', 'n', or 'fs' ")
	 }

      # Check C grid position
      if ( (s_or_i != "s") & (s_or_i !="i") ) {
         stop(" C grid position must be either 's' or 'i' ")
	 }
      
      # Get nz, initialize arrays
      nz = length(z)
      z2halo=numeric(nz+2)
      zhalo=numeric(nz+2)
      A=array(dim=c(nz,nz))
      
      # Initialize A, zhalo, z2halo
      A[ ,]=0 
      zhalo  = c( -z[1],z,2*z[nz]-z[nz-1] )
      zint = zinterp(z)
      z2halo = c( -zint[2],zint,2*z[nz]-zint[nz] )

      # Do scalar calculation
      if (s_or_i == "s"){
         for (k in 1:nz){
            if (k != 1) {
               A[k,k-1] = 2*( (zhalo[k+2]-zhalo[k])*(zhalo[k+1]-zhalo[k]) )^(-1)
             }
             A[k,k] = - 2/(zhalo[k+2]-zhalo[k])*((zhalo[k+2]-zhalo[k+1])^(-1) +
	     	      	 (zhalo[k+1]-zhalo[k])^(-1))
             if (k != nz) {
                A[k,k+1] = 2*( (zhalo[k+2]-zhalo[k])*(zhalo[k+2]-zhalo[k+1]) )^(-1)
             }
          } 
          # Overwrite some values to implement Dirichlet BCs 
          if (bc == "d") {
             A[nz,nz] = -2/(zhalo[nz+2]-zhalo[nz])*
                 ( 2/(zhalo[nz+2]-zhalo[nz+1]) + 1/(zhalo[nz+1]-zhalo[nz]) )
             A[1,1] = -2/(zhalo[3]-zhalo[1])*( 1/(zhalo[3]-zhalo[2]) +
                      2/(zhalo[2]-zhalo[1]) )            
		      }
          # Ditto for Neumann BCs as per my notes, 4/30/14
          else if (bc == "n") {
             A[1,1] = -A[1,2]
             A[nz,nz] = -A[nz,nz-1]
         }
      }

      # Interface calculation  
      else if (s_or_i == "i") {
         for (k in 1:nz){
            if (k != 1) {
               A[k,k-1]<-((zhalo[k+1]-zhalo[k])*(z2halo[k+1]-z2halo[k]) )^(-1)
            }
            A[k,k] <- - 1/( zhalo[k+1]-zhalo[k] )* 
                 ((z2halo[k+2]-z2halo[k+1])**(-1)+(z2halo[k+1]-z2halo[k])^(-1))
            if (k != nz) {
               A[k,k+1]<-( (zhalo[k+1]-zhalo[k])*(z2halo[k+2]-z2halo[k+1]) )^(-1)
             }
         }
        #  BCs
         if ( bc == "d"){ # Dirichlet, top BCs implmntd automatically
            A[1,1] = 0 
            A[2,1] = 0 
	    }
         else if ( bc == "n" ){  # Neumann as per notes 
                                       # 5/2-5/5, 5/15 2014
            A[1,2] = -A[1,1]          
            A[nz,nz] = 1/3*A[nz,nz]
            A[nz,nz-1] = 2/3*A[nz,nz-1]
	    }
         else if (bc == "fs") {
            A[1,1] = 0   # Since d2f/dz2 = 0 at bottom
            A[1,2] = 0
            A[nz,nz-2] = 1/2*A[nz,nz-1] # As per notes 5/14/14
            A[nz,nz-1] = 1/2*A[nz,nz]
            A[nz,nz] = -1/4*A[nz,nz]
	 }
      }
      A
}