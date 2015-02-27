####################
# Open snapshot.nc #
# and retrieve     #
# dimensions       #
####################

source("~/domain/Rtools/calculus_tools.R")
open.ncdf(snpshtpath)->snpsht_nc

get.var.ncdf(snpsht_nc,"x")->x
get.var.ncdf(snpsht_nc,"y")->y
get.var.ncdf(snpsht_nc,"z")->z
try(get.var.ncdf(snpsht_nc,"time"))->time

nx<-length(x)
ny<-length(y)
nz<-length(z)
zint<-zinterp(z)
if (class(time) == "try-error")  {
   } else {
   nt<-length(time)
   dt = time[2]-time[1]
   }
dx = x[2]-x[1]
dy = y[2]-y[1]
dz = z[2]-z[1]

