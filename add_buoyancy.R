#=============#
# Tool to add #
# Archimedean #
# buoyancy to #
# DAM snapshot#
#=============#

ncpath = commandArgs(trailingOnly=TRUE)[1]
library(ncdf)
source("~/Dropbox/Rtools/thermo_tools.R")

nc  = open.ncdf(ncpath,write=TRUE)
rho = get.var.ncdf(nc,"rho")
B   = buoyancy_func(rho,far_field=FALSE)

vars        <- list()
vars[['B']] <- list()
vars[['B']]$longname <- 'Archimedean buoyancy'
vars[['B']]$units    <- 'm/s^2'
vars[['B']]$data     <- B

dims = nc[[8]]
xdim = dims[["x"]]
ydim = dims[["y"]]
zdim = dims[["z"]]
tdim = dims[["time"]]
   
# Make vardef (general, to include RH later)
vardef <- list()
for (name in names(vars)) {
    vardef[[name]] <- var.def.ncdf(name,vars[[name]]$units, 
    		      		     list(xdim,ydim,zdim,tdim), missval=-999,
				     longname=vars[[name]]$longname)
   }

# Fill the NetCDF with data
for (name in names(vars)) {
   var.add.ncdf(nc,vardef[[name]]) 
   close.ncdf(nc)
   nc = open.ncdf(ncpath,write=TRUE)
   put.var.ncdf(nc,name,vars[[name]]$data)
   }

close.ncdf(nc)
