
#=============================#
# David's image.plot function #
#=============================#

my.image.plot <- function(x,y,z,..., add = FALSE, nlevel = 64, horizontal = FALSE, 
    legend.shrink = 0.9, legend.width = 1.2, legend.mar = ifelse(horizontal, 
        3.1, 5.1), legend.lab = NULL, graphics.reset = FALSE, 
    bigplot = NULL, smallplot = NULL, legend.only = FALSE, col = tim.colors(nlevel), 
    lab.breaks = NULL, axis.args = NULL, legend.args = NULL, zlim = NULL, 
    font.main = 2,
    midpoint = FALSE, border = NA, lwd = 1) 
{

    old.par <- par(no.readonly = TRUE)
    info <- image.plot.info(x,y,z,...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            # ----- Start of my modifications -----
            # Draw plot, no legend
            library(png)
            png('temp.png',type='cairo')
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(x,y,z,axes=FALSE,xlab='',ylab='',col=col,zlim=zlim)
            dev.off()
#            plot(range(x),range(y),col='transparent',xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
            image(x,y,array(NA,dim=c(length(x),length(y))),...,font.main=font.main)
            rasterImage(readPNG('temp.png'),min(x),min(y),max(x),max(y))
            #image(..., add = add, col = col)
            # ------ End of my modifications ------



        }
        else {
            poly.image(..., add = add, col = col, midpoint = midpoint, 
                border = border, lwd.poly = lwd)
        }
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    if (is.null(zlim)) {
	    minz <- info$zlim[1]
    	maxz <- info$zlim[2]
		} else {
		minz <- zlim[1]
		maxz <- zlim[2]
		}
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!is.null(breaks) & !is.null(lab.breaks)) {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
            axis.args)
    }
    if (!horizontal) {
        if (is.null(breaks)) {
            # ----- Start of my modifications -----
            library(png)
            png('temp2.png',type='cairo',width=10,height=nlevel-1)
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
               ylab = "", col = col)
            dev.off()
            plot(0:1,c(minz,maxz),col='transparent',xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
            rasterImage(readPNG('temp2.png'),0,minz,1,maxz)
            #image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            #    ylab = "", col = col)
            # ------ End of my modifications ------
        }
        else {
            # ----- Start of my modifications -----
            library(png)
            png('temp2.png',type='cairo',width=10,height=nlevel-1)
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
               ylab = "", col = col, breaks = breaks)
            dev.off()
            plot(0:1,c(minz,maxz),col='transparent',xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
            rasterImage(readPNG('temp2.png'),0,minz,1,maxz)
            #image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            #    ylab = "", col = col, breaks = breaks)
            # ------ End of my modifications ------
        }
    }
    else {
        if (is.null(breaks)) {
            # ----- Start of my modifications -----
            library(png)
            png('temp2.png',type='cairo',height=10,width=nlevel-1)
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
               ylab = "", col = col)
            dev.off()
            plot(c(minz,maxz),0:1,col='transparent',xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
            rasterImage(readPNG('temp2.png'),minz,0,maxz,1)
            #image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            #    ylab = "", col = col)
            # ------ End of my modifications ------
        }
        else {
            # ----- Start of my modifications -----
            library(png)
            png('temp2.png',type='cairo',height=10,width=nlevel-1)
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
               ylab = "", col = col, breaks = breaks)
            dev.off()
            plot(c(minz,maxz),0:1,col='transparent',xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
            rasterImage(readPNG('temp2.png'),minz,0,maxz,1)
            #image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            #    ylab = "", col = col, breaks = breaks)
            # ------ End of my modifications ------
        }
    }
    do.call("axis", axis.args)
    box()
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
            1, 4), line = legend.mar - 2)
    }
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}

#---------------------------#
# Plot set-up functions     #
#---------------------------#
png1x2<-function(filename) {
        png(paste(filename,".png",sep=""),width=12,height=5,units="in",res=72)
	par(mfrow=c(1,2), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf1x1<-function(filename) {
        pdf(filename,width=8,height=6)
	par(mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,1))
}

pdf1x2<-function(filename) {
        pdf(filename,width=10,height=4)
	par(mfrow=c(1,2), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf2x1<-function(filename) {
        pdf(filename,width=8,height=8)
	par(mfrow=c(2,1), mar=c(5,5,5,7), oma=c(0.5,0.5,0.5,3))
}


pdf2x2<-function(filename) {
        pdf(filename,width=10,height=8)
	par(mfrow=c(2,2), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf2x3<-function(filename) {
        pdf(filename,width=12,height=8)
	par(mfrow=c(2,3), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf2x4<-function(filename) {
        pdf(filename,width=15,height=8)
	par(mfrow=c(2,4), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf3x3<-function(filename) {
        pdf(filename,width=12,height=10)
	par(mfrow=c(3,3), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}

pdf3x4<-function(filename) {
        pdf(filename,width=15,height=10)
	par(mfrow=c(3,4), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf4x4<-function(filename) {
        pdf(filename,width=15,height=12)
	par(mfrow=c(4,4), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf5x4<-function(filename) {
        pdf(filename,width=15,height=15)
	par(mfrow=c(5,4), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf6x4<-function(filename) {
        pdf(filename,width=15,height=18)
	par(mfrow=c(6,4), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}



#---------------------------#
#  Basic plotting functions #
#---------------------------#


# Function to make stretched color bar
make_colvec <- function(field,center_value,outer_frac = 0.25){
            fieldlim = c(min(field),max(field))
            x0 = abs(center_value-fieldlim[1])/(fieldlim[2]-fieldlim[1])
            xlow = x0-outer_frac
            xhigh= x0+outer_frac
            nmiddle=30
            nouter=5
            x_col=c( xlow/nouter*(0:nouter),
                      xlow+(xhigh-xlow)/nmiddle*(1:nmiddle),
                      xhigh+ 1:nouter*(1-xhigh)/nouter)
            ncol=length(x_col)
            colvec=designer.colors(64,tim.colors(ncol),x=x_col)
            colvec
            }



plot_linegraph<-function(x,y,xlab,ylab,main,col="black",ylim=c(min(y),max(y))){
	par(mar=c(5,6,5,5))
	plot(x,y,xlab=xlab,ylab=ylab,main=main, type="l",
		 col=col,
		 ylim = ylim,
		 cex.axis = 1.5,
		 cex.lab = 1.5,
		 cex.main = 1.5,
		 lwd = 2)
	}

plot_ts<-function(ts,time, main){
	plot(time, ts, xlab="Time [days]",main= main, type ="l",
		 cex.axis = 1.5,
		 cex.lab = 1.5,
		 cex.main = 1.5,
		 ylab = "",
		 lab = labvec
		 )
	} 

plot_xdist<-function(y,ylab, main){
	length(y)->nx
	plot(1e-3*x[1:nx], y, xlab="x (km)",main= main, type ="l",
		 cex.axis = 1.5,
		 cex.lab = 1.5,
		 cex.main = 1.5,
		 ylab = ylab
		 #lab = labvec
		 )
	} 

plot_zprofile<-function(x,xlab,main,xlim="NULL",color="NULL"){
	length(x)->nz
	if (xlim == "NULL") {
	   xlim=c(min(x),max(x))
	   }
	if (color == "NULL") {
	   color="black"
	   }
	plot(x,1e-3*z[1:nz], ylab="z (km)",main= main, type ="l",
		 cex.axis = 1.35,
		 cex.lab = 1.5,
		 cex.main = 1.5,
                 lwd  = 2, 
		 xlab = xlab,
		 xlim = xlim,
		 col = color
		 )
	} 

plot_p_profile<-function(x,xlab,main,xlim,color,zvec=1:nz){
	p_mod=1e-2*p[zvec]
	if (xlim == "NULL") {
	   xlim=range(x[zvec])
	   }
	if (color == "NULL") {
	   color="black"
	   }
	plot(x[zvec],p_mod, ylab="p [hPa]",main= main, type ="l",
		 cex.axis = 1.35,
		 cex.lab = 1.5,
		 cex.main = 1.5,
                 lwd  = 2,
		 xlab = xlab,
		 xlim = xlim,
		 ylim = rev(range(p_mod)),
		 col = color
		 )
	} 

plot_tprofile<-function(x,xlab,main,xlim="NULL",color="NULL",zvec=1:nz){
        if (xlim == "NULL") {
           xlim=range(x[zvec])
           }
        if (color == "NULL") {
           color="black"
           }
        plot(x[zvec],tabs[zvec], ylab="Temperature [K]",
	         main= main, type ="l",
                 cex.axis = 1.35,
                 cex.lab = 1.5,
                 cex.main = 1.5,
                 lwd  = 2, 
                 xlab = xlab,
                 xlim = xlim,
		 ylim = rev(range(tabs[zvec])),
                 col = color
                 )
        }

plot_anom<-function(var_anom,title){
        image.plot(time[tindvec],1e-3*z[1:zplotind],
                var_anom[tindvec,1:zplotind], 
                ylab = "z (km)", xlab = "Time [days]", 
                main = title, lab= c(10,10,7),
                cex.axis = 1.5,
                cex.lab = 1.5,
                cex.main = 1.5,
                axis.args=list(cex.axis=1.5))
        }


plot_xyfield<-function(field,title,zlimits){
	if (zlimits == "NULL") {
	   zlimits=c(min(field),max(field))
	   }
        image.plot(1e-3*x,1e-3*y,field,
                ylab = "y (km)", xlab = "x (km)", 
                main = title, lab= c(10,10,7),
                cex.axis = 1.5,
                cex.lab = 1.5,
                cex.main = 1.5,
                axis.args=list(cex.axis=1.5),
		zlim= zlimits)
        }
plot_xzfield<-function(field,title,zlimits,col=tim.colors(64)){
	if (zlimits == "NULL") {
	   zlimits=c(min(field),max(field))
	   }
	dim(field)[1]->nx
	dim(field)[2]->nz
        my.image.plot(1e-3*x[1:nx],z[1:nz],field,
                ylab = "z [m]", xlab = "x (km)", 
                main = title, lab= c(10,10,7),
                cex.axis = 1.5,
                cex.lab = 1.5,
                cex.main = 1.5,
                axis.args=list(cex.axis=1.5),
				zlim= zlimits,
				col = col)
        }


plot_rzfield_zlim<-function(field,title,zlimits){
	length(field[1,])->nz
	length(field[,1])->nx
        image.plot(1e-3*(x[(nx/2):nx]-x[nx/2]),z[1:nz],
		field[nx/2:nx,],
                ylab = "z [m]", xlab = "r (km)", 
                main = title, lab= c(10,10,7),
                cex.axis = 1.5,
                cex.lab = 1.5,
                cex.main = 1.5,
                axis.args=list(cex.axis=1.5),
		zlim= zlimits)
        }
plot_rzfield<-function(field,title){
	dim(field)[1]->nx
	dim(field)[2]->nz
        image.plot(1e-3*(x[(nx/2):nx]-x[nx/2]),z[1:nz],
		field[(nx/2):nx,],
                ylab = "z [m]", xlab = "r (km)", 
                main = title, lab= c(10,10,7),
                cex.axis = 1.5,
                cex.lab = 1.5,
                cex.main = 1.5,
                axis.args=list(cex.axis=1.5))
        }


