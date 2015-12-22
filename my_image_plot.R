#=============================#
# David's image.plot function #
# plus my comments and mods   #
# 12/15                       #
#                             #
# Follows original image.plot #s
# code as found at URL below  #
#=============================#

#http://www.image.ucar.edu/~nychka/Fields/Source/image.plot.R

my.image.plot <- function(x,y,z,...,				
	main = "",
	ylim = range(y),
	zlim = range(z),
	xlab = "",
	ylab = "",
	cex.lab = 1,
	cex.axis = 1,
	cex.main = 1,
	cex.legend = 1,
    add = FALSE, 
    nlevel = 64, 
    horizontal = FALSE, 
    legend.shrink = 0.9, 
    legend.width = 1.2, 
    legend.mar = ifelse(horizontal, 3.1, 5.1), # set legend margin for 												   					horizontal or vertical legend
    legend.lab = NULL, 
    graphics.reset = FALSE, 
    bigplot = NULL,   # Plotting region coords for image 
    smallplot = NULL, # plotting region coords for legend
    legend.only = FALSE, 
    col = tim.colors(nlevel), lab.breaks = NULL, axis.args = NULL, legend.args = NULL, 
    font.main = 2, 
    midpoint = FALSE, 
    border = NA, lwd = 1) {

    #================#
    # Begin function #
    #================#

    old.par <- par(no.readonly = TRUE)  # old parameters which can be changed
    info <- imageplot.info(x,y,z)   # returns xlim,ylim,zlim
    if (add) {
        big.plot <- old.par$plt  # Keep old coordinates of image region 
		    		    		 # as fraction of figure region (for adding legend)
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    # temp produces $smallplot, $bigplot, plotting regions for image and legend
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)			   # Get size of image region
        }

        #=======================#
        # Plot image, no legend #
		#=======================#

 		# plot image-only as png, save as temp.png
        library(png)
        png('temp.png')  # deleted 'cairo' type -- doesn't work!
        par(xaxs='i',yaxs='i',mar=c(0,0,0,0))  # Style for axis finding, no margin
        image(x,y,z,axes=FALSE,xlab='',ylab='',col=col,zlim=zlim) # Plot image!
        dev.off()
        # Plot axes, ticks, and labels for image
        image(x,y,array(NA,dim=c(length(x),length(y))),
        	 xlab=xlab,ylab=ylab,main=main,ylim=ylim,
        	 cex.lab=cex.lab,cex.axis = cex.axis, 
        	 cex.main=cex.main, font.main=font.main)
		# Add image to above skeleton
        rasterImage(readPNG('temp.png'),min(x),min(y),max(x),max(y))
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
	minz <- zlim[1]
	maxz <- zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!is.null(breaks) & !is.null(lab.breaks)) {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
            at = breaks, labels = lab.breaks), cex.axis=cex.legend)
    }
    else {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
            cex.axis=cex.legend)
    }
    if (!horizontal) {
        if (is.null(breaks)) {
            # ----- Start of my modifications -----
            library(png)
            png('temp2.png',width=10,height=nlevel-1)
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
               ylab = "", col = col)
            dev.off()
            plot(0:1,c(minz,maxz),col='transparent',
            	xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
            rasterImage(readPNG('temp2.png'),0,minz,1,maxz)
            #image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            #    ylab = "", col = col)
            # ------ End of my modifications ------
        }
        else {
            # ----- Start of my modifications -----
            library(png)
            png('temp2.png',width=10,height=nlevel-1)
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
               ylab = "", col = col, breaks = breaks)
            dev.off()
            plot(0:1,c(minz,maxz),col='transparent',
            	xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
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
            png('temp2.png',height=10,width=nlevel-1)
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
               ylab = "", col = col)
            dev.off()
            plot(c(minz,maxz),0:1,col='transparent',
            	xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
            rasterImage(readPNG('temp2.png'),minz,0,maxz,1)
            #image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            #    ylab = "", col = col)
            # ------ End of my modifications ------
        }
        else {
            # ----- Start of my modifications -----
            library(png)
            png('temp2.png',height=10,width=nlevel-1)
            par(xaxs='i',yaxs='i',mar=c(0,0,0,0))
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
               ylab = "", col = col, breaks = breaks)
            dev.off()
            plot(c(minz,maxz),0:1,col='transparent',
            	xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
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

