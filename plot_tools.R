

#---------------------------#
# Plot set-up functions     #
#---------------------------#
png1x2<-function(filename) {
        png(paste(filename,".png",sep=""),width=12,height=5,units="in",res=72)
	par(mfrow=c(1,2), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
}
pdf1x1<-function(filename) {
        pdf(filename,width=8,height=6,bg="white")
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
pdf4x3<-function(filename) {
        pdf(filename,width=12,height=12)
	par(mfrow=c(4,3), mar=c(5,6,4,7), oma=c(0.5,0.5,0.5,3))
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

plot_zprofile<-function(x,xlab,main,xlim="NULL",color="NULL",zvec=1:nz){
	length(x)->nz
	if (xlim == "NULL") {
	   xlim=range(x[zvec])
	   }
	if (color == "NULL") {
	   color="black"
	   }
	plot(x[zvec],1e-3*z[zvec], ylab="z (km)",main= main, type ="l",
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

plot_tprofile<-function(x,xlab,main,xlim="NULL",ylim="NULL",log="",
			color="NULL",zvec=1:nz){
        if (xlim == "NULL") {
           xlim=range(x[zvec])
           }
        if (ylim == "NULL") {
           ylim=rev(range(tabs[zvec]))
           }
        if (color == "NULL") {
           color="black"
           }
        plot(x[zvec],tabs[zvec], ylab="Temperature [K]",
	         main= main, type ="l",
		 log = log,
                 cex.axis = 1.35,
                 cex.lab = 1.5,
                 cex.main = 1.5,
                 lwd  = 2, 
                 xlab = xlab,
                 xlim = xlim,
		 ylim = ylim,
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


