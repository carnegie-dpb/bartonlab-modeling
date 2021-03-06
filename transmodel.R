source("~/modeling/rhoc.R")
source("~/modeling/rhon.R")
source("~/modeling/rhop.R")

source("~/modeling/Rsquared.R")
source("~/modeling/errorMetric.R")

##
## plot linear transcription model for a direct target in one condition
##

transmodel = function(turnOff=0, rhon0, rhoc0, nu, rhop0, etap, gammap, dataTimes, dataValues, dataLabel=NA, plotBars=FALSE, main="") {

    ## set rhop0 = mean of minimum time points
    if (!hasArg(rhop0) && hasArg(dataValues) && hasArg(dataTimes)) {
        tMin = min(dataTimes)
        rhop0 = mean(dataValues[dataTimes==tMin])
    }

    ## calculation interval
    if (hasArg(dataTimes)) {
        t = seq(0, max(dataTimes), by=0.01)
    } else {
        t = seq(0, 2, by=0.01)
    }

    ## cytoplasmic GR-TF concentration
    rhoc_t = rhoc(t=t, rhoc0=rhoc0, nu=nu)
    
    ## nuclear GR-TF concentration
    rhon_t = rhon(t=t, rhoc0=rhoc0,rhon0=rhon0,nu=nu)
    
    ## transcript concentration
    rhop_t = rhop(t=t, turnOff=turnOff, rhoc0=rhoc0,nu=nu, rhop0=rhop0,etap=etap,gammap=gammap)
    
    ## plot TF concentrations on right axis
    plot(t, rhon_t, type="l", col="blue", axes=FALSE, xlab=NA, ylab=NA, ylim=c(0,max(rhon_t)))
    axis(side=4) 
    par(new=TRUE)

    ## plot transcript concentration on left axis
    if (hasArg(dataValues)) {
        ymin = min(dataValues,rhop_t)
        ymax = max(dataValues,rhop_t)
        plot(t, rhop_t, type="l", xlab="time (h)", ylab="model concentration (lines), mRNA level (points)", col="red", ylim=c(0,ymax), main=main)
    } else {
        ymin = min(rhop_t)
        ymax = max(rhop_t)
        plot(t, rhop_t, type="l", xlab="time (h)", ylab="nuclear concentration", col="red", ylim=c(0,ymax), main=main)
    }

    ## compare with provided data on left axis
    if (hasArg(dataTimes) && hasArg(dataValues)) {
        if (plotBars) {
            ## plot mean and error bars
            for (ti in unique(dataTimes)) {
                y = mean(dataValues[dataTimes==ti])
                sd = sd(dataValues[dataTimes==ti])
                points(ti, y, pch=19, col="red")
                segments(ti, (y-sd), ti, (y+sd), col="red")
            }
        } else {
            ## plot each point
            points(dataTimes, dataValues, pch=19, col="red")
        }
        ## get R-squared and error metric
        fitValues = rhop(t=dataTimes, turnOff=turnOff, rhoc0=rhoc0,nu=nu, rhop0=rhop0,etap=etap,gammap=gammap)
        R2 = Rsquared(fitValues,dataValues)
        error = errorMetric(fitValues,dataValues)
    }

    ## other metrics
    logFCinf = log2(1 + etap/gammap*rhoc0/rhop0)
    kappa = nu*etap*rhoc0/rhop0
    etap.hat = etap*rhon0/rhop0

    ## LIN
    maxRight = ymax*0.8
    step = ymax*0.05
    xtext = par()$usr[2]*0.80
    ylegend = 0

    if (hasArg(dataLabel) && !is.na(dataLabel)) {
        legend(max(t), ylegend, xjust=1, yjust=0, lty=c(1,1,0), pch=c(-1,-1,19), col=c("blue","red","red"),
               c(
                   expression(paste(rho[n]," (right axis)")),
                   expression(paste(rho[p])),
                   dataLabel
                   )
               )
    } else {
        legend(max(t), ylegend, xjust=0, yjust=0, lty=1, col=c("blue","red"),
               c(
                   expression(paste(rho[n],"  (right axis)")),
                   expression(rho[p])
                   )
               )
    }

    text(xtext, maxRight-step*1, bquote(rho[c0]==.(round(rhoc0,1))), pos=3, col="blue")
    text(xtext, maxRight-step*2, bquote(rho[n0]==.(round(rhon0,1))), pos=3, col="blue")
    text(xtext, maxRight-step*3, bquote(paste(nu==.(signif(nu,3))," ",h^-1)), pos=3, col="blue")

    text(xtext, maxRight-step*5, bquote(rho[p0]==.(signif(rhop0,3))), pos=3, col="red")
    text(xtext, maxRight-step*6, bquote(paste(eta[p]==.(signif(etap,3))," ",h^-1)), pos=3, col="red")
    text(xtext, maxRight-step*7, bquote(paste(hat(eta)[p]==.(signif(etap.hat,3))," ",h^-1)), pos=3, col="red")
    text(xtext, maxRight-step*8, bquote(paste(gamma[p]==.(round(gammap,3))," ",h^-1)), pos=3, col="red")

    ## flag suspect fits
    if (rhop0<0) {
        text(xtext+0.2, maxRight-step*5+step*0.2, "!", pos=3, col="red", font=2)
    }
    if (abs(etap.hat)>5) {
        text(xtext+0.2, maxRight-step*7+step*0.2, "!", pos=3, col="red", font=2)
    }
    if (gammap<0.1 || gammap>10) {
        text(xtext+0.2, maxRight-step*8+step*0.2, "!", pos=3, col="red", font=2)
    }

    ## derived fit metrics
    if (hasArg(dataTimes) && hasArg(dataValues)) {
        text(xtext, maxRight-step*10, bquote(paste(kappa==.(signif(kappa,3))," ",h^-2)), pos=3, col="black")
        text(xtext, maxRight-step*11, bquote(logFC(inf)==.(round(logFCinf,2))), pos=3, col="black")
        text(xtext, maxRight-step*12, bquote(r^2==.(round(R2,5))), pos=3, col="black")
    }

    par(new=FALSE)

}
