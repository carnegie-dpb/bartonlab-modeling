##
## plot both RNA-seq and microarray results for a given gene
##

plot.both = function(condition="GR-REV", gene="ZPR1", turnOff=0.0, ylim=c(.1,10)) {

    rhoc0 = 25
    rhon0 = 1
    nu = 10

    gse30703.t = getTimes(schema="gse30703",condition=condition)/60
    gse70796.t = getTimes(schema="gse70796",condition=condition)/60

    gse30703 = getExpression(schema="gse30703",condition=condition,gene=gene)
    gse70796 = getExpression(schema="gse70796",condition=condition,gene=gene)

    gse30703.base = mean(gse30703[gse30703.t==0])
    gse70796.base = mean(gse70796[gse70796.t==0])

    gse30703.fit = transmodel.fit(schema="gse30703", condition=condition, gene=gene, turnOff=turnOff, doPlot=F)
    gse70796.fit = transmodel.fit(schema="gse70796", condition=condition, gene=gene, turnOff=turnOff, doPlot=F)

    gse30703.rhop0 = gse30703.fit$estimate[1]
    gse30703.etap = gse30703.fit$estimate[2]
    gse30703.gammap = gse30703.fit$estimate[3]
    gse30703.etap.hat = gse30703.etap*rhon0/gse30703.rhop0

    gse70796.rhop0 = gse70796.fit$estimate[1]
    gse70796.etap = gse70796.fit$estimate[2]
    gse70796.gammap = gse70796.fit$estimate[3]
    gse70796.etap.hat = gse70796.etap*rhon0/gse70796.rhop0

    ## hack for RNA-seq ZPR1 fit having negative rhop0
    if (gene=="ZPR1") {
        gse70796.rhop0 = gse70796.base*0.9
        gse70796.fit = transmodel.fit(schema="gse70796", condition=condition, gene=gene, turnOff=turnOff, doPlot=F, fitTerms="etap.gammap", rhop0=gse70796.rhop0)
        gse70796.etap = gse70796.fit$estimate[1]
        gse70796.gammap = gse70796.fit$estimate[2]
        gse70796.etap.hat = gse70796.etap*rhon0/gse70796.rhop0
    }

    ## get R-squared and error metric
    gse30703.fitValues = rhop(rhoc=rhoc0, nu=nu, t=gse30703.t, rhop0=gse30703.rhop0, etap=gse30703.etap, gammap=gse30703.gammap, turnOff=turnOff)
    gse30703.R2 = Rsquared(gse30703.fitValues,gse30703)
    gse30703.error = errorMetric(gse30703.fitValues,gse30703)

    gse70796.fitValues = rhop(rhoc=rhoc0, nu=nu, t=gse70796.t, rhop0=gse70796.rhop0, etap=gse70796.etap, gammap=gse70796.gammap, turnOff=turnOff)
    gse70796.R2 = Rsquared(gse70796.fitValues,gse70796)
    gse70796.error = errorMetric(gse70796.fitValues,gse30703)

    ## the model
    t = 0:200/100
    gse30703.model = rhop(rhoc=rhoc0,nu=nu, t=t, rhop0=gse30703.rhop0, etap=gse30703.etap, gammap=gse30703.gammap, turnOff=turnOff)
    gse70796.model = rhop(rhoc=rhoc0,nu=nu, t=t, rhop0=gse70796.rhop0, etap=gse70796.etap, gammap=gse70796.gammap, turnOff=turnOff)

    ## the plot
    plot(t, gse30703.model/gse30703.rhop0, type="l", xlim=c(0,2), xlab="time (h)", ylab="relative transcript level", ylim=ylim, col="blue", log="y")
    plot.bars(gse30703.t, gse30703/gse30703.rhop0, over=T, pch=19, cex=1.5, bg="blue", col="blue")
    lines(t, gse70796.model/gse70796.rhop0, col="red")
    plot.bars(gse70796.t, gse70796/gse70796.rhop0, over=T, pch=19, cex=1.5, bg="red", col="red")

    legend(0, ylim[2], xjust=0, yjust=1, pch=19, pt.cex=1.5, col=c("blue","red"), text.col=c("blue","red"), cex=1.2,
           c("microarray","RNA-seq")
           )

    text(1.5, ylim[2], paste(condition,":",gene), pos=1, cex=1.5)

    xtext = 1.6
    ystart = sqrt(ylim[1]*ylim[2])
    frac = (ylim[1]/ylim[2])^0.08
    text(xtext, ystart*frac^0, col="blue", pos=3, bquote(paste("microarray: ",hat(eta)[p]==.(round(gse30703.etap.hat,2))," ",h^-1)))
    text(xtext, ystart*frac^1, col="blue", pos=3, bquote(paste("                    ",gamma[p]==.(round(gse30703.gammap,2))," ",h^-1)))
    text(xtext, ystart*frac^2, col="blue", pos=3, bquote(paste("              ",r^2==.(round(gse30703.R2,2)))))

    text(xtext, ystart*frac^4, col="red", pos=3, bquote(paste("   RNA-seq: ",hat(eta)[p]==.(round(gse70796.etap.hat,2))," ",h^-1)))
    text(xtext, ystart*frac^5, col="red", pos=3, bquote(paste("                     ",gamma[p]==.(round(gse70796.gammap,2))," ",h^-1)))
    text(xtext, ystart*frac^6, col="red", pos=3, bquote(paste("                ",r^2==.(round(gse70796.R2,2)))))

    


}
