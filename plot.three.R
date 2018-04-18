##
## overplot three different genes and fits for the same schema, condition and turn-off
##

source("~/R/plot.bars.R")

plot.three = function(schema="gse70796", condition="GR-REV", gene1="ZPR1", gene2="CPuORF8", gene3="PKS4", turnOff=0.0, ylim=c(.1,10)) {

    ## fixed params
    rhoc0 = 25
    rhon0 = 1
    nu = 10

    t.data = getTimes(schema=schema, condition=condition)/60

    gene1.data = getExpression(schema=schema, condition=condition, gene=gene1)
    gene2.data = getExpression(schema=schema, condition=condition, gene=gene2)
    gene3.data = getExpression(schema=schema, condition=condition, gene=gene3)

    gene1.base = mean(gene1.data[t.data==0])
    gene2.base = mean(gene2.data[t.data==0])
    gene3.base = mean(gene3.data[t.data==0])

    gene1.fit = transmodel.fit(schema=schema, condition=condition, gene=gene1, turnOff=turnOff, doPlot=F)
    gene2.fit = transmodel.fit(schema=schema, condition=condition, gene=gene2, turnOff=turnOff, doPlot=F)
    gene3.fit = transmodel.fit(schema=schema, condition=condition, gene=gene3, turnOff=turnOff, doPlot=F)

    ## hack for RNA-seq ZPR1 fit having negative rhop0
    if (gene1=="ZPR1") {
        gene1.rhop0 = gene1.base*0.9
        gene1.fit = transmodel.fit(schema=schema, condition=condition, gene=gene1, turnOff=turnOff, doPlot=F, rhop0=gene1.rhop0)
        gene1.etap = gene1.fit$estimate[1]
        gene1.gammap = gene1.fit$estimate[2]
    } else {
        gene1.rhop0 = gene1.fit$estimate[1]
        gene1.etap = gene1.fit$estimate[2]
        gene1.gammap = gene1.fit$estimate[3]
    }
    gene1.etap.hat = gene1.etap*rhon0/gene1.rhop0
    gene1.logFCinf = log2(1 + gene1.etap.hat/gene1.gammap*rhoc0/rhon0)

    gene2.rhop0 = gene2.fit$estimate[1]
    gene2.etap = gene2.fit$estimate[2]
    gene2.gammap = gene2.fit$estimate[3]
    gene2.etap.hat = gene2.etap*rhon0/gene2.rhop0
    gene2.logFCinf = log2(1 + gene2.etap.hat/gene2.gammap*rhoc0/rhon0)

    gene3.rhop0 = gene3.fit$estimate[1]
    gene3.etap = gene3.fit$estimate[2]
    gene3.gammap = gene3.fit$estimate[3]
    gene3.etap.hat = gene3.etap*rhon0/gene3.rhop0
    gene3.logFCinf = log2(1 + gene3.etap.hat/gene3.gammap*rhoc0/rhon0)

    ## get R-squared and error metric
    gene1.fitValues = rhop(rhoc=rhoc0, nu=nu, t=t.data, rhop0=gene1.rhop0, etap=gene1.etap, gammap=gene1.gammap, turnOff=turnOff)
    gene1.R2 = Rsquared(gene1.fitValues, gene1.data)
    gene1.error = errorMetric(gene1.fitValues, gene1.data)

    gene2.fitValues = rhop(rhoc=rhoc0, nu=nu, t=t.data, rhop0=gene2.rhop0, etap=gene2.etap, gammap=gene2.gammap, turnOff=turnOff)
    gene2.R2 = Rsquared(gene2.fitValues,gene2.data)
    gene2.error = errorMetric(gene2.fitValues, gene2.data)

    gene3.fitValues = rhop(rhoc=rhoc0, nu=nu, t=t.data, rhop0=gene3.rhop0, etap=gene3.etap, gammap=gene3.gammap, turnOff=turnOff)
    gene3.R2 = Rsquared(gene3.fitValues, gene3.data)
    gene3.error = errorMetric(gene3.fitValues, gene1.data)

    ## the model
    t.model = 0:200/100
    gene1.model = rhop(rhoc=rhoc0,nu=nu, t=t.model, rhop0=gene1.rhop0, etap=gene1.etap, gammap=gene1.gammap, turnOff=turnOff)
    gene2.model = rhop(rhoc=rhoc0,nu=nu, t=t.model, rhop0=gene2.rhop0, etap=gene2.etap, gammap=gene2.gammap, turnOff=turnOff)
    gene3.model = rhop(rhoc=rhoc0,nu=nu, t=t.model, rhop0=gene3.rhop0, etap=gene3.etap, gammap=gene3.gammap, turnOff=turnOff)

    ## the colors
    colors = c("blue","red","darkgreen")

    ## the plot
    plot(t.model, gene1.model/gene1.rhop0, type="l", xlim=c(0,2), xlab="time (h)", ylab="relative transcript level", ylim=ylim, col=colors[1], log="y")
    plot.bars(t.data, gene1.data/gene1.rhop0, over=T, pch=19, cex=1.5, bg=colors[1], col=colors[1])
    
    lines(t.model, gene2.model/gene2.rhop0, col=colors[2])
    plot.bars(t.data, gene2.data/gene2.rhop0, over=T, pch=19, cex=1.5, bg=colors[2], col=colors[2])
    
    lines(t.model, gene3.model/gene3.rhop0, col=colors[3])
    plot.bars(t.data, gene3.data/gene3.rhop0, over=T, pch=19, cex=1.5, bg=colors[3], col=colors[3])
    
    ## the legend
    legend(0, ylim[2], xjust=0, yjust=1, pch=19, pt.cex=1.5, col=colors, text.col=colors, cex=1.2,
           c(gene1,gene2,gene3)
           )

    text(1.5, ylim[2], paste(condition), pos=1, cex=1.5)

    xtext = 2.0
    ystart = (ylim[1]*ylim[2]^3)^0.25
    frac = (ylim[1]/ylim[2])^0.07

    text(xtext, ystart*frac^0, col="black",   pos=2, bquote(paste(hat(eta)[p], "        ",             gamma[p], "      ",              logFC[infinity],"    ",                r^2)))

    text(xtext, ystart*frac^1, col=colors[1], pos=2, bquote(paste(.(round(gene1.etap.hat,3)), "   ", .(round(gene1.gammap,3)), "   ", .(round(gene1.logFCinf,2)), "   ", .(round(gene1.R2,2)) )))
    text(xtext, ystart*frac^2, col=colors[2], pos=2, bquote(paste(.(round(gene2.etap.hat,3)), "   ", .(round(gene2.gammap,3)), "   ", .(round(gene2.logFCinf,2)), "   ", .(round(gene2.R2,2)) )))
    text(xtext, ystart*frac^3, col=colors[3], pos=2, bquote(paste(.(round(gene3.etap.hat,3)), "   ", .(round(gene3.gammap,3)), "   ", .(round(gene3.logFCinf,2)), "   ", .(round(gene3.R2,2)) )))


}
