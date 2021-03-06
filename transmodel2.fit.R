source("~/R/getTimes.R")
source("~/R/getExpression.R")
source("~/modeling/transmodel.fit.R")
source("~/modeling/transmodel2.R")
source("~/modeling/transmodel2.error.R")

##
## use nlm to find the best fit to a given set of parameter guesses and a given time,data series
##

transmodel2.fit = function(host="localhost",
                           fitTerms="rhos0.etas.gammas",
                           schema, condition, gene1, gene2,
                           rhon0=1, rhoc0=25, nu=10,
                           rhop0=1, etap=1, gammap=4,
                           rhos0=1, etas=1, gammas=4,
                           dataTimes, data1Values, data2Values, 
                           data1Label=NA, data2Label=NA, plotBars=FALSE, doPlot=TRUE, main=""
                           ) {

    ## default plot title if not given
    if (doPlot && main=="") {
        main = paste(schema,condition,":",gene1,"->",gene2)
    }

    ## if gammas supplied, leave it fixed in fits
    if (hasArg(gammas)) fitTerms = "rhos0.etas"

    ## if rhos0 supplied, leave it fixed in fits
    if (hasArg(rhos0)) fitTerms = "etas.gammas"
    
    ## get time (in hours) and expression arrays for the given schema and gene IDs
    if (!hasArg(dataTimes)) {
        dataTimes = getTimes(schema=schema, condition=condition, host=host)
        if (max(dataTimes)>5) dataTimes = dataTimes/60
        data1Values = getExpression(schema=schema, condition=condition, gene=toupper(gene1), host=host)
        data2Values = getExpression(schema=schema, condition=condition, gene=toupper(gene2), host=host)
        if (is.na(data1Label)) data1Label = paste(condition,gene1,sep=":")
        if (is.na(data2Label)) data2Label = paste(condition,gene2,sep=":")
    }

    ## estimate rhos0 from minimum t data points
    tMin = min(dataTimes)
    if (!hasArg(rhos0)) {
        rhos0 = mean(data2Values[dataTimes==tMin])
    }
    
    ## do minimization of primary transcript params
    fit1 = transmodel.fit(turnOff=0, rhoc0=rhoc0, rhon=rhon0, nu=nu,
                          dataTimes=dataTimes, dataValues=data1Values, doPlot=FALSE)
    
    ## new params from fit
    rhop0 = fit1$estimate[1]
    etap = fit1$estimate[2]
    gammap = fit1$estimate[3]

    if (fit1$code==4) print("*** fit1 ITERATION LIMIT EXCEEDED ***")
    if (fit1$code==5) print("*** fit1 MAXIMUM STEP SIZE EXCEEDED FIVE CONSECUTIVE TIMES ***")

    ## do minimization of secondary transcript params
    if (fitTerms=="rhos0.etas.gammas") {

        fit2 = nlm(p=c(rhos0,etas,gammas), fitTerms=fitTerms,
                   rhoc0=rhoc0, rhon0=rhon0, nu=nu, etap=etap, gammap=gammap, f=transmodel2.error, gradtol=1e-5, iterlim=1000, dataTimes=dataTimes, data2Values=data2Values)
        rhos0 = fit2$estimate[1]
        etas = fit2$estimate[2]
        gammas = fit2$estimate[3]
        
    } else if (fitTerms=="rhos0.etas") {

        fit2 = nlm(p=c(rhos0,etas), fitTerms=fitTerms, gammas=gammas,
                   rhoc0=rhoc0, rhon0=rhon0, nu=nu, etap=etap, gammap=gammap, f=transmodel2.error, gradtol=1e-5, iterlim=1000, dataTimes=dataTimes, data2Values=data2Values)
        rhos0 = fit2$estimate[1]
        etas = fit2$estimate[2]
        
    } else if (fitTerms=="etas.gammas") {
        
        fit2 = nlm(p=c(etas,gammas), fitTerms=fitTerms, rhos0=rhos0,
                   rhoc0=rhoc0, rhon0=rhon0, nu=nu, etap=etap, gammap=gammap, f=transmodel2.error, gradtol=1e-5, iterlim=1000, dataTimes=dataTimes, data2Values=data2Values)
        etas = fit2$estimate[1]
        gammas = fit2$estimate[2]

    }
    
    if (fit2$code==4) print("*** fit2 ITERATION LIMIT EXCEEDED ***")
    if (fit2$code==5) print("*** fit2 MAXIMUM STEP SIZE EXCEEDED FIVE CONSECUTIVE TIMES ***")

    ## plot it
    if (doPlot) {
        transmodel2(turnOff=0,
                    rhoc0=rhoc0, rhon0=rhon0, nu=nu,
                    rhop0=rhop0, etap=etap, gammap=gammap,
                    rhos0=rhos0, etas=etas, gammas=gammas,
                    dataTimes=dataTimes, data1Values=data1Values,data1Label=data1Label, data2Values=data2Values,data2Label=data2Label, plotBars=plotBars, main=main)
    }

    return(fit2)

}

