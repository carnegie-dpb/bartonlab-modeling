source("~/R/getExpression.R")
source("~/R/getTimes.R")
source("transmodel.R")
source("transmodel.error.R")

##
## use nlm to find the best fit to a given set of parameter guesses and a given time,data series
##
## set nuFixed=TRUE to hold nu fixed, i.e. not be a fitted parameter

transmodel.fit = function(
                          host="localhost", 
                          fitTerms="rhop0.etap.gammap", turnOff=0,
                          rhoc0=25, rhon0=1, nu=10, rhop0=1, etap=1, gammap=4,
                          schema, gene, condition,
                          dataTimes, dataValues, dataLabel=NA,
                          main="", plotBars=FALSE,  doPlot=TRUE
                          ) {

    ## if gammap supplied, leave it fixed in fits
    if (hasArg(gammap)) fitTerms = "rhop0.etap"

    ## get time (in hours) and expression arrays for the given schema and gene ID from the database
    if (!hasArg(dataTimes)) {
        dataTimes = getTimes(schema=schema, condition=condition, host=host)
        if (max(dataTimes)>5) dataTimes = dataTimes/60
        dataValues = getExpression(schema=schema, condition=condition, gene=toupper(gene), host=host)
        if (is.null(dataValues)) {
            print("No data - aborting.")
            return(NULL)
        }
        if (is.na(dataLabel)) dataLabel = paste(condition,":",gene,sep="")
    }
    
    ## estimate rhop0 from minimum t data points
    tMin = min(dataTimes)
    if (!hasArg(rhop0)) {
        rhop0 = mean(dataValues[dataTimes==tMin])
    }

    ## estimate etap from rhop0, assume up-regulated
    if (!hasArg(etap)) etap = 0.1*rhop0

    if (fitTerms=="rhop0") {

        ## fix etap and gammap, just adjust rhop0 (offset)
        fit = try(nlm(p=c(rhop0), f=transmodel.error, gradtol=1e-5, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff,
            rhoc0=rhoc0,nu=nu,etap=etap, gammap=gammap, dataTimes=dataTimes, dataValues=dataValues))
        if (class(fit)=="try-error") return(NULL)
        ## new params from fit
        rhop0 = fit$estimate[1]
        
    } else if (fitTerms=="etap.gammap") {
        
        ## two-param fit with rhop0 fixed, typically mean of t=0 points set above
        fit = try(nlm(p=c(etap,gammap), f=transmodel.error, gradtol=1e-5, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff,
            rhoc0=rhoc0,nu=nu,rhop0=rhop0, dataTimes=dataTimes, dataValues=dataValues))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, gradtol=1e-5, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, rhoc0=rhoc0, nu=nu, rhop0=rhop0, dataTimes=dataTimes, dataValues=dataValues))
            if (class(fit)=="try-error") return(NULL)
        }
        ## new params from fit
        etap = fit$estimate[1]
        gammap = fit$estimate[2]
        
    } else {

        ## two-param fit with gammap fixed; used to refine rhop0 and etap estimates below as well
        fit = try(nlm(p=c(rhop0,etap), f=transmodel.error, gradtol=1e-5, iterlim=1000, fitTerms="rhop0.etap", turnOff=turnOff,
            rhoc0=rhoc0, nu=nu, gammap=gammap, dataTimes=dataTimes, dataValues=dataValues))
        if (class(fit)=="try-error") return(NULL)
        rhop0 = fit$estimate[1]
        etap = fit$estimate[2]
        
    }

    if (fitTerms=="rhop0.etap.gammap") {

        ## three parameter fit: rhop0, etap, gammap
        fit = try(nlm(p=c(rhop0,etap,gammap), f=transmodel.error, gradtol=1e-5, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff,
            rhoc0=rhoc0, nu=nu, dataTimes=dataTimes, dataValues=dataValues))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, gradtol=1e-5, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, rhoc0=rhoc0, nu=nu, dataTimes=dataTimes, dataValues=dataValues))
            if (class(fit)=="try-error") return(NULL)
        }
        ## new params from fit
        rhop0 = fit$estimate[1]
        etap = fit$estimate[2]
        gammap = fit$estimate[3]

    } else if (fitTerms=="nu.rhop0.etap.gammap") {

        ## fourparameter fit: nu, rhop0, etap, gammap
        fit = try(nlm(p=c(nu,rhop0,etap,gammap), f=transmodel.error, gradtol=1e-5, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff,
            rhoc0=rhoc0, dataTimes=dataTimes, dataValues=dataValues))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, gradtol=1e-5, iterlim=1000, fitTerms=fitTerms, rhoc0=rhoc0, dataTimes=dataTimes, dataValues=dataValues))
            if (class(fit)=="try-error") return(NULL)
        }
        ## new params from fit
        nu = fit$estimate[1]
        rhop0 = fit$estimate[2]
        etap = fit$estimate[3]
        gammap = fit$estimate[4]
        
    }

    ## get R-squared and error metric
    fitValues = rhop(t=dataTimes, turnOff=turnOff, rhoc0=rhoc0,nu=nu, rhop0=rhop0,etap=etap,gammap=gammap)

    R2 = Rsquared(fitValues,dataValues)
    fit$R2 = R2
    
    ## plot it
    if (doPlot) {
        transmodel(turnOff=turnOff,
                   rhoc0=rhoc0, rhon0=rhon0, nu=nu,
                   rhop0=rhop0, etap=etap, gammap=gammap,
                   dataTimes=dataTimes, dataValues=dataValues, dataLabel=dataLabel,
                   main=main, plotBars=plotBars)
    }

    ## return fit in case we want it
    return(fit)

}

