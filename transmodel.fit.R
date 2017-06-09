source("~/R/getExpression.R")
source("~/R/getTimes.R")

source("~/modeling/transmodel.R")
source("~/modeling/transmodel.error.R")

##
## use nlm to find the best fit to a given set of parameter guesses and a given time,data series
##
## set nuFixed=TRUE to hold nu fixed, i.e. not be a fitted parameter

transmodel.fit = function(host="localhost", 
                          fitTerms="rhop0.etap.gammap", turnOff=0,
                          rhoc0=25, rhon0=1, nu=10, rhop0=1, etap=1, gammap=4,
                          schema, gene, condition,
                          dataTimes, dataValues, dataLabel=NA,
                          main="", plotBars=FALSE,  doPlot=TRUE) {

    ## default fitTerms if certain args are supplied
    if (!hasArg(fitTerms)) {
        if (hasArg(rhop0) && hasArg(gammap)) {
            fitTerms = "etap"
        } else if (hasArg(rhop0) && hasArg(etap)) {
            fitTerms = "gammap"
        } else if (hasArg(gammap)) {
            fitTerms = "rhop0.etap"
        } else if (hasArg(rhop0)) {
            fitTerms = "etap.gammap"
        } else if (hasArg(etap)) {
            fitTerms = "rhop0.gammap"
        }
    }

    ## get time (in hours) and expression arrays for the given schema and gene from the database
    if (!hasArg(dataTimes)) {
        dataTimes = getTimes(schema=schema, condition=condition, host=host)
        if (max(dataTimes)>5) dataTimes = dataTimes/60
        dataValues = getExpression(schema=schema, condition=condition, gene=gene, host=host)
        if (is.null(dataValues)) {
            print("No data - aborting.")
            return(NULL)
        }
        ## zero or nearly zero expression leads to crazy etap values, so limit values from below
        dataValues[dataValues<0.1] = 0.1
        ## data label for legend
        if (is.na(dataLabel)) dataLabel = paste(condition,":",gene,sep="")
    }
    
    ## estimate rhop0 from minimum t data points
    tMin = min(dataTimes)
    if (!hasArg(rhop0)) {
        rhop0 = mean(dataValues[dataTimes==tMin])
    }

    ## estimate etap from rhop0, get sign from sign of last-first points
    if (!hasArg(etap)) {
        tMax = max(dataTimes)
        rhope = mean(dataValues[dataTimes==tMax])
        if (rhope>rhop0) {
            etap = 0.1*rhop0
        } else {
            etap = -0.1*rhop0
        }
    }

    if (fitTerms=="nu") {

        ## fix rhop0, etap, gammap; adjust nu
        fit = try(nlm(p=c(nu), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes,dataValues=dataValues,
                      rhoc0=rhoc0, rhop0=rhop0,etap=etap,gammap=gammap ))
        if (class(fit)=="try-error") return(NULL)
        nu = fit$estimate[1]

    } else if (fitTerms=="rhop0") {

        ## fix nu, etap and gammap; adjust rhop0
        fit = try(nlm(p=c(rhop0), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes,dataValues=dataValues,
                      rhoc0=rhoc0,nu=nu, etap=etap,gammap=gammap ))
        if (class(fit)=="try-error") return(NULL)
        rhop0 = fit$estimate[1]

    } else if (fitTerms=="etap") {

        ## fix nu, rhop0 and gammap; adjust etap
        fit = try(nlm(p=c(etap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes,dataValues=dataValues,
                      rhoc0=rhoc0,nu=nu, rhop0=rhop0,gammap=gammap ))
        if (class(fit)=="try-error") return(NULL)
        etap = fit$estimate[1]

    } else if (fitTerms=="gammap") {

        ## fix nu, rhop0, etap; adjust gammap
        fit = try(nlm(p=c(gammap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes,dataValues=dataValues,
                      rhoc0=rhoc0,nu=nu, rhop0=rhop0,etap=etap ))
        if (class(fit)=="try-error") return(NULL)
        gammap = fit$estimate[1]

    } else if (fitTerms=="nu.etap") {

        ## fix rhop0, gammap; adjust nu, etap
        fit = try(nlm(p=c(nu,etap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes,dataValues=dataValues,
                      rhoc0=rhoc0, rhop0=rhop0,gammap=gammap ))
        if (class(fit)=="try-error") return(NULL)
        nu = fit$estimate[1]
        etap = fit$estimate[2]
        
    } else if (fitTerms=="etap.gammap") {
        
        ## fix nu, rhop0; adjust etap, gammap
        fit = try(nlm(p=c(etap,gammap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                      rhoc0=rhoc0,nu=nu, rhop0=rhop0 ))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                          rhoc0=rhoc0, nu=nu, rhop0=rhop0 ))
            if (class(fit)=="try-error") return(NULL)
        }
        etap = fit$estimate[1]
        gammap = fit$estimate[2]

    } else if (fitTerms=="rhop0.gammap") {
        
        ## fix nu, etap; adjust rhop0, gammap
        fit = try(nlm(p=c(rhop0,gammap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                      rhoc0=rhoc0,nu=nu, etap=etap ))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                          rhoc0=rhoc0,nu=nu, etap=etap ))
            if (class(fit)=="try-error") return(NULL)
        }
        rhop0 = fit$estimate[1]
        gammap = fit$estimate[2]

    } else if (fitTerms=="rhop0.etap") {
        
        ## fix nu, gammap; adjust rhop0, etap -- default to refine rhop0 and etap for larger parameter estimates below as well
        fit = try(nlm(p=c(rhop0,etap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms="rhop0.etap", turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                      rhoc0=rhoc0, nu=nu, gammap=gammap ))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                          rhoc0=rhoc0, nu=nu, gammap=gammap ))
            if (class(fit)=="try-error") return(NULL)
        }
        rhop0 = fit$estimate[1]
        etap = fit$estimate[2]
        
    }

    if (fitTerms=="rhop0.etap.gammap") {

        ## fix nu; adjust rhop0, etap, gammap
        fit = try(nlm(p=c(rhop0,etap,gammap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                      rhoc0=rhoc0, nu=nu ))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                          rhoc0=rhoc0, nu=nu ))
            if (class(fit)=="try-error") return(NULL)
        }
        rhop0 = fit$estimate[1]
        etap = fit$estimate[2]
        gammap = fit$estimate[3]

    } else if (fitTerms=="nu.etap.gammap") {

        ## fix rhop0; adjust nu, etap, gammap
        fit = try(nlm(p=c(nu,etap,gammap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                      rhoc0=rhoc0, rhop0=rhop0 ))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                          rhoc0=rhoc0, rhop0=rhop0 ))
            if (class(fit)=="try-error") return(NULL)
        }
        nu = fit$estimate[1]
        etap = fit$estimate[2]
        gammap = fit$estimate[3]

    } else if (fitTerms=="nu.rhop0.etap.gammap") {

        ## adjust nu, rhop0, etap, gammap
        fit = try(nlm(p=c(nu,rhop0,etap,gammap), f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                      rhoc0=rhoc0 ))
        if (class(fit)=="try-error") return(NULL)
        if (fit$code==4 || fit$code==5) {
            ## try again with failed fit params
            fit = try(nlm(p=fit$estimate, f=transmodel.error, steptol=1e-6, gradtol=1e-6, iterlim=1000, fitTerms=fitTerms, turnOff=turnOff, dataTimes=dataTimes, dataValues=dataValues,
                          rhoc0=rhoc0 ))
            if (class(fit)=="try-error") return(NULL)
        }
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

