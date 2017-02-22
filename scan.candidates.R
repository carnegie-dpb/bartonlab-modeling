source("transmodel.fit.R")
source("Rsquared.R")
source("rhop.R")

source("~/R/getTimes.R")
source("~/R/getExpression.R")
source("~/R/getSymbol.R")

##
## Scan over candidate genes to assemble a dataframe of fit parameters
## Input an array of candidate gene IDs
## Returns a dataframe of fits.
##
scan.candidates = function(schema, condition, ids, rhon0=1, rhoc0=19, nu=10, R2min=0.60) {

    ## loop over candidates, doing fit to three turnOff choices, saving if one exceeds R2min
    fits = data.frame(id=character(), symbol=character(), 
                      rhop0=numeric(),     etap=numeric(),     gammap=numeric(),     R2=numeric(), 
                      rhop0.30=numeric(), etap.30=numeric(), gammap.30=numeric(), R2.30=numeric(), 
                      ## rhop0.60=numeric(),etap.60=numeric(),gammap.60=numeric(),R2.60=numeric(),
                      check.rows=T
                      )
    dataTimes = getTimes(schema=schema, condition=condition)/60
    for (i in 1:length(ids)) {
        ## expression values
        dataValues = getExpression(schema=schema, condition=condition, gene=ids[i])
        ## fits
        fit = transmodel.fit(turnOff=0.0, doPlot=F,schema=schema,condition=condition, rhon0=rhon0,rhoc0=rhoc0,nu=nu, dataTimes=dataTimes, dataValues=dataValues)
        fit.30 = transmodel.fit(turnOff=30, doPlot=F,schema=schema,condition=condition, rhon0=rhon0,rhoc0=rhoc0,nu=nu, dataTimes=dataTimes, dataValues=dataValues)
        ## fit.60 = transmodel.fit(turnOff=60, doPlot=F,schema=schema,condition=condition, rhon0=rhon0,rhoc0=rhoc0,nu=nu, dataTimes=dataTimes, dataValues=dataValues)
        ## save if one fit is significant
        R2 = 0
        R2.30 = 0
        ## R2.60 = 0
        save = FALSE
        if (fit$code!=4 && fit$code!=5) {
            fitValues = rhop(t=dataTimes, turnOff=0.0, rhoc0=rhoc0, nu=nu, rhop0=fit$estimate[1], etap=fit$estimate[2], gammap=fit$estimate[3])
            R2 = Rsquared(fitValues,dataValues)
            gammap = fit$estimate[3]
            if (R2>R2min && gammap>0) save = TRUE
        }
        if (fit.30$code!=4 && fit.30$code!=5) {
            fitValues = rhop(t=dataTimes, turnOff=30, rhoc0=rhoc0, nu=nu, rhop0=fit.30$estimate[1], etap=fit.30$estimate[2], gammap=fit.30$estimate[3])
            R2.30 = Rsquared(fitValues,dataValues)
            gammap.30 = fit.30$estimate[3]
            if (R2.30>R2min && gammap.30>0) save = TRUE
        }
        ## if (fit.60$code!=4 && fit.60$code!=5) {
        ##   fitValues = rhop(t=dataTimes, turnOff=60, rhoc0=rhoc0, nu=nu, rhop0=fit.60$estimate[1], etap=fit.60$estimate[2], gammap=fit.60$estimate[3])
        ##   R2.60 = Rsquared(fitValues,dataValues)
        ##   gammap.60 = fit.60$estimate[3]
        ##   if (R2.60>R2min && gammap.60>0) save = TRUE
        ## }
        if (save) {
            symbol = getSymbol(ids[i])
            if (is.null(symbol)) symbol = ""
            ## print(paste(ids[i],symbol,R2,R2.30,R2.60))
            print(paste(ids[i],symbol,R2,R2.30))
            fits = rbind(fits, data.frame(id=ids[i],symbol=symbol,
                                          rhop0=fit$estimate[1], etap=fit$estimate[2], gammap=fit$estimate[3], R2=R2,
                                          rhop0.30=fit.30$estimate[1], etap.30=fit.30$estimate[2], gammap.30=fit.30$estimate[3], R2.30=R2.30
                                          ## rhop0.60=fit.60$estimate[1], etap.60=fit.60$estimate[2], gammap.60=fit.60$estimate[3], R2.60=R2.60
                                          ))
        }
    }

    ## move id over to row name
    rownames(fits) = fits$id
    fits$id = NULL
    
    return(fits)

}
