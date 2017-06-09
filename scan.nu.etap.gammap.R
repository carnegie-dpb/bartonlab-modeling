source("transmodel.fit.R")
source("Rsquared.R")
source("rhop.R")

source("~/R/getExpression.R")
source("~/R/getAllIDs.R")
source("~/R/getName.R")

## Scan ALL genes doing nu.etap.gammap fits to assemble a dataframe of fit parameters, choosing best fit from turnOff options.
## Genes are fit if:
##   meanExpr>minExpr
## Genes are stored if:
##   R2>minR2 (L fit trumps E fit if both acceptable)

scan.nu.etap.gammap = function(schema, condition, minR2=0.60, minExpr=1.0, host="localhost") {

    ## fixed values
    rhon0 = 1
    rhoc0 = 25

    ## get an array of all gene ids
    ids = getAllIDs(schema=schema, host=host)

    fits = data.frame(id=character(),name=character(),group=character(),
                      nu=numeric(), rhop0=numeric(), etap=numeric(), etap.hat=numeric(), gammap=numeric(),
                      kappa=numeric(), logFCinf=numeric(), R2=numeric(), check.rows=T)
    dataTimes = getTimes(schema=schema, condition=condition, host=host)/60
    tMin = min(dataTimes)
    
    ## loop over all ids, doing L and E fits
    for (i in 1:length(ids)) {
        dataValues = getExpression(schema=schema, condition=condition, host=host, gene=ids[i])
        meanExpr = mean(dataValues)
        name = getName(ids[i])
        if (meanExpr>=minExpr) {
            rhop0 = mean(dataValues[dataTimes==tMin])
            fit.L = transmodel.fit(turnOff=0.0, doPlot=F, schema=schema, condition=condition, dataTimes=dataTimes, dataValues=dataValues, host=host, rhop0=rhop0, fitTerms="nu.etap.gammap")
            if (!is.null(fit.L)) {
                fitValues.L = rhop(turnOff=0.0, t=dataTimes, rhon0=rhon0, rhoc0=rhoc0, rhop0=rhop0, nu=fit.L$estimate[1], etap=fit.L$estimate[2], gammap=fit.L$estimate[3])
                R2.L = Rsquared(fitValues.L,dataValues)
                fit.L.good = !is.nan(R2.L) && R2.L>minR2
                if (fit.L.good) {
                    nu = fit.L$estimate[1]
                    etap = fit.L$estimate[2]
                    gammap = fit.L$estimate[3]
                    logFCinf = log2(1 + etap/gammap*rhoc0/rhop0)
                    kappa = nu*etap*rhoc0/rhop0
                    etap.hat = etap*rhon0/rhop0
                    print(paste(ids[i], name, "L", round(nu,2), round(rhop0,2), round(etap,2), round(etap.hat,2), round(gammap,2), round(kappa,2), round(logFCinf,2), round(R2.L,2)), quote=F)
                    fits = rbind(fits, data.frame(id=ids[i],name=name,group="L",nu=nu,rhop0=rhop0,etap=etap,etap.hat=etap.hat,gammap=gammap,kappa=kappa,logFCinf=logFCinf,R2=R2.L))
                }
            }
            if (!fit.L.good) {
                fit.E = transmodel.fit(turnOff=0.5, doPlot=F, schema=schema, condition=condition, dataTimes=dataTimes, dataValues=dataValues, host=host, rhop0=rhop0, fitTerms="nu.etap.gammap")
                if (!is.null(fit.E)) {
                    fitValues.E = rhop(turnOff=0.5, t=dataTimes, rhon0=rhon0, rhoc0=rhoc0, rhop0=rhop0, nu=fit.E$estimate[1], etap=fit.E$estimate[2], gammap=fit.E$estimate[3])
                    R2.E = Rsquared(fitValues.E,dataValues)
                    fit.E.good = !is.nan(R2.E) && R2.E>minR2
                    if (fit.E.good) {
                        nu = fit.E$estimate[1]
                        etap = fit.E$estimate[2]
                        gammap = fit.E$estimate[3]
                        logFCinf = log2(1 + etap/gammap*rhoc0/rhop0)
                        kappa = nu*etap*rhoc0/rhop0
                        etap.hat = etap*rhon0/rhop0
                        print(paste(ids[i], name, "E", round(nu,2), round(rhop0,2), round(etap,2), round(etap.hat,2), round(gammap,2), round(kappa,2), round(logFCinf,2), round(R2.E,2)), quote=F)
                        fits = rbind(fits, data.frame(id=ids[i],name=name,group="E", nu=nu,rhop0=rhop0,etap=etap,etap.hat=etap.hat,gammap=gammap,kappa=kappa,logFCinf=logFCinf,R2=R2.E))
                    }
                }
            }
            if (!fit.L.good && !fit.E.good) {
                fit.M = transmodel.fit(turnOff=1.0, doPlot=F, schema=schema, condition=condition, dataTimes=dataTimes, dataValues=dataValues, host=host, rhop0=rhop0, fitTerms="nu.etap.gammap")
                if (!is.null(fit.M)) {
                    fitValues.M = rhop(turnOff=1.0, t=dataTimes, rhon0=rhon0, rhoc0=rhoc0, rhop0=rhop0, nu=fit.E$estimate[1], etap=fit.E$estimate[2], gammap=fit.E$estimate[3])
                    R2.M = Rsquared(fitValues.M,dataValues)
                    fit.M.good = !is.nan(R2.M) && R2.M>minR2
                    if (fit.M.good) {
                        nu = fit.M$estimate[1]
                        etap = fit.M$estimate[2]
                        gammap = fit.M$estimate[3]
                        logFCinf = log2(1 + etap/gammap*rhoc0/rhop0)
                        kappa = nu*etap*rhoc0/rhop0
                        etap.hat = etap*rhon0/rhop0
                        print(paste(ids[i], name, "M", round(nu,2), round(rhop0,2), round(etap,2), round(etap.hat,2), round(gammap,2), round(kappa,2), round(logFCinf,2), round(R2.M,2)), quote=F)
                        fits = rbind(fits, data.frame(id=ids[i],name=name,group="E", nu=nu,rhop0=rhop0,etap=etap,etap.hat=etap.hat,gammap=gammap,kappa=kappa,logFCinf=logFCinf,R2=R2.M))
                    }
                }
            }
        }
    }

    ## move id to rowname
    rownames(fits)=fits$id
    fits$id = NULL

    return(fits)
}
