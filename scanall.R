source("~/modeling/transmodel.fit.R")
source("~/modeling/Rsquared.R")
source("~/modeling/rhop.R")

source("~/R/getExpression.R")
source("~/R/getAllIDs.R")
source("~/R/getName.R")

## Scan ALL genes to assemble a dataframe of fit parameters, choosing early turnOff if it fits and no turnOff doesn't.
## Genes are stored if:
##   R2>minR2
##   meanExpr>minExpr

scanall = function(schema, condition, minR2=0.60, minExpr=1, nu=10, host="localhost") {

    ## fixed parameters
    rhon0 = 1
    rhoc0 = 25

    ## our output data frame
    fits = data.frame(id=character(),name=character(),group=character(),
                      rhop0=numeric(),etap=numeric(),etap.hat=numeric(),gammap=numeric(),
                      kappa=numeric(),logFCinf=numeric(),R2=numeric(),check.rows=T)

    ## get an array of all gene ids
    ids = getAllIDs(schema=schema, host=host)
    dataTimes = getTimes(schema=schema, condition=condition, host=host)/60
    
    ## loop over all, doing fit
    for (i in 1:length(ids)) {
        dataValues = getExpression(schema=schema, condition=condition, host=host, gene=ids[i])
        meanExpr = mean(dataValues)
        name = getName(ids[i])
        if (meanExpr>=minExpr) {
            fit.L.good = FALSE
            fit.L = transmodel.fit(turnOff=0.0, doPlot=F, schema=schema, condition=condition, dataTimes=dataTimes, dataValues=dataValues, host=host, nu=nu)
            if (!is.null(fit.L) && !is.nan(fit.L$estimate[1]) && !is.nan(fit.L$estimate[2]) && !is.nan(fit.L$estimate[3])) {
                fitValues.L = rhop(turnOff=0.0, t=dataTimes, rhon0=rhon0, rhoc0=rhoc0, nu=nu, rhop0=fit.L$estimate[1], etap=fit.L$estimate[2], gammap=fit.L$estimate[3])
                R2.L = Rsquared(fitValues.L,dataValues)
                fit.L.good = !is.nan(R2.L) && R2.L>minR2
            }
            fit.E.good = FALSE
            fit.E = transmodel.fit(turnOff=0.5, doPlot=F, schema=schema, condition=condition, dataTimes=dataTimes, dataValues=dataValues, host=host, nu=nu)
            if (!is.null(fit.E) && !is.nan(fit.E$estimate[1]) && !is.nan(fit.E$estimate[2]) && !is.nan(fit.E$estimate[3])) {
                fitValues.E = rhop(turnOff=0.5, t=dataTimes, rhon0=rhon0, rhoc0=rhoc0, nu=nu, rhop0=fit.E$estimate[1], etap=fit.E$estimate[2], gammap=fit.E$estimate[3])
                R2.E = Rsquared(fitValues.E,dataValues)
                fit.E.good = !is.nan(R2.E) && R2.E>minR2
            }
            if (fit.L.good && (!fit.E.good || R2.L>R2.E-0.1)) {
                rhop0 = fit.L$estimate[1]
                etap = fit.L$estimate[2]
                gammap = fit.L$estimate[3]
                logFCinf = log2(1 + etap/gammap*rhoc0/rhop0)
                kappa = nu*etap*rhoc0/rhop0
                etap.hat = etap*rhon0/rhop0
                R2 = R2.L
                print(paste(ids[i], name, "L", round(meanExpr,2), round(etap,2), round(etap.hat,2), round(gammap,2), round(kappa,2), round(logFCinf,2), round(R2,2)), quote=F)
                fits = rbind(fits, data.frame(id=ids[i],name=name,group="L", rhop0=rhop0,etap=etap,etap.hat=etap.hat,gammap=gammap, kappa=kappa,logFCinf=logFCinf,R2=R2))
            } else if (fit.E.good) {
                rhop0 = fit.E$estimate[1]
                etap = fit.E$estimate[2]
                gammap = fit.E$estimate[3]
                logFCinf = log2(1 + etap/gammap*rhoc0/rhop0)
                kappa = nu*etap*rhoc0/rhop0
                etap.hat = etap*rhon0/rhop0
                R2 = R2.E
                print(paste(ids[i], name, "E", round(meanExpr,2), round(etap,2), round(etap.hat,2), round(gammap,2), round(kappa,2), round(logFCinf,2), round(R2,2)), quote=F)
                fits = rbind(fits, data.frame(id=ids[i],name=name,group="E", rhop0=rhop0,etap=etap,etap.hat=etap.hat,gammap=gammap, kappa=kappa,logFCinf=logFCinf,R2=R2))
            }
        }
    }

    ## move id to rowname
    rownames(fits)=fits$id
    fits$id = NULL

    return(fits)
}
