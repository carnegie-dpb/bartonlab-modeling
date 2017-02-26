source("transmodel.fit.R")
source("Rsquared.R")
source("rhop.R")

source("~/R/getExpression.R")
source("~/R/getAllIDs.R")
source("~/R/getName.R")

##
## Scan ALL genes to assemble a dataframe of fit parameters, choosing best fit from turnOff options.
## Genes are stored if:
##   R2>minR2
##   meanExpr>minExpr
##   etap>minEtap
##   etaphat>minEtap

scanall = function(schema, condition, minR2=0.60, minExpr=10, minEtap=1, minEtapHat=0.1, host="localhost") {

    ## defaults
    nu = 10
    rhon0 = 1
    rhoc0 = 25

    ## get an array of all gene ids
    ids = getAllIDs(schema=schema, host=host)

    fits = data.frame(id=character(),group=character(),rhop0=numeric(),etap=numeric(),gammap=numeric(),minimum=numeric(),R2=numeric(),etap.hat=numeric(),kappa=numeric(),check.rows=T)
    dataTimes = getTimes(schema=schema, condition=condition, host=host)/60
    
    ## loop over all, doing fit
    for (i in 1:length(ids)) {
        dataValues = getExpression(schema=schema, condition=condition, host=host, gene=ids[i])
        meanExpr = mean(dataValues)
        if (meanExpr>=minExpr) {
            fit.E = transmodel.fit(turnOff=0.5, doPlot=F, schema=schema, condition=condition, dataTimes=dataTimes, dataValues=dataValues, host=host)
            fit.L = transmodel.fit(turnOff=0.0, doPlot=F, schema=schema, condition=condition, dataTimes=dataTimes, dataValues=dataValues, host=host)
            if (!is.null(fit.E) && !is.null(fit.L)) {
                ## los fit values
                fitValues.E = rhop(turnOff=0.5, t=dataTimes, rhon0=rhon0, rhoc0=rhoc0, nu=nu, rhop0=fit.E$estimate[1], etap=fit.E$estimate[2], gammap=fit.E$estimate[3])
                fitValues.L = rhop(turnOff=0,   t=dataTimes, rhon0=rhon0, rhoc0=rhoc0, nu=nu, rhop0=fit.L$estimate[1], etap=fit.L$estimate[2], gammap=fit.L$estimate[3])
                ## los r-squared
                R2.L = Rsquared(fitValues.L,dataValues)
                R2.E = Rsquared(fitValues.E,dataValues)
                ## default to L if it's good enough
                if (!is.nan(R2.L) && R2.L>minR2) {
                    rhop0 = fit.L$estimate[1]
                    etap = fit.L$estimate[2]
                    gammap = fit.L$estimate[3]
                    kappa = nu*etap*rhoc0/rhop0
                    etap.hat = etap*rhon0/rhop0
                    if (abs(etap)>minEtap && abs(etap.hat)>minEtapHat) {
                        print(paste(ids[i], getName(ids[i]), "L", round(R2.L,2), round(meanExpr,2), round(etap,2), round(etap.hat,2), round(kappa,2)), quote=F)
                        fits = rbind(fits, data.frame(id=ids[i],group="L",rhop0=rhop0,etap=etap,gammap=gammap,minimum=fit.L$minimum,R2=R2.L,etap.hat=etap.hat,kappa=kappa))
                    }
                } else if (!is.nan(R2.E) && R2.E>minR2) {
                    rhop0 = fit.E$estimate[1]
                    etap = fit.E$estimate[2]
                    gammap = fit.E$estimate[3]
                    kappa = nu*etap*rhoc0/rhop0
                    etap.hat = etap*rhon0/rhop0
                    if (abs(etap)>minEtap && abs(etap.hat)>minEtapHat) {
                        print(paste(ids[i], getName(ids[i]), "E", round(R2.E,2), round(meanExpr,2), round(etap,2), round(etap.hat,2), round(kappa,2)), quote=F)
                        fits = rbind(fits, data.frame(id=ids[i],group="E", rhop0=rhop0,etap=etap,gammap=gammap,minimum=fit.E$minimum, R2=R2.E,etap.hat=etap.hat,kappa=kappa))
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
