source("~/modeling/rhop.R")
source("~/modeling/errorMetric.R")

## return the error metric for modeling a direct target, for nlm usage
##
## assumes no nuclear loss of TF (gamman=0)
## gammapMax is upper bound on allowable gammap value; if above that, fit is set to zero

transmodel.error = function(p, fitTerms, turnOff, rhoc0,rhon0,nu, rhop0,etap,gammap, dataTimes, dataValues, gammapMax=5) {

    ## have do do this before gammap is set below
    gammapSupplied = hasArg(gammap)

    if (fitTerms=="nu") {
        nu = p[1]
    } else if (fitTerms=="rhop0") {
        rhop0 = p[1]
    } else if (fitTerms=="etap") {
        etap = p[1]
    } else if (fitTerms=="gammap") {
        gammap = p[1]
    } else if (fitTerms=="nu.etap") {
        nu = p[1]
        etap = p[2]
    } else if (fitTerms=="etap.gammap") {
        etap = p[1]
        gammap = p[2]
    } else if (fitTerms=="rhop0.etap") {
        rhop0 = p[1]
        etap = p[2]
    } else if (fitTerms=="rhop0.gammap") {
        rhop0 = p[1]
        gammap = p[2]
    } else if (fitTerms=="rhop0.etap.gammap") {
        rhop0 = p[1]
        etap = p[2]
        gammap = p[3]
    } else if (fitTerms=="nu.etap.gammap") {
        nu = p[1]
        etap = p[2]
        gammap = p[3]
    } else if (fitTerms=="nu.rhop0.etap.gammap") {
        nu = p[1]
        rhop0 = p[2]
        etap = p[3]
        gammap = p[4]
    }

    if (!gammapSupplied && gammap>gammapMax) {
        ## upper limit on gamma required
        fitValues = dataTimes*0
    } else {
        fitValues = rhop(t=dataTimes, rhoc0=rhoc0, rhon0=rhon0, nu=nu, rhop0=rhop0, etap=etap, gammap=gammap, turnOff=turnOff)
    }

    ## return error metric
    return(errorMetric(fitValues,dataValues))

}
