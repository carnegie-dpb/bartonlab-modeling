source("~/modeling/rhon.R")
source("~/modeling/rhop.R")
source("~/modeling/rhos.R")

##
## return the error metric for modeling an indirect target, given parameters for the TF and primary target
##

transmodel2.error = function(p, rhoc0, rhon0, nu, etap, gammap, dataTimes, data2Values, gammasMax=8) {

    ## secondary target fit parameters
    rhos0 = p[1]
    etas = p[2]
    gammas = p[3]

    ## set fit=0 if parameters out of bounds
    if (gammas>gammasMax) {
        fitValues = dataTimes*0
    } else {
        fitValues = rhos(t=dataTimes, rhoc0=rhoc0, rhon0=rhon0, nu=nu, etap=etap, gammap=gammap, rhos0=rhos0, etas=etas, gammas=gammas)
    }
    
    return(errorMetric(fitValues,data2Values))

}
