##
## returns the error metric used for transmodel fits
##

errorMetric = function(fitValues, dataValues) {

    ## ## sum of square deviation normalized by sum(data^2)
    ## numer = 0
    ## denom = 0
    ## for (i in 1:length(dataValues)) {
    ##     numer = numer + (dataValues[i]-fitValues[i])^2
    ##     denom = denom + dataValues[i]^2
    ## }
    ## err = numer/denom

    ## least squares - seems to converge the fastest
    err = sum( (dataValues-fitValues)^2 )

    return(err)
    
}
