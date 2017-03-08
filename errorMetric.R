##
## returns the error metric used for transmodel fits
##

errorMetric = function(fitValues, dataValues) {

    ## ## sum of squared deviation scaled by sum(data^2)
    ## numer = 0
    ## denom = 0
    ## for (i in 1:length(dataValues)) {
    ##     numer = numer + (dataValues[i]-fitValues[i])^2
    ##     denom = denom + dataValues[i]^2
    ## }
    ## err = numer/denom

    ## 1 - r^2 : seems to produce sharper gradients, "better" fits (and maximizes r^2)
    dataMean = mean(dataValues)
    numer = sum( (dataValues-fitValues)^2 )
    denom = sum( (dataValues-dataMean)^2 )
    err = numer/denom

    ## ## sum of squared deviation - seems to converge the fastest
    ## err = sum( (dataValues-fitValues)^2 )

    return(err)
    
}
