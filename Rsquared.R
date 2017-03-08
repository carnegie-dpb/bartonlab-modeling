##
## return R-squared for a set of fit and data values
##

Rsquared = function(fitValues, dataValues) {
  dataMean = mean(dataValues)
  numer = sum( (dataValues-fitValues)^2 )
  denom = sum( (dataValues-dataMean)^2 )
  r2 = 1.0 - numer/denom
  return(r2)
}
