##
## return the group membership (E, M, L) of the given schema, condition, gene
##

source("~/R/getExpression.R")
source("~/R/getSymbol.R")
source("~/R/getTimes.R")

getGroup = function(schema="gse70796", condition="GR-REV", gene="HAT22") {

    t = getTimes(schema, condition)
    expr = getExpression(schema, condition, gene)

    if (length(expr)==0) {
        return()
    }

    t.0 = (t==0)
    t.early = (t==30)
    t.middle = (t==60)
    t.late = (t==120)

    expr.0 = mean(expr[t.0])
    
    expr.early = abs(log2(mean(expr[t.early])/expr.0))
    expr.middle = abs(log2(mean(expr[t.middle])/expr.0))
    expr.late = abs(log2(mean(expr[t.late])/expr.0))

    if (expr.early>expr.middle && expr.early>expr.late) {
        return("E")
    } else if (expr.middle>expr.early && expr.middle>expr.late) {
        return("M")
    } else {
        return("L")
    }

}

    
