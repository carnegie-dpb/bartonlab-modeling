##
## plot REV vs t for the RNA-seq data to show that REV level is constant in GR-REV samples

source("~/R/getExpression.R")
source("~/R/getTimes.R")

t = getTimes(schema="gse70796", condition="ALL")
REV = getExpression(schema="gse70796", condition="ALL", gene="REV")

plot.bars(t[REV<10]/60,REV[REV<10], log="y", ylim=c(1,200), pch=1, cex=1.5, col="black",
     xlab="time after DEX application (h)", ylab="REV mRNA level (FPKM)"
     )
plot.bars(t[REV>10]/60,REV[REV>10], pch=19, cex=1.5, over=TRUE)

legend(x=0, y=25, pch=c(19,1), pt.cex=1.5, y.intersp=1.5, legend=c("GR-REV", "WT, GR-AS2, GR-KAN1, GR-STM"))

