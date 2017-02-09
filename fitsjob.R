##
## load the DE data write to files
##

source("~/R/getAllIDs.R")
source("compare.fits.R")

micro.genes = getAllIDs(schema="gse30703")
micro.genes = micro.genes[substring(micro.genes,3,3)=="1"|substring(micro.genes,3,3)=="2"|substring(micro.genes,3,3)=="3"|substring(micro.genes,3,3)=="4"|substring(micro.genes,3,3)=="5"]

rnaseq.genes = getAllIDs(schema="gse70796")
rnaseq.genes = rnaseq.genes[substring(rnaseq.genes,3,3)=="1"|substring(rnaseq.genes,3,3)=="2"|substring(rnaseq.genes,3,3)=="3"|substring(rnaseq.genes,3,3)=="4"|substring(rnaseq.genes,3,3)=="5"]

micro.rev.fits = compare.fits(schema="gse30703", condition="GR-REV", genes=micro.genes)
rownames(micro.rev.fits) = micro.rev.fits$gene
micro.rev.fits$gene = NULL
write.table(file="micro.rev.fits.txt", micro.rev.fits)

rnaseq.rev.fits = compare.fits(schema="gse70796", condition="GR-REV", genes=rnaseq.genes)
rownames(rnaseq.rev.fits) = rnaseq.rev.fits$gene
rnaseq.rev.fits$gene = NULL
write.table(file="rnaseq.rev.fits.txt", rnaseq.rev.fits)

