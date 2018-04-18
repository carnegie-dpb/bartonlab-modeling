source("~/R/getAllIDs.R")
source("compare.fits.R")

##
## load the DE data write to files
##

## ## only chromosomal genes
## micro.genes = getAllIDs(schema="gse30703")
## micro.genes = micro.genes[substring(micro.genes,3,3)=="1"|substring(micro.genes,3,3)=="2"|substring(micro.genes,3,3)=="3"|substring(micro.genes,3,3)=="4"|substring(micro.genes,3,3)=="5"]

## only chromosomal genes
rnaseq.genes = getAllIDs(schema="gse70796")
rnaseq.genes = rnaseq.genes[substring(rnaseq.genes,3,3)=="1"|substring(rnaseq.genes,3,3)=="2"|substring(rnaseq.genes,3,3)=="3"|substring(rnaseq.genes,3,3)=="4"|substring(rnaseq.genes,3,3)=="5"]

## ## get microarray data
## micro.kan.fits = compare.fits(schema="gse30703", condition="GR-KAN", fitTerms="etap.gammap", genes=micro.genes)
## rownames(micro.kan.fits) = micro.kan.fits$gene
## micro.kan.fits$gene = NULL
## write.table(file="micro.kan.fits.txt", micro.kan.fits)

## get RNA-seq data
rnaseq.kan.fits = compare.fits(schema="gse70796", condition="GR-KAN", fitTerms="etap.gammap", genes=rnaseq.genes)
rownames(rnaseq.kan.fits) = rnaseq.kan.fits$gene
rnaseq.kan.fits$gene = NULL
write.table(file="rnaseq.kan.fits.txt", rnaseq.kan.fits)

