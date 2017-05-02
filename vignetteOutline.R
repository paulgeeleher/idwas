### R code from vignette source 'vignetteOutline.Snw'
### Encoding: UTF-8

###################################################
### code chunk number 1: vignetteOutline.Snw:11-15
###################################################
library(idwas)
library(ggplot2)
library(pRRophetic)
set.seed(12345)


###################################################
### code chunk number 2: vignetteOutline.Snw:20-24
###################################################
data("testIdwas") # clinicalExpressionData_tcga, clinicalGenomicsData_tcga, clinicalCoVariatesMatrix_tcga: This is a subset of the TCGA expression and genomics data.

dataOut <- idwas(clinicalExpressionData_tcga, clinicalGenomicsData_tcga, clinicalCoVariatesMatrix_tcga, predictableDrugs=c("Nutlin.3a", "PD.0332991", "PD.0325901"))
print(dataOut[1:10,])


