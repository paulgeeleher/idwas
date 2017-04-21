#' Given an input of expression data and genomics data for a set of tumors, run IDWAS.
#' 
#' This function takes as its input gene expression and genomics data from a set of cancer
#' patients' tumors. The purpose of this function is to impute response in these tumors
#' using the pRRophetic prediction library, then to identify variables in the genomics data
#' supplied that are correlated with predicted drug response.
#'
#' @param clinicalExpressionData Numeric matrix of tumor gene expresison data. Colum names must be patient IDs, row names must be HUGO gene symbols.
#' @param clinicalGenomicsData Numeric matrix representing another source of data that will be comapred to imputed drug response. E.g. genomics data, whereby a mutated gene as been coded as a 1 and wt as 0. Rownames must identify the abberation (e.g. gene symbols) and column names must represent patients. Column names must be the same as the gene expression data.
#' @param clinicalCoVariatesMatrix A matrix of co-variates to included in the association analysis, e.g. cancer type or GLDS.
#' @param predictionLibrary What should be used to impute drug response. Currently, pRRophetic is supported, but in theory any prediction method could be implemented.
#' @param predictionFunction What function within the predictionLibrary should be used, currently pRRopheticPredict is supported.
#' @param predictableDrugs A vector of drugs for which we want to predict drug response. By default it is a list of drug with Spearman correlation of greater than 0.3 in cross-validation on the GDSC 2014 dataset.
#'
#' @return a list of matrices of p-values and effect sizes, for each drug against each aberration in clinicalGenomicsData.
#'
#' @import pRRophetic
#'
#' @keywords perform idwas
#'
#' @export
idwas <- function(clinicalExpressionData, clinicalGenomicsData, clinicalCoVariatesMatrix=NULL, predictionLibrary="pRRophetic", predictionFunction="pRRopheticPredict", predictableDrugs=c("ABT.263", "ABT.888", "AG.014699", "AICAR", "ATRA", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "BAY.61.3606", "BIBW2992", "Bicalutamide", "BI.D1870", "Bleomycin", "BMS.536924", "BMS.754807", "Bortezomib", "Bosutinib", "BX.795", "Camptothecin", "CEP.701", "CGP.082996", "CHIR.99021", "CI.1040", "Cisplatin", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Elesclomol", "Erlotinib", "FH535", "FTI.277", "Gefitinib", "Gemcitabine", "IPA.3", "Lapatinib", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "Nilotinib", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "PAC.1", "PD.0325901", "PD.0332991", "PD.173074", "PLX4720", "RDEA119", "SB590885", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vorinostat", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "ZM.447439"), ...)
{
  # check if both row and column names have been specified
  if(is.null(rownames(clinicalExpressionData)) || is.null(rownames(clinicalExpressionData)) || is.null(rownames(clinicalGenomicsData)) || is.null(rownames(clinicalGenomicsData)))
  {
    stop("ERROR: Gene identifiers must be specified as \"rownames()\" on both matrices. ")
  }
  
  # Are the expression data, the genomics data and the coVariates data (if supplied) of class "matrix"
  if((class(clinicalExpressionData) != "matrix") || (class(clinicalGenomicsData) != "matrix"))
  {
    stop("ERROR: clinicalExpressionData and clinicalGenomicsData must be of classl \"matrix\".")
  }
  
  # If clinical covariates were supplied, are they also of class "matrix"
  if(!is.null(clinicalCoVariatesMatrix))
  {
    if(class(clinicalCoVariatesMatrix) != "matrix")
      stop("ERROR: clinicalCoVariatesMatrix must be of classl \"matrix\". See function model.matrix() for information on converting co-variates (e.g. categorical variables) to a design matrix format. Also discussed in package Vignette.")
  }
  
  # subset and order the datasets provided
  if(!is.null(clinicalCoVariatesMatrix))
  {
    samplesCommon <- intersect(intersect(colnames(clinicalExpressionData), colnames(clinicalGenomicsData)), colnames(clinicalCoVariatesMatrix))
    clinicalExpressionData <- clinicalExpressionData[, samplesCommon]
    clinicalGenomicsData <- clinicalGenomicsData[, samplesCommon]
    clinicalCoVariatesMatrix <- clinicalCoVariatesMatrix[, samplesCommon]
  }
  else
  {
    samplesCommon <- intersect(colnames(clinicalExpressionData), colnames(clinicalGenomicsData))
    clinicalExpressionData <- clinicalExpressionData[, samplesCommon]
    clinicalGenomicsData <- clinicalGenomicsData[, samplesCommon]
  }
 
  # Stop if there are no samples common to all of these datasets.
  if(length(samplesCommon) < 1)
  {
    stop("ERROR: No samples found common to all data matrices provided. Please ensure that the column names of each matrix is set to a unique patient id and that the same ids are used across expression, genomic and co-variate matrices.")
  }
  
  # first create the predictions.
  if(predictionLibrary == "pRRophetic")
  {
    if(predictionFunction == "pRRopheticPredict")
    {
      # create  prediction for each drug in each sample
      drugPredList <- list()
      cat("Imputing drug response scores for all drugs:\n\n")
      for(i in 1:length(predictableDrugs))
      {
	drugPredList[[i]] <-  pRRopheticPredict(clinicalExpressionData, predictableDrugs[i], selection=1, batchCorrect="standardize", removeLowVaryingGenes=0.2, removeLowVaringGenesFrom="rawData")
      }
      drugPredMat <- do.call(rbind, drugPredList)
      rownames(drugPredMat) <- predictableDrugs
    }
  }
  else
  {
    stop("Unrecognized prediction library")
  }
  
  # Now do the IDWAS, between the predicted responses and the imputed drug data.
  pValList <- list()
  betaValList <- list()
  for(i in 1:nrow(drugPredMat))
  {
    cat(paste("Calculating P-values for ", predictableDrugs[i], "...\n", sep=""))
    pValList[[i]] <- numeric()
    betaValList[[i]] <- numeric()
    for(j in 1:nrow(clinicalGenomicsData))
    {
      thecoefs <- coef(summary(lm(drugPredMat[i,]~clinicalGenomicsData[j,]+t(clinicalCoVariatesMatrix))))
      pValList[[i]][[j]] <- thecoefs[2,4]
      betaValList[[i]][[j]] <- thecoefs[2,1]
    }
  }
  
  # Create FDDR, adjusted p-values for the number of models, output all as a matrix.
  cat("Calculating adjusted P-values and outputting data...\n")
  pAdjList <- list()
  for(i in 1:length(pValList))
  {
    names(pValList[[i]]) <- rownames(clinicalGenomicsData)
    names(betaValList[[i]]) <- rownames(clinicalGenomicsData)
    pAdjList[[i]] <- p.adjust(pValList[[i]], method="BH")
  }
  names(pValList) <- rownames(drugPredMat)
  names(betaValList) <- rownames(drugPredMat)
  names(pAdjList) <- rownames(drugPredMat)

  # We only want to report the p-values for the "predictable" drugs....
  pValList_predictable <- unlist(pValList[predictableDrugs])
  ord <- order(pValList_predictable)
  pValList_predictable_ord <- pValList_predictable[ord]
  betaValList_predictable <- unlist(betaValList[predictableDrugs])[ord]
  pAdjListCantype_predictable <- unlist(pAdjList[predictableDrugs])[ord]
  pAdjListCantype_predictable_forNumModels <- (pAdjListCantype_predictable*length(predictableDrugs))
  pAdjListCantype_predictable_forNumModels[pAdjListCantype_predictable_forNumModels > 1] <- 1 # adjust the FDRs for the number of models and cap this at 1. This is equivalent to a Bonferroni correction for the number of drugs/models.
  outTab <- cbind(pValList_predictable_ord, pAdjListCantype_predictable, betaValList_predictable, pAdjListCantype_predictable_forNumModels)
  colnames(outTab) <- c("P-value", "FDR", "Beta", "FDR Adjusted for Number of Drug-Models")
  return(outTab)
}
