#' CalcSpearmanMatrix() Function
#'
#' Function to calculate spearman correlation of OTU abundance to DNA concentration
#'
#' Returns a list of 3 matrices:
#'  * retMatrix -- full matrix of results (Rho, P.value, FDR)
#'  * negSigMatrix -- subset to stat sig (p < 0.05) and negatively correlated OTUs
#'  * posSigMatrix -- subset to stat sig (p < 0.05) and positively correlated OTUs
#'
#' @param otuMatrix (otu_table(phyloseq_object)) works fine
#' @param intensityVector list of concentrations in the same order (sample_data(phyloseq_object)$Concentration)
#' @keywords spearman, correlation, jervis-bardy, contamination
#' @export
#' @examples
#' JB_corr <- CalcSpearmanMatrix(phyloseq::otu_table(ps$r),phyloseq::sample_data(psc$r)$DNA_Concentration)
CalcSpearmanMatrix = function(otuMatrix, intensityVector) {

  #check whether the matrix needs to be transposed
  nsamples <- length(intensityVector)
  if(nsamples != ncol(otuMatrix)){
    otuMatrix<- t(otuMatrix)
  }

  retMatrix = matrix(ncol=3, nrow=0)
  rowNames = c()


  for(i in 1:nrow(otuMatrix)) {
    otu = rownames(otuMatrix)[i]
    otuVector = as.vector(unlist(otuMatrix[i,]))

    temp = cor.test(intensityVector, otuVector, method='spearman', conf.level=0.95, exact=FALSE)

    rho = as.numeric(temp$estimate)
    pval = temp$p.value

    retMatrix = rbind(retMatrix, c(rho, pval, NA))
    rowNames = append(rowNames, otu)
  }

  #assign row & colnames and correct the P-value
  rownames(retMatrix) = rowNames
  colnames(retMatrix) = c('Rho','P.value','FDR')
  retMatrix[,3] = p.adjust(retMatrix[,2], method='BH')

  #Subset to pos and negative matrices
  neg <-  retMatrix[ ,1] < 0
  negMatrix <- retMatrix[neg, ] #subset to the negatively correlated ones
  pos <-  retMatrix[ ,1] > 0
  posMatrix <- retMatrix[pos, ] #subset to the negatively correlated ones

  #then subset to jsut the statsig ones ( p < 0.05)
  negsig <- negMatrix[ ,2] < 0.05
  negSigMatrix <- negMatrix[negsig, ] #subset to the statsig ones
  possig <- posMatrix[ ,2] < 0.05
  posSigMatrix <- posMatrix[possig, ] #subset to the statsig ones

  retBundle = list()
  retBundle$retMatrix <- retMatrix
  retBundle$negSigMatrix <- negSigMatrix
  retBundle$posSigMatrix <- posSigMatrix

  return(retBundle)
}
