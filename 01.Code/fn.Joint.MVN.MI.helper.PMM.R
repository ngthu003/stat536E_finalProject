#' Joint Multivariate Normal MI Analysis
#'
#' fn.Joint.MVN.MI.helper.PMM is a helper subfunction to impute via 
#' - PMM: predictive mean meatching
#' instead of the default
#' - PPD: posterior predictive distribution
#' in Joint MVN MI
#'
#' @param dgp          # randomly generated dataset w/ missing data
#' @param m            #(imputation)
#' @param imputed.Data imputed data
#' @param col.idx.Y    column index of Y


fn.Joint.MVN.MI.helper.PMM <- function(
    dgp,
    m, 
    imputed.Data,
    col.idx.Y
) {
  
  cmeans <- vector("list", m)
  idx    <- vector("list", m)
  draws  <- vector("list", m)
  
  for(k in seq_along(imputed.Data$imputations)){
    data.cmu <- imputed.Data$imputations[[k]][,-col.idx.Y]
    mu       <- imputed.Data$mu[,k]
    sigma    <- imputed.Data$covMatrices[,,k]	
    B        <- sigma[ col.idx.Y,  col.idx.Y]
    C        <- sigma[ col.idx.Y, -col.idx.Y,drop=FALSE]
    D        <- sigma[-col.idx.Y, -col.idx.Y]
    CDinv    <- C%*%solve(D)
    cMu      <- apply(data.cmu, 1, function(x){c(mu[col.idx.Y] + CDinv %*% (x - mu[-col.idx.Y]))})
    cmeans[[k]] <- t(cMu)  #for each of 1000 rows, conditional for choice B, conditional for choice C this is also way there are 500 marks
    idx[[k]]    <- sapply(cmeans[[k]][dgp$dgp.Sample$ind.Miss.Y], 
                          function(x){
                            which.min(abs(x - cmeans[[k]][!dgp$dgp.Sample$ind.Miss.Y]))
                          }) 
    data2.cmu   <- dgp$dgp.Sample$rdf.Obs$y_1[!dgp$dgp.Sample$ind.Miss.Y]
    draws[[k]]  <- data2.cmu[idx[[k]]]
    imputed.Data$imputations[[k]]$y_1                            <- dgp$dgp.Sample$rdf.Obs$y_1
    imputed.Data$imputations[[k]]$y_1[dgp$dgp.Sample$ind.Miss.Y] <- draws[[k]]
  }
  imputed.Data$idx <- idx
  return(imputed.Data)
  
}

