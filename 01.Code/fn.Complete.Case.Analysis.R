#' Complete Case Analysis
#'
#' fn.Complete.Case.Analysis is used to analyze MI-simulated data via complete case analysis
#'
#' @param dgp # randomly generated dataset w/ missing data


fn.Complete.Case.Analysis <- function(dgp) {
  
  # ======================================================!
  # 1) Meta Info ------------------------------------------
  # ------------------------------------------------------!
  var.Type   <- dgp$dgp.Settings$var.Type 
  imp.Method <- dgp$dgp.Settings$imp.Method
  
  # ======================================================!
  # 2) Build Model w/ Observed Data Only ------------------
  # ------------------------------------------------------!
  
  ## f1: Target, X1, is fully observed --------------------
  ## Same for all var. Types
  tmp.mdl.X1 <- bayesglm(formula = dgp$dgp.Models$f1.X1, data = dgp$dgp.Sample$rdf.Obs)
  tmp.mdl.X1.fitted <- fitted(tmp.mdl.X1)
  
  ## f2: Target, Y, has missing data ----------------------
  if (var.Type == "continuous") {
    ### Continuous --------------------
    tmp.mdl.Y <- bayesglm(formula = dgp$dgp.Models$f2.Y, data = dgp$dgp.Sample$rdf.Obs)
    # --------------------------------!
  } else if (var.Type == 'ordinal') {
    ### Ordinal -----------------------
    tmp.mdl.Y  <- bayespolr(formula = dgp$dgp.Models$f2.Y, data = dgp$dgp.Sample$rdf.Obs, drop.unused.levels=FALSE)
    # --------------------------------!
  } else {
    ### ==============================!
    ### Else --------------------------
    stop('var.Type must be continuous or ordinal.')
  }
  
  
  # ======================================================!
  # 3) Acc: Coefficients ----------------------------------
  # ------------------------------------------------------!
  
  ## f1: X1 ----------------------------
  mdl.X1.coef.RMSE <- sqrt(mean((coef(tmp.mdl.X1) - coef(dgp$dgp.Models$mdl.X1))^2))
  mdl.X1.coef.Bias <- abs(mean(coef(tmp.mdl.X1)) - mean(coef(dgp$dgp.Models$mdl.X1)))
  mdl.X1.coef.MNS  <- mahalanobis(x      = coef(dgp$dgp.Models$mdl.X1), 
                                  center = coef(tmp.mdl.X1), 
                                  cov    = vcov(dgp$dgp.Models$mdl.X1))
  ## f2: Y -----------------------------
  mdl.Y.coef.RMSE <- sqrt(mean((coef(tmp.mdl.Y) - coef(dgp$dgp.Models$mdl.Y))^2))
  mdl.Y.coef.Bias <- abs(mean(coef(tmp.mdl.Y)) - mean(coef(dgp$dgp.Models$mdl.Y)))
  mdl.Y.coef.MNS  <- mahalanobis(x      = coef(dgp$dgp.Models$mdl.Y)[1:(dgp$dgp.Settings$no.Var.Total-1)], 
                                 center = coef(tmp.mdl.Y)[1:(dgp$dgp.Settings$no.Var.Total-1)], 
                                 cov    = vcov(dgp$dgp.Models$mdl.Y)[1:(dgp$dgp.Settings$no.Var.Total-1),
                                                                     1:(dgp$dgp.Settings$no.Var.Total-1)])
  
  
  # ======================================================!
  # 4) Acc: Fitted Values ---------------------------------
  # ------------------------------------------------------!
  
  ## f1: X1 ----------------------------
  mdl.X1.fitted.RMSE <- sqrt(mean((tmp.mdl.X1.fitted - 
                                     dgp$dgp.Models$mdl.X1.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY])^2))
  mdl.X1.fitted.Bias <- abs(mean(tmp.mdl.X1.fitted) - 
                              mean(dgp$dgp.Models$mdl.X1.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY]))
  
  ## f2: Y -----------------------------
  if (var.Type == "continuous") {
    ### ==============================!
    ### Continuous --------------------
    mdl.Y.fitted.RMSE <- sqrt(mean((fitted(tmp.mdl.Y)  - 
                                      dgp$dgp.Models$mdl.Y.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY])^2))
    mdl.Y.fitted.Bias <- abs(  mean(fitted(tmp.mdl.Y)) - 
                                 mean(dgp$dgp.Models$mdl.Y.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY]))	
    # --------------------------------!
  } else if (var.Type == 'ordinal') {
    ### ==============================!
    ### Ordinal -----------------------
    ld                <- sapply(tmp.mdl.Y$zeta, 
                                function(x){x - model.matrix(tmp.mdl.Y)[,-1]%*%coef(tmp.mdl.Y)}, 
                                simplify=TRUE)
    mdl.Y.fitted.RMSE <- sqrt(mean((ld - 
                                      dgp$dgp.Models$mdl.Y.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY])^2))	
    mdl.Y.fitted.Bias <- abs( mean(ld[dgp$dgp.Sample$ind.Miss.Y]) - 
                                mean(dgp$dgp.Models$mdl.Y.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY]))		
    # --------------------------------!
  } else {
    ### ==============================!
    ### Else --------------------------
    stop('var.Type must be continuous or ordinal.')
  }
  
  
  # ======================================================!
  # 99) Return Results ------------------------------------
  # ------------------------------------------------------!
  results <- list()
  # Acc. Coefficients
  results$mdl.X1.coef.RMSE <- mdl.X1.coef.RMSE     
  results$mdl.X1.coef.Bias <- mdl.X1.coef.Bias
  results$mdl.X1.coef.MNS  <- mdl.X1.coef.MNS
  results$mdl.Y.coef.RMSE <- mdl.Y.coef.RMSE     
  results$mdl.Y.coef.Bias <- mdl.Y.coef.Bias
  results$mdl.Y.coef.MNS  <- mdl.Y.coef.MNS
  # Acc. Fitted Values
  results$mdl.X1.fitted.RMSE <- mdl.X1.fitted.RMSE     
  results$mdl.X1.fitted.Bias <- mdl.X1.fitted.Bias
  results$mdl.Y.fitted.RMSE <- mdl.Y.fitted.RMSE     
  results$mdl.Y.fitted.Bias <- mdl.Y.fitted.Bias
  
  return(results)
  
}