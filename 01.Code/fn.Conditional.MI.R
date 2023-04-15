#' Conditional MI Analysis
#'
#' fn.Conditional.MI is used to analyze MI-simulated data via Conditional MI approach w/ mi package
#'
#' @param dgp         # randomly generated dataset w/ missing data
#' @param m           #(imputation)
#' @param imp.Method  imputation method: PPD or PMM
#' @param chains      for calling mi fn, default to 4
#' @param iter        #(iteration) per imputation, default to 30
#' @param no.Category #(category) per discrete var.
#' @param mcar        Indicator for if MCAR, default to FALSE


fn.Conditional.MI <- function(
    dgp,
    m, 
    imp.Method ="ppd", chains=5, iter=30, no.Category=3, mcar=FALSE
) {

  
  # ======================================================!
  # 1) Meta Info ------------------------------------------
  # ------------------------------------------------------!
  var.Type   <- dgp$dgp.Settings$var.Type 
  imp.Method <- dgp$dgp.Settings$imp.Method
  # Meta data on each var.
  mdf        <- missing_data.frame(dgp$dgp.Sample$rdf.Obs)
  # Check for MCAR
  if(!mcar){
    mdf@variables[["y_1"]]@imputation_method <- imp.Method
  } else {
    iter <- 0
  } 
  
  
  # ======================================================!
  # 2) MI Data --------------------------------------------
  # ------------------------------------------------------!
  
  ## Impute Data ----------------------
  time.Start   <- proc.time()
  imputed.Data <- mi(mdf, n.chains = chains, n.iter = iter, verbose = FALSE, save_models = TRUE); 
  time.End     <- proc.time()
  time.to.MI   <- as.numeric((time.End - time.Start)[3])
  # Get completely imputed data
  imputed.Data.Complete <- complete(imputed.Data, m)
  
  ## Choice Prob. ---------------------
  choice.Prob <- vector("list", chains)
  for(k in 1:chains){
    if (var.Type == 'ordinal') {
      choice.Prob[[k]] <- (imputed.Data@data[[k]]@variables$y_1@fitted)[dgp$dgp.Sample$ind.Miss.Y,]
    }
    if(var.Type!="continuous") {
      choice.Prob[[k]] <- choice.Prob[[k]][dgp$dgp.Sample$col]  
    }
  }
  # choice.Prob
  
  
  
  # ======================================================!
  # 3) Acc: Imp. Values -----------------------------------
  # 4) Acc: Imp. Choice Prob. -----------------------------
  # ------------------------------------------------------!
  
  if (var.Type == "continuous") {
    ### ==============================!
    ### Continuous --------------------
    Choice.Prob.RMSE <- NA
    Choice.Prob.Bias <- NA
    Imputed.Values.RMSE <- mean(sapply(imputed.Data.Complete, 
                                       function(d) { 
                                         sqrt(mean(((d$y_1 - dgp$dgp.Sample$rdf.True.Y)^2)[dgp$dgp.Sample$ind.Miss.Y])) 
                                       }))
    Imputed.Values.Bias <- mean(sapply(imputed.Data.Complete, 
                                       function(d) {
                                         abs(mean(d$y_1[dgp$dgp.Sample$ind.Miss.Y]) - 
                                               mean(dgp$dgp.Sample$rdf.True.Y[dgp$dgp.Sample$ind.Miss.Y]))
                                       }))
    
    # mvmatch      <- mean(sapply(complete, FUN = function(d){sqrt(mean(((d$y_1 - dgp$true)^2)[dgp$miss]))}))
    # mvmatch.bias <- mean(sapply(complete, FUN = function(d){abs(mean(d$y_1[dgp$miss]) - mean(dgp$true[dgp$miss]))}))
    # --------------------------------!
  } else if (var.Type == 'ordinal') {
    ### ==============================!
    ### Ordinal -----------------------
    Choice.Prob.RMSE    <- sqrt(mean(sapply(choice.Prob, 
                                            function(x){ (x - dgp$dgp.Models$choice.Prob)^2 })))
    Choice.Prob.Bias    <- abs( mean(sapply(choice.Prob, 
                                            function(x){ abs(mean(x)-mean(dgp$dgp.Models$choice.Prob)) })))
    Imputed.Values.RMSE <- mean(sapply(imputed.Data.Complete, 
                                       function(d){ 
                                         mean((d$y_1 ==
                                                 dgp$dgp.Sample$rdf.True.Y)[dgp$dgp.Sample$ind.Miss.Y]) 
                                       })) 
    Imputed.Values.Bias <- NA
    # --------------------------------!
  } else {
    ### ==============================!
    ### Else --------------------------
    stop('var.Type must be continuous or ordinal.')
  }
  
  
  
  # ======================================================!
  # 5) Build Model w/ Observed Data Only ------------------
  # ------------------------------------------------------!
  
  ## f1: Target, X1, is fully observed --------------------
  ## Same for all var. Types
  tmp.mdl.X1 <- pool(dgp$dgp.Models$f1.X1, data=imputed.Data, m=m, FUN=bayesglm)
  
  ### ====================================================!
  ### 6A) Acc: Coefficients -------------------------------
  mdl.X1.coef.RMSE <- sqrt(mean((tmp.mdl.X1@coefficients  - coef(dgp$dgp.Models$mdl.X1))^2   ))
  mdl.X1.coef.Bias <- abs( mean( tmp.mdl.X1@coefficients) - mean(coef(dgp$dgp.Models$mdl.X1) ))
  mdl.X1.coef.MNS  <- mahalanobis(x      = coef(dgp$dgp.Models$mdl.X1),
                                  center = tmp.mdl.X1@coefficients, 
                                  cov    = vcov(dgp$dgp.Models$mdl.X1))
  
  ### ====================================================!
  ### 6B) Acc: Fitted Values ------------------------------
  eta                    <- model.matrix(dgp$dgp.Models$mdl.X1) %*% tmp.mdl.X1@coefficients 
  mdl.X1.fitted.RMSE     <- sqrt(mean((eta[dgp$dgp.Sample$ind.Fully.Obs.XY] - 
                                         dgp$dgp.Models$mdl.X1.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY])^2))
  mdl.X1.fitted.Bias     <- abs(mean(eta[dgp$dgp.Sample$ind.Fully.Obs.XY]) - 
                                  mean(dgp$dgp.Models$mdl.X1.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY]))
  mdl.X1.fitted.All.RMSE <- sqrt(mean((eta - dgp$dgp.Models$mdl.X1.fitted)^2))
  mdl.X1.fitted.All.Bias <- abs( mean(eta) - mean(dgp$dgp.Models$mdl.X1.fitted) )
  
  
  
  ## f2: Target, Y, has missing data ----------------------
  if (var.Type == "continuous") {
    ### Continuous --------------------
    tmp.mdl.Y <- pool(dgp$dgp.Models$f2.Y, data=imputed.Data, m=m, FUN=bayesglm)
    #### ===================================================!
    #### 6A) Acc: Coefficients ------------------------------
    mdl.Y.coef.RMSE <- sqrt(mean((tmp.mdl.Y@coefficients  - coef(dgp$dgp.Models$mdl.Y))^2  ))
    mdl.Y.coef.Bias <- abs( mean( tmp.mdl.Y@coefficients) - mean(coef(dgp$dgp.Models$mdl.Y)))
    mdl.Y.coef.MNS  <- mahalanobis(x      = coef(dgp$dgp.Models$mdl.Y), 
                                   center = tmp.mdl.Y@coefficients, 
                                   cov    = vcov(dgp$dgp.Models$mdl.Y))
    
    #### ===================================================!
    #### 6B) Acc: Fitted Values -----------------------------
    eta                   <- model.matrix(dgp$dgp.Models$mdl.Y) %*% tmp.mdl.Y@coefficients 
    mdl.Y.fitted.RMSE     <- sqrt(mean((eta[dgp$dgp.Sample$ind.Fully.Obs.XY] - 
                                          dgp$dgp.Models$mdl.Y.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY])^2))
    mdl.Y.fitted.Bias     <- abs(mean(eta[dgp$dgp.Sample$ind.Fully.Obs.XY]) - 
                                   mean(dgp$dgp.Models$mdl.Y.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY]))	
    mdl.Y.fitted.All.RMSE <- sqrt(mean((eta - dgp$dgp.Models$mdl.Y.fitted)^2))
    mdl.Y.fitted.All.Bias <- abs( mean(eta) - mean(dgp$dgp.Models$mdl.Y.fitted))
    
    # --------------------------------!
  } else if (var.Type == 'ordinal') {
    ### Ordinal -----------------------
    tmp.mdl.Y <- pool(dgp$dgp.Models$f2.Y, data=imputed.Data, m=m, FUN=bayespolr)
    coef      <- tmp.mdl.Y@coefficients[1:(dgp$dgp.Settings$no.Var.Total-1)]
    zeta      <- tmp.mdl.Y@coefficients[(dgp$dgp.Settings$no.Var.Total):length(tmp.mdl.Y@coefficients)]
    
    #### ===================================================!
    #### 6A) Acc: Coefficients ------------------------------
    mdl.Y.coef.RMSE <- sqrt(mean((coef - coef(dgp$dgp.Models$mdl.Y))^2)) 
    mdl.Y.coef.Bias <- abs(mean(coef) - mean(coef(dgp$dgp.Models$mdl.Y)))
    mdl.Y.coef.MNS  <- mahalanobis(x      = coef(dgp$dgp.Models$mdl.Y)[1:(dgp$dgp.Settings$no.Var.Total-1)], 
                                   center = tmp.mdl.Y@coefficients[1:(dgp$dgp.Settings$no.Var.Total-1)], 
                                   cov    = vcov(dgp$dgp.Models$mdl.Y)[1:(dgp$dgp.Settings$no.Var.Total-1),
                                                                       1:(dgp$dgp.Settings$no.Var.Total-1)])
    
    #### ===================================================!
    #### 6B) Acc: Fitted Values -----------------------------
    eta                   <- sapply(zeta, 
                                    function(x) {x - model.matrix(dgp$dgp.Models$mdl.Y)[,-1]%*%coef}, 
                                    simplify=TRUE)
    mdl.Y.fitted.RMSE     <- sqrt(mean((eta[dgp$dgp.Sample$ind.Fully.Obs.XY,] - 
                                          dgp$dgp.Models$mdl.Y.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY,])^2))
    mdl.Y.fitted.Bias     <- abs(mean(eta[dgp$dgp.Sample$ind.Fully.Obs.XY,]) - 
                                   mean(dgp$dgp.Models$mdl.Y.fitted[dgp$dgp.Sample$ind.Fully.Obs.XY,]))	
    mdl.Y.fitted.All.RMSE <- sqrt(mean((eta - dgp$dgp.Models$mdl.Y.fitted)^2))
    mdl.Y.fitted.All.Bias <- abs( mean(eta) - mean(dgp$dgp.Models$mdl.Y.fitted))
    
    # --------------------------------!
  } else {
    ### ==============================!
    ### Else --------------------------
    stop('var.Type must be continuous or ordinal.')
  }
  
  
  
  # ======================================================!
  # 99) Return Results ------------------------------------
  # ------------------------------------------------------!
  results       <- list()
  results$time  <- time.to.MI
  # Acc. Imputed Values
  results$Imputed.Values.RMSE <- Imputed.Values.RMSE
  results$Imputed.Values.Bias <- Imputed.Values.Bias
  # Acc. Choice Prob.
  results$Choice.Prob.RMSE <- Choice.Prob.RMSE
  results$Choice.Prob.Bias <- Choice.Prob.Bias
  # Acc. Coefficients
  results$mdl.X1.coef.RMSE <- mdl.X1.coef.RMSE     
  results$mdl.X1.coef.Bias <- mdl.X1.coef.Bias
  results$mdl.X1.coef.MNS  <- mdl.X1.coef.MNS
  results$mdl.Y.coef.RMSE  <- mdl.Y.coef.RMSE     
  results$mdl.Y.coef.Bias  <- mdl.Y.coef.Bias
  results$mdl.Y.coef.MNS   <- mdl.Y.coef.MNS
  # Acc. Fitted Values
  results$mdl.X1.fitted.RMSE     <- mdl.X1.fitted.RMSE     
  results$mdl.X1.fitted.Bias     <- mdl.X1.fitted.Bias
  results$mdl.X1.fitted.All.RMSE <- mdl.X1.fitted.All.RMSE
  results$mdl.X1.fitted.All.Bias <- mdl.X1.fitted.All.Bias
  results$mdl.Y.fitted.RMSE      <- mdl.Y.fitted.RMSE     
  results$mdl.Y.fitted.Bias      <- mdl.Y.fitted.Bias
  results$mdl.Y.fitted.All.RMSE  <- mdl.Y.fitted.All.RMSE
  results$mdl.Y.fitted.All.Bias  <- mdl.Y.fitted.All.Bias
  
  
  return(results)

}
  