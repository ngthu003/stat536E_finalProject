#' Joint Multivariate Normal MI Analysis
#'
#' fn.Joint.MVN.MI is used to analyze MI-simulated data via 
#' Joint Multivariate Normal MI approach w/ amelia packages
#'
#' @param dgp         # randomly generated dataset w/ missing data
#' @param m           #(imputation)
#' @param command     amelia or norm, both differ in modeling variance of imputed values


fn.Joint.MVN.MI <- function(
    dgp,
    m, 
    command = 'amelia'
) {
  
  
  # ======================================================!
  # 1) Meta Info ------------------------------------------
  # ------------------------------------------------------!
  var.Type    <- dgp$dgp.Settings$var.Type 
  imp.Method  <- dgp$dgp.Settings$imp.Method
  rdf.Obs     <- dgp$dgp.Sample$rdf.Obs
  no.Category <- dgp$dgp.Settings$no.Category 
  
  # Meta data on target var.
  v <- which(names(rdf.Obs) == 'y_1')
  if (var.Type == 'ordinal') {
    rdf.Obs$y_1 <- as.numeric(rdf.Obs$y_1)
  }
  
  
  
  # ======================================================!
  # 2) MI Data --------------------------------------------
  # ------------------------------------------------------!
  
  ## Impute Data ----------------------
  time.Start   <- proc.time()
  if (command == 'amelia') {
    imputed.Data <- amelia(rdf.Obs, m=m, p2s = 0)
  } else {
    imputed.Data             <- list()
    imputed.Data$imputations <- vector("list", m)
    imputed.Data$mu          <- matrix(, nrow = ncol(rdf.Obs), ncol=m)
    imputed.Data$covMatrices <- array(, dim = c(ncol(rdf.Obs), ncol(rdf.Obs), m))
    if (command == 'mlest'){
      mvn <- mlest(rdf.Obs)
      for(k in 1:m){
        imputed.Data$imputations[[k]]                 <- rdf.Obs
        draws                                         <- mvrnorm(n=nrow(rdf.Obs), mu=mvn$muhat, Sigma=mvn$sigmahat)
        imputed.Data$imputations[[k]][is.na(rdf.Obs)] <- draws[is.na(rdf.Obs)]
        imputed.Data$mu[,k]                           <- colMeans(imputed.Data$imputations[[k]])
        imputed.Data$covMatrices[,,k]                 <- cov(imputed.Data$imputations[[k]])
      }
    }
    if (command == 'norm') {
      for(k in 1:m){
        prelim                          <- prelim.norm(as.matrix(rdf.Obs))
        text                            <- capture.output(thetahat <- em.norm(prelim))
        rngseed(round(runif(1, min=0, max=10000000)))
        draws                           <- imp.norm(prelim, thetahat, rdf.Obs)
        imputed.Data$imputations[[k]]   <- data.frame(draws)
        imputed.Data$mu[,k]             <- colMeans(draws)
        imputed.Data$covMatrices[,,k]   <- cov(draws)
      }
    }
  }
  time.End     <- proc.time()
  time.to.MI   <- as.numeric((time.End - time.Start)[3])
  class(imputed.Data$imputations) <- "list"
  
  
  
  
  
  # ======================================================!
  # 3) Fix Categorical Var --------------------------------
  # ------------------------------------------------------!
  
  if (var.Type == 'ordinal') {
    res <- lapply(
      imputed.Data$imputations, 
      function(x) {
        p <- x$y_1
        # Normalize the probability
        if (var.Type == 'ordinal') {
          p <- (p - min(rdf.Obs$y_1, na.rm=T)) / (max(rdf.Obs$y_1, na.rm=T) - min(rdf.Obs$y_1, na.rm=T))
        }
        # Standardize the prob. to 
        p  <- p*(p>0)*(p<1) + (p>1)
        # Sample from Bino(p.vector)
        d  <- sapply(p, function(y){rbinom(1, (no.Category-1), y)+1})
        choice.Prob <- t(sapply(p, 
                                function(y){ dbinom(0:(no.Category-1), (no.Category-1), y)}, 
                                simplify=TRUE))
        res       <- list()
        res$draws <- d[dgp$dgp.Sample$ind.Miss.Y]
        res$choice.Prob    <- choice.Prob
        return(res)
      })
  }
  
  
  
  # ======================================================!
  # 4) Extra: PMM (instead of PPD) ------------------------
  # ------------------------------------------------------!
  
  if (imp.Method == 'pmm'){
    imputed.Data <- .amelia.pmm(dgp, imputed.Data, var.Type=var.Type, v=v, command=command)
    if (type == 'ordinal'){
      for (k in seq_along(res$choice.Prob)){
        res$choice.Prob[[k]][dgp$dgp.Sample$ind.Miss.Y,] <- res$choice.Prob[[k]][!dgp$dgp.Sample$ind.Miss.Y,][imputed.Data$mark,]
      }	
    }
  }	
  
  
  
  
  # ======================================================!
  # 5) Fix Ordering in Ordinal Var ------------------------
  # ------------------------------------------------------!
  if (var.Type == 'ordinal') {
    for (k in seq_along(imputed.Data$imputations)) {
      if (imp.Method == 'pmm') {
        imputed.Data$imputations[[k]]$y_1[dgp$dgp.Sample$ind.Miss.Y] <- toupper(letters[res[[k]]$draws]) 
      } else {
        imputed.Data$imputations[[k]]$y_1[dgp$dgp.Sample$ind.Miss.Y] <- res[[k]]$draws
        imputed.Data$imputations[[k]]$y_1                            <- ordered(imputed.Data$imputations[[k]]$y_1, 
                                                                                labels=toupper(letters[1:no.Category]))
      }
    }
  }
  
  
  
  
  # ======================================================!
  # 6) Acc: Imp. Values -----------------------------------
  # ------------------------------------------------------!
  if (var.Type == 'continuous') {
    ### ==============================!
    ### Continuous --------------------
    Imputed.Values.RMSE <- mean(sapply(imputed.Data$imputations, 
                                       function(d) {
                                         sqrt(mean(((d$y_1 - 
                                                       dgp$dgp.Sample$rdf.True.Y)^2)[dgp$dgp.Sample$ind.Miss.Y]))
                                       }))
    Imputed.Values.Bias <- mean(sapply(imputed.Data$imputations, 
                                       function(d) { 
                                         abs(mean(d$y_1[dgp$dgp.Sample$ind.Miss.Y]) - 
                                               mean(dgp$dgp.Sample$rdf.True.Y[dgp$dgp.Sample$ind.Miss.Y]))
                                       }))
    # --------------------------------!
  } else if (var.Type == 'ordinal') {
    ### ==============================!
    ### Ordinal -----------------------
    Imputed.Values.RMSE <- mean(sapply(imputed.Data$imputations, 
                                       function(d) {
                                         mean((d$y_1 == dgp$dgp.Sample$rdf.True.Y)[dgp$dgp.Sample$ind.Miss.Y])
                                       }))
    Imputed.Values.Bias <- NA	
    # --------------------------------!
  } else {
    ### ==============================!
    ### Else --------------------------
    stop('var.Type must be continuous or ordinal.')
  }
  
  
  
  # ======================================================!
  # 7) Acc: Imp. Choice Prob. -----------------------------
  # ------------------------------------------------------!
  if (var.Type == 'ordinal') {
    ### ==============================!
    ### Ordinal -----------------------
    for(k in seq_along(res)) {
      res[[k]]$choice.Prob <- res[[k]]$choice.Prob[dgp$dgp.Sample$ind.Miss.Y,][dgp$dgp.Sample$col]
    }	
    Choice.Prob.RMSE <- sqrt(mean(sapply(res, 
                                         function(x) {
                                           (x$choice.Prob - dgp$dgp.Models$choice.Prob)^2
                                         })))
    Choice.Prob.Bias <- abs(mean(sapply(res, 
                                        function(x) {
                                          abs(mean(x$choice.Prob)-mean(dgp$dgp.Models$choice.Prob))
                                        })))
  } 
  
  
  
  
  
  # ======================================================!
  # 8) Build Model w/ Observed Data Only ------------------
  # ------------------------------------------------------!
  
  ## f1: Target, X1, is fully observed --------------------
  ## Same for all var. Types
  tmp.mdl.X1 <- pool(dgp$dgp.Models$f1.X1, data=imputed.Data$imputations, m=m, FUN=bayesglm)
  
  
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
    Choice.Prob.RMSE <- Choice.Prob.Bias <- NA
    tmp.mdl.Y <- pool(dgp$dgp.Models$f2.Y, data=imputed.Data$imputations, m=m, FUN=bayesglm)
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
    tmp.mdl.Y <- pool(dgp$dgp.Models$f2.Y, data=imputed.Data$imputations, m=m, FUN=bayespolr)
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