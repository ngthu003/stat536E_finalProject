#' (Missing) Data Generation Process
#'
#' fn.dgp is used to generate random (w/ missing) data using rdata.frame from mi
#'
#' A) Fixed --------------------------
#' @param N               # #(obs)
#' @param no.Var.Observed # #(var. fully observed)
#' B) Adjustable ---------------------
#' @param var.Type        # Type of variables
#' @param no.Category     # #(category for discrete var.)
#' @param no.Var.Missing  # #(var. w/ missing data)
#' @param miss.Rate.X     # missing data rate per X explanatory var.
#' @param miss.Rate.Y     # missing data rate per Y response var.
#' @param restriction     # missing data mechanism
#' @param strong          # applicable only when "triangular"
#' 
#' @param imp.Method      # imputation method, PPD or PMM
#' 
#' C) Other --------------------------
#' @param estimate.CPCs   # F to save time, default is T
#'
#' @examples
#' ### A) Fixed --------------------------
#' N               <- 10           # #(obs)
#' no.Var.Observed <- 3            # #(var. fully observed)
#' 
#' ### B) Adjustable ---------------------
#' var.Type        <- 'continuous' # Type of variables
#' no.Var.Missing  <- 2            # #(var. w/ missing data)
#' miss.Rate.X     <- .1           # missing data rate per X explanatory var.
#' miss.Rate.Y     <- .25          # missing data rate per Y response var.
#' restriction     <- 'triangular' # missing data mechanism
#' strong <- 0                     # applicable only when "triangular"
#' 
#' no.Category     <- 3
#' imp.Method      <- 'ppd'
#'
#' ### C) Other --------------------------
#' estimate.CPCs <- F              # F to save time, default is T


fn.DGP <- function(N, 
                   var.Type,
                   no.Category = 0,
                   no.Var.Observed, 
                   no.Var.Missing,
                   miss.Rate.X, 
                   miss.Rate.Y,
                   restriction, 
                   strong = 0,
                   imp.Method = 'ppd') {
  
  estimate.CPCs <- (restriction == 'none')
  
  # ======================================================!
  # 1) Generate Data --------------------------------------
  # ------------------------------------------------------!
  rdf.Colnames <- c(paste("x_", 1:(no.Var.Observed + no.Var.Missing), sep=""), "y_1")
  
  if (var.Type == 'continuous') {
    ## ===============================!
    ## Continuous ---------------------
    rdf <- rdata.frame(
      N             = N,
      n_full        = no.Var.Observed,
      n_partial     = no.Var.Missing + 1,
      restrictions  = restriction,
      types         = c(rep('continuous', no.Var.Observed + no.Var.Missing), var.Type),
      pr_miss       = c(rep(miss.Rate.X, no.Var.Missing), miss.Rate.Y),
      strong        = strong,
      estimate_CPCs = estimate.CPCs
    ) # ------------------------------!
  } else if (var.Type == 'ordinal') {
    ## ===============================!
    ## Ordinal ------------------------
    rdf <- rdata.frame(
      N             = N,
      n_full        = no.Var.Observed,
      n_partial     = no.Var.Missing + 1,
      restrictions  = restriction,
      types         = c(rep('continuous', no.Var.Observed + no.Var.Missing), var.Type),
      pr_miss       = c(rep(miss.Rate.X, no.Var.Missing), miss.Rate.Y),
      n_cat         = no.Category,
      strong        = strong,
      estimate_CPCs = estimate.CPCs
    ) # ------------------------------!
  } else {
    ## ===============================!
    ## Else ---------------------------
    stop('var.Type must be continuous or ordinal.')
  }
  
  
  ## Extract true & observed dfs
  rdf.True <- rdf$true
  rdf.Obs  <- rdf$obs
  colnames(rdf.Obs) <- colnames(rdf.True) <- rdf.Colnames
  rdf.Obs
  
  ## Row indicators for fully observed rows and row w/ missing Y
  ind.Fully.Obs.XY <- !apply(rdf.Obs, 1, function(x) any(is.na(x))) 
  ind.Miss.Y       <- is.na(rdf.Obs[,ncol(rdf.Obs)])
  
  ## Other Stuff
  # cpc <- (restriction == "none")
  # if(cpc) nmar <- sqrt(mean(rdf$empirical_CPCs^2)) else nmar <- NA
  nmar <- sqrt(mean(rdf$empirical_CPCs^2))
  
  ## Other stuff
  if (var.Type != "continuous") {
    col <- cbind(1:sum(ind.Miss.Y), 
                 sapply(rdf.True[["y_1"]][ind.Miss.Y], 
                        function(x) { which(x==levels(rdf.True[["y_1"]][ind.Miss.Y])) } ))
  }
  
  
  
  # ======================================================!
  # 2) Build True Models ----------------------------------
  # ------------------------------------------------------!
  
  ## f1: Target, X1, is fully observed --------------------
  ## Same for all var. Types
  f1.X1 <- as.formula(
    paste("x_1 ~ y_1", 
          paste(paste("x", 2:(no.Var.Observed + no.Var.Missing), sep="_"), 
                collapse="+"), sep="+"))
  mdl.X1 <- bayesglm(f1.X1, data=rdf.True)
  mdl.X1.fitted <- fitted(mdl.X1)
  
  ## f2: Target, Y, has missing data ----------------------
  f2.Y <- as.formula(
    paste("y_1", 
          paste(paste("x", 1:(no.Var.Observed + no.Var.Missing), sep = "_"), 
                collapse = " + "), sep = " ~ "))
  
  if (var.Type == "continuous") {
    ### ==============================!
    ### Continuous --------------------
    mdl.Y <- bayesglm(f2.Y, data=rdf.True)
    mdl.Y.fitted <- fitted(mdl.Y) 	
    ev <- eigen(vcov(mdl.Y), symmetric = TRUE)	
    params <- (coef(mdl.Y) + (ev$vectors %*% (sqrt(ev$values) * rnorm(length(coef(mdl.Y)))))[,1])
    
    # Imputation Method ----------------------------------------------!
    if (imp.Method == 'ppd') {
      # Posterior Predictive Distribution ----------------------------!
      draws              <- model.matrix(mdl.Y)%*%params
      acc.Imputed.Values <- sqrt(mean((draws[ind.Miss.Y] - rdf.True$y_1[ind.Miss.Y])^2))
    } else if (imp.Method == 'pmm') {
      # Predictive Mean Matching -------------------------------------!
      idx                <- sapply((model.matrix(mdl.Y)%*%params)[ind.Miss.Y], 
                                   function(m) which.min(abs(m - rdf.True$y_1[!ind.Miss.Y])) )
      acc.Imputed.Values <- sqrt(mean((rdf.True$y_1[ind.Miss.Y] - rdf.True$y_1[!ind.Miss.Y][idx])^2))
    } else {
      stop('Imputation method must be either ppd or pmm.')
    }
    
    choice.Prob <- NA
    
    # --------------------------------!
  } else if (var.Type == 'ordinal') {
    ### ==============================!
    ### Ordinal -----------------------
    mdl.Y        <- bayespolr(f2.Y, data=rdf.True, drop.unused.levels=FALSE)
    mdl.Y.fitted <- sapply(mdl.Y$zeta, 
                           function(x) { x - model.matrix(mdl.Y)[,-1]%*%coef(mdl.Y) }, 
                           simplify=TRUE)
    ev     <- eigen(vcov(mdl.Y), symmetric = TRUE)
    params <- (c(coef(mdl.Y), 
                 mdl.Y$zeta) + (ev$vectors %*% (sqrt(ev$values) * rnorm(length(c(coef(mdl.Y), mdl.Y$zeta)))))[,1])
    coef   <- params[1:(length(params)-no.Category+1)] 
    zeta   <- params[  (length(params)-no.Category+2):length(params)]
    eta    <- sapply(zeta, 
                     function(x) { x - model.matrix(mdl.Y)[,-1]%*%coef }, 
                     simplify=TRUE)
    eta    <- cbind(0, plogis(eta), 1)
    
    choice.Prob <- t(diff(t(eta)))
    
    # Imputation Method ----------------------------------------------!
    if (imp.Method == 'ppd') {
      # Posterior Predictive Distribution ----------------------------!
      draws              <- apply(choice.Prob, 1, function(p) which(rmultinom(1, 1, p) == 1))
      draws              <- factor(draws, ordered=TRUE)
      acc.Imputed.Values <- mean(factor(toupper(letters[draws]), ordered=TRUE)[ind.Miss.Y] == rdf.True$y_1[ind.Miss.Y])
      choice.Prob        <- choice.Prob[ind.Miss.Y,]
    } else if (imp.Method == 'pmm') {
      # Predictive Mean Matching -------------------------------------!
      idx                <- sapply((model.matrix(mdl.Y)[,-1]%*%coef)[ind.Miss.Y], 
                                   function(m) {
                                     which.min(abs(m - (model.matrix(mdl.Y)[,-1]%*%coef)[!ind.Miss.Y]))
                                   })
      acc.Imputed.Values <- mean(rdf.True$y_1[ind.Miss.Y] == rdf.True$y_1[!ind.Miss.Y][idx])
      choice.Prob        <- cbind((1-t(diff(t(eta)))), 
                                  t(diff(t(eta))))[!ind.Miss.Y,][idx,]			
    } else {
      stop('Imputation method must be either ppd or pmm.')
    }
    
    choice.Prob <- choice.Prob[col]
    
    # --------------------------------!
  } else {
    ### ==============================!
    ### Else --------------------------
    stop('var.Type must be continuous or ordinal.')
  }
  

  
  
  # ======================================================!
  # 3) Return Everything ----------------------------------
  # ------------------------------------------------------!
  dgp.Sample <- list()
  # --
  dgp.Sample$rdf              <- rdf
  # --
  dgp.Sample$rdf.True.Y       <- rdf.True$y_1
  dgp.Sample$rdf.True         <- rdf.True
  dgp.Sample$rdf.Obs          <- rdf.Obs
  # --
  dgp.Sample$ind.Miss.Y       <- ind.Miss.Y
  dgp.Sample$ind.Fully.Obs.XY <- rowSums(is.na(rdf.Obs))==0
  # --
  dgp.Sample$col              <- col
  # ------------------------------------!
  dgp.Models <- list()
  # --
  dgp.Models$f1.X1         <- f1.X1
  dgp.Models$mdl.X1        <- mdl.X1
  dgp.Models$mdl.X1.fitted <- mdl.X1.fitted
  # --
  dgp.Models$f2.Y          <- f2.Y
  dgp.Models$mdl.Y         <- mdl.Y
  dgp.Models$mdl.Y.fitted  <- mdl.Y.fitted
  # --
  dgp.Models$acc.Imputed.Values <- acc.Imputed.Values
  dgp.Models$choice.Prob   <- choice.Prob
  # ------------------------------------!
  dgp.Settings <- list()
  # --
  dgp.Settings$no.Var.Observed <- no.Var.Observed
  dgp.Settings$no.Var.Missing  <- no.Var.Missing
  dgp.Settings$no.Var.Total    <- no.Var.Observed + no.Var.Missing +1
  dgp.Settings$var.Type        <- var.Type
  dgp.Settings$no.Category     <- no.Category
  dgp.Settings$imp.Method      <- imp.Method
  dgp.Settings$nmar            <- nmar
  # ------------------------------------!
  # ------------------------------------!
  dgp <- list(
    dgp.Sample   = dgp.Sample,
    dgp.Models   = dgp.Models,
    dgp.Settings = dgp.Settings
  )
  
  return(dgp)
}
