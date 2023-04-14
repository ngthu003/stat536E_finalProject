
##############################################################################|
#
# Stat 536E - Final Project
# - Author: Thu Nguyen
# - Topic:  MI - Join MVN vs. Condition
#
##############################################################################|



# ============================================================================!
# 0) Libraries ----------------------------------------------------------------
# ----------------------------------------------------------------------------!
libs <- c('tidyverse', 'mi', 'arm', 'Amelia', 'mvnmle', 'MASS')
lapply(libs, require, character.only = TRUE)
rm(libs)






# ============================================================================!
# 1) Data Generation ----------------------------------------------------------
# ----------------------------------------------------------------------------!


## 1.1) Sample settings -----------------------------------

### A) Fixed --------------------------
N               <- 10          # #(obs)
no.Var.Observed <- 3            # #(var. fully observed)

### B) Adjustable ---------------------
var.Type        <- 'continuous' # Type of variables
no.Var.Missing  <- 2            # #(var. w/ missing data)
miss.Rate.X     <- .1           # missing data rate per X explanatory var.
miss.Rate.Y     <- .25          # missing data rate per Y response var.
restriction     <- 'triangular' # missing data mechanism
strong <- 0                     # applicable only when "triangular"

imp.Method      <- 'ppd'

### C) Other --------------------------
estimate.CPCs <- F              # F to save time, default is T




## 1.2) Sample data ---------------------------------------
rdf.Colnames <- c(paste("x_", 1:(no.Var.Observed + no.Var.Missing), sep=""), "y_1")

rdf <- rdata.frame(
  N             = N,
  n_full        = no.Var.Observed,
  n_partial     = no.Var.Missing + 1,
  restrictions  = restriction,
  types         = c(rep('continuous', no.Var.Observed + no.Var.Missing), var.Type),
  pr_miss       = c(rep(miss.Rate.X, no.Var.Missing), miss.Rate.Y),
  strong        = strong,
  estimate_CPCs = estimate.CPCs
)

## Extract true & observed dfs
rdf.True <- rdf$true
rdf.Obs  <- rdf$obs
colnames(rdf.Obs) <- colnames(rdf.True) <- rdf.Colnames
rdf.Obs

## Row indicators for fully observed rows and row w/ missing Y
ind.Fully.Obs.XY <- !apply(rdf.Obs, 1, function(x) any(is.na(x))) 
ind.Miss.Y       <- is.na(rdf.Obs[,ncol(rdf.Obs)])

## Other Stuff
cpc <- (restrict=="none")
# if(cpc) nmar <- sqrt(mean(rdf$empirical_CPCs^2)) else nmar <- NA
nmar <- sqrt(mean(rdf$empirical_CPCs^2))




# ============================================================================!
# 1) Models -------------------------------------------------------------------
# ----------------------------------------------------------------------------!

f1.X1 <- as.formula(
  paste("x_1 ~ y_1", 
        paste(paste("x", 2:(no.Var.Observed + no.Var.Missing), sep="_"), 
              collapse="+"), sep="+"))
f2.Y <- as.formula(
  paste("y_1", 
        paste(paste("x", 1:(no.Var.Observed + no.Var.Missing), sep = "_"), 
              collapse = " + "), sep = " ~ "))



mdl.X1 <- bayesglm(f1.X1, data=rdf.True)
mdl.X1.fitted <- fitted(mdl.X1)


if (var.Type == "continuous") {
  mdl.Y <- bayesglm(f2.Y, data=rdf.True)
  mdl.Y.fitted <- fitted(mdl.Y) 	
  ev <- eigen(vcov(mdl.Y), symmetric = TRUE)	
  params <- (coef(mdl.Y) + (ev$vectors %*% (sqrt(ev$values) * rnorm(length(coef(mdl.Y)))))[,1])
  
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
} 

if (var.Type!="continuous") {
  col <- cbind(1:sum(miss), 
               sapply(true[["y_1"]][miss], function(x) {
                 which(x==levels(true[["y_1"]][miss]))
               }
              ))
}


dgp <- list()
# ------------------------------------!
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
dgp.Settings$no.Var.Total <- no.Var.Observed + no.Var.Missing +1
dgp.Settings$var.Type     <- var.Type
dgp.Settings$imp.Method   <- imp.Method
dgp.Settings$nmar         <- nmar
# --

dgp <- list(
  dgp.Sample = dgp.Sample,
  dgp.Models = dgp.Models,
  dgp.Settings = dgp.Settings
)

dgp



