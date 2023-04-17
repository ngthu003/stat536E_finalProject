
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
libs <- c('tidyverse', 'mi', 'arm', 'Amelia', 'mvnmle', 'MASS', 'norm2', 'nnet')
lapply(libs, require, character.only = TRUE)
rm(libs)



# ============================================================================!
# 1) DGP ----------------------------------------------------------------------
# ----------------------------------------------------------------------------!

source('01.Code/fn.DGP.R')
#' @examples
#' ### A) Fixed --------------------------
#' tmp.N               <- 100           # #(obs)
#' tmp.no.Var.Observed <- 3            # #(var. fully observed)
#' 
#' ### B) Adjustable ---------------------
#' tmp.var.Type        <- 'continuous' # Type of variables
#' tmp.no.Var.Missing  <- 2            # #(var. w/ missing data)
#' tmp.miss.Rate.X     <- .1           # missing data rate per X explanatory var.
#' tmp.miss.Rate.Y     <- .25          # missing data rate per Y response var.
#' tmp.restriction     <- 'triangular' # missing data mechanism
#' tmp.strong <- 0                     # applicable only when "triangular"
#' 
#' tmp.no.Category     <- 3
#' tmp.imp.Method      <- 'pmm'
#' tmp.no.Imputation   <- 5
#'
#' ### C) Other --------------------------
#' tmp.estimate.CPCs <- F              # F to save time, default is T

## Continuous ---------------------------------------------
dgp.cts <- fn.DGP(
  N               = tmp.N,
  var.Type        = tmp.var.Type,
  no.Category     = tmp.no.Category,
  no.Var.Observed = tmp.no.Var.Observed,
  no.Var.Missing  = tmp.no.Var.Missing,
  miss.Rate.X     = tmp.miss.Rate.X,
  miss.Rate.Y     = tmp.miss.Rate.Y,
  restriction     = tmp.restriction,
  strong          = tmp.strong,
  imp.Method      = tmp.imp.Method)



## Ordinal ------------------------------------------------
tmp.var.Type = 'ordinal'
dgp.ord <- fn.DGP(
  N               = tmp.N,
  var.Type        = tmp.var.Type,
  no.Category     = tmp.no.Category,
  no.Var.Observed = tmp.no.Var.Observed,
  no.Var.Missing  = tmp.no.Var.Missing,
  miss.Rate.X     = tmp.miss.Rate.X,
  miss.Rate.Y     = tmp.miss.Rate.Y,
  restriction     = tmp.restriction,
  strong          = tmp.strong,
  imp.Method      = tmp.imp.Method)

rm(list = setdiff(ls(), c('fn.DGP', 'dgp.cts', 'dgp.ord')))




# ============================================================================!
# 2) Complete Case Analysis ---------------------------------------------------
# ----------------------------------------------------------------------------!

source('01.Code/fn.Complete.Case.Analysis.R')

cca.cts <- fn.Complete.Case.Analysis(dgp.cts)
data.frame(t(sapply(cca.cts,c)))
cca.ord <- fn.Complete.Case.Analysis(dgp.ord)
data.frame(t(sapply(cca.ord,c)))






# ============================================================================!
# 3) Conditional MI Analysis --------------------------------------------------
# ----------------------------------------------------------------------------!

source('01.Code/fn.Conditional.MI.R')

tmp.no.Imputation <- 5
cond.MI.cts <- fn.Conditional.MI(dgp.cts, m = tmp.no.Imputation)
data.frame(t(sapply(cond.MI.cts,c)))
cond.MI.ord <- fn.Conditional.MI(dgp.ord, m = tmp.no.Imputation, no.Category = dgp.ord$dgp.Settings$no.Category)
data.frame(t(sapply(cond.MI.ord,c)))





# ============================================================================!
# 4) Joint MVN MI Analysis ----------------------------------------------------
# ----------------------------------------------------------------------------!

source('01.Code/fn.Joint.MVN.MI.helper.PMM.R')
source('01.Code/fn.Joint.MVN.MI.R')


tmp.no.Imputation <- 5
joint.MI.amelia.cts <- fn.Joint.MVN.MI(dgp.cts, m = tmp.no.Imputation, command = 'amelia')
data.frame(t(sapply(joint.MI.amelia.cts,c)))
joint.MI.amelia.ord <- fn.Joint.MVN.MI(dgp.ord, m = tmp.no.Imputation, command = 'amelia')
data.frame(t(sapply(joint.MI.amelia.ord,c)))

joint.MI.amelia.cts <- fn.Joint.MVN.MI(dgp.cts, m = tmp.no.Imputation, command = 'norm')
data.frame(t(sapply(joint.MI.amelia.cts,c)))
joint.MI.amelia.ord <- fn.Joint.MVN.MI(dgp.ord, m = tmp.no.Imputation, command = 'norm')
data.frame(t(sapply(joint.MI.amelia.ord,c)))






# ============================================================================!
# 5) Overall Results ----------------------------------------------------------
# ----------------------------------------------------------------------------!

## Initialize df ------------------------------------------
overall.Results <- data.frame(
  Iteration = integer()
  # Adjustable hyper-parameters
  , restriction    = character()
  , var.Type       = character()
  , no.Var.Missing = integer()
  , miss.Rate.Y    = numeric()
  , miss.Rate.X    = numeric()
  # Fied hyper-parameters
  , N               = integer()
  , no.Var.Observed = integer()
  , no.Category     = integer()
  , no.Imputation   = integer()
  # Analysis method
  , analysis.Method = character()
  , imp.Method      = character()
  # Results
  # -- Acc: Imputed values
  , Imputed.Values.RMSE = numeric()
  , Imputed.Values.Bias = numeric()
  # -- Acc: Choice probabilities
  , Choice.Prob.RMSE = numeric()
  , Choice.Prob.Bias = numeric()
  # -- Acc: Coefficients
  , mdl.X1.coef.RMSE = numeric()
  , mdl.X1.coef.Bias = numeric()
  , mdl.X1.coef.MNS  = numeric()
  , mdl.Y.coef.RMSE  = numeric()
  , mdl.Y.coef.Bias  = numeric()
  , mdl.Y.coef.MNS   = numeric()
  # -- Acc: Fitted values
  , mdl.X1.fitted.RMSE     = numeric()
  , mdl.X1.fitted.Bias     = numeric()
  , mdl.X1.fitted.All.RMSE = numeric()
  , mdl.X1.fitted.All.Bias = numeric()
  , mdl.Y.fitted.RMSE      = numeric()
  , mdl.Y.fitted.Bias      = numeric()
  , mdl.Y.fitted.All.RMSE  = numeric()
  , mdl.Y.fitted.All.Bias  = numeric()
  # -- Time to Convergence
  , time.to.MI = numeric()
)



fn.null.to.NA <- function(x) {ifelse(is.null(x), NA, x)}


## Small Simulations --------------------------------------
iter <- 1

Options.analysis.Method <- c(
  'Complete_Case_Analysis',
  'Joint_MVN_MI_Amelia',
  'Joint_MVN_MI_Norm',
  'Conditional_MI'
)

### 1: var.Type -------------------------------------------
for (tmp.var.Type in c('continuous', 'ordinal')) {
  ### 2: imp.Method ---------------------------------------
  for (tmp.imp.Method in c('ppd', 'pmm')) {
    
    #### Generate data ====================================
    dgp <- fn.DGP(
      N               = tmp.N,
      var.Type        = tmp.var.Type,
      no.Category     = tmp.no.Category,
      no.Var.Observed = tmp.no.Var.Observed,
      no.Var.Missing  = tmp.no.Var.Missing,
      miss.Rate.X     = tmp.miss.Rate.X,
      miss.Rate.Y     = tmp.miss.Rate.Y,
      restriction     = tmp.restriction,
      strong          = tmp.strong,
      imp.Method      = tmp.imp.Method)
    
    ## 3: analysis.Method ---------------------------------
    for (tmp.analysis.Method in Options.analysis.Method) {
      
      #### A) Diff. Analysis Methods ======================
      if (tmp.analysis.Method == 'Complete_Case_Analysis') {
        ##### CCA -------------------------|
        results <- fn.Complete.Case.Analysis(dgp)
        # --------------------------------!
      } else if (tmp.analysis.Method == 'Joint_MVN_MI_Amelia') {
        ##### Joint-Amelia ----------------|
        results <- fn.Joint.MVN.MI(dgp, m = tmp.no.Imputation, command = 'amelia')
        # --------------------------------!
      } else if (tmp.analysis.Method == 'Joint_MVN_MI_Norm') {
        ##### Joint-Norm ------------------|
        results <- fn.Joint.MVN.MI(dgp, m = tmp.no.Imputation, command = 'norm')
        # --------------------------------!
      } else if (tmp.analysis.Method == 'Conditional_MI') {
        ##### Conditional -----------------|
        results <- fn.Conditional.MI(dgp, m = tmp.no.Imputation)
        # --------------------------------!
      } else {
        ##### Else ------------------------|
        stop('analysis.Method is not recognized.')
      }
      
      ### B) Add to Tracking Table =======================
      overall.Results[nrow(overall.Results)+1,] <- c(
        iter
        # Adjustable hyper-parameters
        , tmp.restriction
        , tmp.var.Type
        , tmp.no.Var.Missing
        , tmp.miss.Rate.Y
        , tmp.miss.Rate.X
        # Fied hyper-parameters
        , tmp.N
        , tmp.no.Var.Observed
        , tmp.no.Category
        , tmp.no.Imputation
        # Analysis method
        # , 'Complete_Case_Analysis'
        , tmp.analysis.Method
        , tmp.imp.Method
        # Results
        # -- Acc: Imputed values
        , fn.null.to.NA(results$Imputed.Values.RMSE)
        , fn.null.to.NA(results$Imputed.Values.Bias)
        # -- Acc: Choice probabilities
        , fn.null.to.NA(results$Choice.Prob.RMSE)
        , fn.null.to.NA(results$Choice.Prob.Bias)
        # -- Acc: Coefficients
        , fn.null.to.NA(results$mdl.X1.coef.RMSE)
        , fn.null.to.NA(results$mdl.X1.coef.Bias)
        , fn.null.to.NA(results$mdl.X1.coef.MNS)
        , fn.null.to.NA(results$mdl.Y.coef.RMSE)
        , fn.null.to.NA(results$mdl.Y.coef.Bias)
        , fn.null.to.NA(results$mdl.Y.coef.MNS)
        # -- Acc: Fitted values
        , fn.null.to.NA(results$mdl.X1.fitted.RMSE)
        , fn.null.to.NA(results$mdl.X1.fitted.Bias)
        , fn.null.to.NA(results$mdl.X1.fitted.All.RMSE)
        , fn.null.to.NA(results$mdl.X1.fitted.All.Bias)
        , fn.null.to.NA(results$mdl.Y.fitted.RMSE)
        , fn.null.to.NA(results$mdl.Y.fitted.Bias)
        , fn.null.to.NA(results$mdl.Y.fitted.All.RMSE)
        , fn.null.to.NA(results$mdl.Y.fitted.All.Bias)
        # -- Time to Convergence
        , fn.null.to.NA(results$time)
      )
    }
  }
}







