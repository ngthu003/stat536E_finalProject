
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
#' tmp.N               <- 40           # #(obs)
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
#' tmp.imp.Method      <- 'ppd'
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
  imp.Method      = tmp.imp.Method,
  estimate.CPCs   = tmp.estimate.CPCs)



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
  imp.Method      = tmp.imp.Method,
  estimate.CPCs   = tmp.estimate.CPCs)

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

m <- 3
cond.MI.cts <- fn.Conditional.MI(dgp.cts, m = m)
data.frame(t(sapply(cond.MI.cts,c)))
cond.MI.ord <- fn.Conditional.MI(dgp.ord, m = m, no.Category = dgp.ord$dgp.Settings$no.Category)
data.frame(t(sapply(cond.MI.ord,c)))





# ============================================================================!
# 3) Joint MVN MI Analysis ----------------------------------------------------
# ----------------------------------------------------------------------------!

source('01.Code/fn.Joint.MVN.MI.R')

m <- 3
joint.MI.amelia.cts <- fn.Joint.MVN.MI(dgp.cts, m = m, command = 'amelia')
data.frame(t(sapply(joint.MI.amelia.cts,c)))
joint.MI.amelia.ord <- fn.Joint.MVN.MI(dgp.ord, m = m, command = 'amelia')
data.frame(t(sapply(joint.MI.amelia.ord,c)))
