
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
libs <- c('tidyverse', 'mi', 'arm', 'Amelia', 'mvnmle', 'MASS', 'norm2', 'norm', 'nnet')
lapply(libs, require, character.only = TRUE)
rm(libs)



# ============================================================================!
# 1) Helper fns ---------------------------------------------------------------
# ----------------------------------------------------------------------------!

# Generate random data
source('01.Code/fn.DGP.R')
# Complete Case Analysis
source('01.Code/fn.Complete.Case.Analysis.R')
# MI: 1) Conditional MI
source('01.Code/fn.Conditional.MI.R')
# MI: 2) Joint MI
source('01.Code/fn.Joint.MVN.MI.helper.PMM.R')
source('01.Code/fn.Joint.MVN.MI.R')




# ============================================================================!
# 2) Hyper-Parameters ---------------------------------------------------------
# ----------------------------------------------------------------------------!


## Fixed --------------------------------------------------
global.N               <- 1000  #(obs)
global.no.Var.Observed <- 5    #(var. fully observed)
global.no.Category     <- 5    #(categories) per discrete var.
global.strong          <- 0    # applicable only when "triangular"
global.estimate.CPCs   <- F    # F to save time, default is T
global.no.Imputation   <- 5    # m, the number of imputation per 1 dataset


## Adjustable ---------------------------------------------

### A) Methodology --------------------
Options.imp.Method     <- c('ppd', 'pmm')

### B) Scenario -----------------------
Options.restriction     <- c('triangular', 'none')
Options.var.Type        <- c('continuous', 'ordinal')
Options.no.Var.Missing  <- seq(1, 5, 2)
Options.miss.Rate.X     <- seq(.1, .3, .1)
Options.miss.Rate.Y     <- c(.25, .5)
Options.analysis.Method <- c(
  'Complete_Case_Analysis',
  'Joint_MVN_MI_Amelia',
  'Joint_MVN_MI_Norm',
  'Conditional_MI'
)




# ============================================================================!
# 3) Sim. Tracking DF ---------------------------------------------------------
# ----------------------------------------------------------------------------!

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
  , time.to.MI.Convergence = numeric()
  , time.to.MI.Analyze     = numeric()
)

fn.null.to.NA <- function(x) {ifelse(is.null(x), NA, x)}



# ============================================================================!
# 4) Simulations --------------------------------------------------------------
# ----------------------------------------------------------------------------!

counter <- 1
no.Iteration <- 5

## 0: Iteration -------------------------------------------
for (iter in 1:no.Iteration) {
  
  ## 1: Missing Mechanism -----------------------------------
  for (tmp.restriction in Options.restriction) {
    
    ## 2: Var. Type -----------------------------------------
    for (tmp.var.Type in Options.var.Type) {
      
      ## 3: No. Var. Missing --------------------------------
      for (tmp.no.Var.Missing in Options.no.Var.Missing) {
        
        ## 4: Miss Rate Y -----------------------------------
        for (tmp.miss.Rate.Y in Options.miss.Rate.Y) {
          
          ## 5: Miss Rate X ---------------------------------
          for (tmp.miss.Rate.X in Options.miss.Rate.X) {
            
            ## 6: imp.Method --------------------------------
            for (tmp.imp.Method in c('ppd', 'pmm')) {
            
              ## =================================================================!
              ### 91: Simulate Data -----------------------
              dgp <- fn.DGP(
                N               = global.N,
                var.Type        = tmp.var.Type,
                no.Category     = global.no.Category,
                no.Var.Observed = global.no.Var.Observed,
                no.Var.Missing  = tmp.no.Var.Missing,
                miss.Rate.X     = tmp.miss.Rate.X,
                miss.Rate.Y     = tmp.miss.Rate.Y,
                restriction     = tmp.restriction,
                strong          = global.strong,
                imp.Method      = tmp.imp.Method)
              
              ## =================================================================!
              ## 92: Analyze w. diff. Methods -------------
              for (tmp.analysis.Method in Options.analysis.Method) {

                ## =================================================================!
                ### A) Analyze ============================
                time.Start   <- proc.time()
                if (tmp.analysis.Method == 'Complete_Case_Analysis') {
                  results <- fn.Complete.Case.Analysis(dgp)
                  # --------------------------------!
                } else if (tmp.analysis.Method == 'Joint_MVN_MI_Amelia') {
                  results <- fn.Joint.MVN.MI(dgp, m = global.no.Imputation, command = 'amelia')
                  # --------------------------------!
                } else if (tmp.analysis.Method == 'Joint_MVN_MI_Norm') {
                  results <- fn.Joint.MVN.MI(dgp, m = global.no.Imputation, command = 'norm')
                  # --------------------------------!
                } else if (tmp.analysis.Method == 'Conditional_MI') {
                  results <- fn.Conditional.MI(dgp, m = global.no.Imputation)
                  # --------------------------------!
                } else {
                  ##### Else ------------------------|
                  stop('analysis.Method is not recognized.')
                }
                time.End        <- proc.time()
                time.to.Analyze <- as.numeric((time.End - time.Start)[3])

                ## =================================================================!
                ### B) Add to Tracking Table ================
                overall.Results[nrow(overall.Results)+1,] <- c(
                  iter
                  # Adjustable hyper-parameters
                  , tmp.restriction
                  , tmp.var.Type
                  , tmp.no.Var.Missing
                  , tmp.miss.Rate.Y
                  , tmp.miss.Rate.X
                  # Fied hyper-parameters
                  , global.N
                  , global.no.Var.Observed
                  , global.no.Category
                  , global.no.Imputation
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
                  , time.to.Analyze
                )
                
                counter <- counter+1
              }
            }
          }
        }
      }
    }
  }
}







# ============================================================================!
# 5) Fix df str ---------------------------------------------------------------
# ----------------------------------------------------------------------------!

str(overall.Results)
idx <- which(colnames(overall.Results) %in% c('Iteration', 'no.Var.Missing', 'N', 
                                              'no.Var.Observed', 'no.Category','no.Imputation'))
overall.Results[idx] <- sapply(overall.Results[idx], as.integer)

idx <- c(which(colnames(overall.Results) %in% c('miss.Rate.Y', 'miss.Rate.X')),
         13:ncol(overall.Results))
overall.Results[idx] <- sapply(overall.Results[idx], as.numeric)
str(overall.Results)


save(overall.Results, file = 'overall.Results.1000.RData')
write.table(overall.Results, file='overall.Results.1000.csv', sep=",", col.names=TRUE, row.names=FALSE)
