
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
libs <- c('dplyr', 'ggplot2', 'ggthemes', 'ggpubr', 'reshape2')
lapply(libs, require, character.only = TRUE)
rm(libs)





# ============================================================================!
# 1) Load Simulation Data -----------------------------------------------------
# ----------------------------------------------------------------------------!

load('01.Code/overall.Results.RData')



colnames(overall.Results)


# ============================================================================!
# 2) fn -----------------------------------------------------------------------
# ----------------------------------------------------------------------------!


# ========================================================!
## MAR vs. NMAR -------------------------------------------
fn.plot.restriction <- function(tmp.var.Type, 
                                tmp.evaluation.Metric, label.evaluation.Metric,
                                tmp.imp.Method, legend.text.Indicator) {
  # Get data
  tmp <- overall.Results %>% 
    mutate(
      `Missing Mechanism` = ifelse(restriction == 'triangular', 'MAR', 'NMAR')
    ) %>% 
    filter(
      miss.Rate.Y    == .25,
      miss.Rate.X    == .1,
      no.Var.Missing == 3,
    ) %>% 
    filter(
      imp.Method == tmp.imp.Method,
      var.Type   == tmp.var.Type
    ) %>% 
    group_by(
      `Missing Mechanism`,
      analysis.Method
    ) %>% 
    summarize(
      value = mean(!!sym(tmp.evaluation.Metric))
    )
  # Plot
  p <- tmp %>% 
    ggplot() +
    geom_bar(aes(x = analysis.Method, y = value, fill = analysis.Method),
             width = 0.75, position = 'dodge', stat = 'identity') +
    xlab('Missing Mechanism') +
    ylab(label.evaluation.Metric) +
    labs(title = label.evaluation.Metric) +
    guides(color = 'none') + 
    scale_fill_brewer(palette = "Set2",
                      direction = -1,
                      name = 'Missing Mechanism',
                      labels = c('Complete Case Analysis', 'Conditional MI',
                                 'Joint MVN MI - Amelia', 'Joint MVN MI - norm')) +
    theme(plot.background    = element_rect(fill = 'white'),
          panel.background   = element_rect(fill = 'white'),
          panel.grid.major   = element_line(linetype = 'dashed',
                                            size = .5, color = 'lightgray'),
          axis.text.x = element_blank(),
          # Legend info
          legend.background = element_rect(fill = "white"),
          legend.position   = 'bottom',
          legend.title      = element_blank()
    ) +
    facet_wrap(~`Missing Mechanism`) +
    theme(panel.spacing.x = unit(2, "lines"))
  return(p)
}







# ========================================================!
## Miss Rate Y --------------------------------------------
fn.plot.miss.Rate.Y <- function(tmp.var.Type, 
                                tmp.evaluation.Metric, label.evaluation.Metric,
                                tmp.imp.Method, legend.text.Indicator) {
  # Get data
  tmp <- overall.Results %>% 
    # mutate(
    #   `Missing Mechanism` = ifelse(restriction == 'triangular', 'MAR', 'NMAR')
    # ) %>% 
    filter(
      restriction    == 'triangular',
      # miss.Rate.Y    == .25,
      miss.Rate.X    == .1,
      no.Var.Missing == 3,
    ) %>% 
    filter(
      imp.Method == tmp.imp.Method,
      var.Type   == tmp.var.Type
    ) %>% 
    group_by(
      miss.Rate.Y,
      analysis.Method
    ) %>% 
    summarize(
      value = mean(!!sym(tmp.evaluation.Metric))
    )
  # Plot
  p <- tmp %>% 
    ggplot() +
    geom_bar(aes(x = analysis.Method, y = value, fill = analysis.Method),
             width = 0.75, position = 'dodge', stat = 'identity') +
    xlab('Missing Mechanism') +
    ylab(label.evaluation.Metric) +
    labs(title = label.evaluation.Metric) +
    guides(color = 'none') + 
    scale_fill_brewer(palette = "Set2",
                      direction = -1,
                      name = 'Miss Rate in Y, the variable of interest',
                      labels = c('Complete Case Analysis', 'Conditional MI',
                                 'Joint MVN MI - Amelia', 'Joint MVN MI - norm')) +
    theme(plot.background    = element_rect(fill = 'white'),
          panel.background   = element_rect(fill = 'white'),
          panel.grid.major   = element_line(linetype = 'dashed',
                                            size = .5, color = 'lightgray'),
          axis.text.x = element_blank(),
          # Legend info
          legend.background = element_rect(fill = "white"),
          legend.position   = 'bottom',
          legend.title      = element_blank()
    ) +
    facet_wrap(~ miss.Rate.Y) +
    theme(panel.spacing.x = unit(2, "lines"))
  return(p)
}





colnames(overall.Results)

# ============================================================================!
# 3) Plots --------------------------------------------------------------------
# ----------------------------------------------------------------------------!

Options.value.Type <- c(
  'Imputed Values', 'Choice Probabilities',
  'GLM Coefficients', 'GLM Fitted Values'
)



# ========================================================!
## MAR vs. NMAR -------------------------------------------
dimension <- 'MAR.vs.NMAR'

imp.Method <- 'ppd'

### PPD-Continuous --------------------
var.Type   <- 'continuous'

for (value.Type in Options.value.Type) {
  
  if (value.Type == 'Imputed Values') {
    evaluation.Metric <- 'Imputed.Values'
  } else if (value.Type == 'Choice Probabilities') {
    evaluation.Metric <- 'Choice.Prob'
  } else if (value.Type == 'GLM Coefficients') {
    evaluation.Metric <- 'mdl.Y.coef'
  } else if (value.Type == 'GLM Fitted Values') {
    evaluation.Metric <- 'mdl.Y.fitted.All'
  } else {
    stop('value.Type not recognized.')
  }
  
  p_RMSE <- fn.plot.restriction(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.RMSE'),
    label.evaluation.Metric = paste0(value.Type, ' - RMSE'),
    legend.text.Indicator   = 'Yes'
  )
  p_Bias <- fn.plot.restriction(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.Bias'),
    label.evaluation.Metric = paste0(value.Type, ' - Bias'),
    legend.text.Indicator   = 'Yes'
  )  
  p <- ggarrange(plotlist = list(p_RMSE, p_Bias), 
                 ncol = 1, common.legend = TRUE , legend = "bottom")
  png(paste0('03.Figures/01.MAR_vs_NMAR/', dimension, '.', var.Type, '.', imp.Method, 
             '.', gsub(' ', '\\.', value.Type), '.png'), 
      width = 800, height = 400, pointsize = 15)
  print(p)
  dev.off()
}


### PPD-Ordinal -----------------------
var.Type   <- 'ordinal'

for (value.Type in Options.value.Type) {
  
  if (value.Type == 'Imputed Values') {
    evaluation.Metric <- 'Imputed.Values'
  } else if (value.Type == 'Choice Probabilities') {
    evaluation.Metric <- 'Choice.Prob'
  } else if (value.Type == 'GLM Coefficients') {
    evaluation.Metric <- 'mdl.Y.coef'
  } else if (value.Type == 'GLM Fitted Values') {
    evaluation.Metric <- 'mdl.Y.fitted.All'
  } else {
    stop('value.Type not recognized.')
  }
  
  p_RMSE <- fn.plot.restriction(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.RMSE'),
    label.evaluation.Metric = paste0(value.Type, ' - RMSE'),
    legend.text.Indicator   = 'Yes'
  )
  p_Bias <- fn.plot.restriction(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.Bias'),
    label.evaluation.Metric = paste0(value.Type, ' - Bias'),
    legend.text.Indicator   = 'Yes'
  )  
  p <- ggarrange(plotlist = list(p_RMSE, p_Bias), 
                 ncol = 1, common.legend = TRUE , legend = "bottom")
  png(paste0('03.Figures/01.MAR_vs_NMAR/', dimension, '.', var.Type, '.', imp.Method, 
             '.', gsub(' ', '\\.', value.Type), '.png'), 
      width = 800, height = 400, pointsize = 15)
  print(p)
  dev.off()
}




imp.Method <- 'pmm'

### PMM-Continuous --------------------
var.Type   <- 'continuous'

for (value.Type in Options.value.Type) {
  
  if (value.Type == 'Imputed Values') {
    evaluation.Metric <- 'Imputed.Values'
  } else if (value.Type == 'Choice Probabilities') {
    evaluation.Metric <- 'Choice.Prob'
  } else if (value.Type == 'GLM Coefficients') {
    evaluation.Metric <- 'mdl.Y.coef'
  } else if (value.Type == 'GLM Fitted Values') {
    evaluation.Metric <- 'mdl.Y.fitted.All'
  } else {
    stop('value.Type not recognized.')
  }
  
  p_RMSE <- fn.plot.restriction(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.RMSE'),
    label.evaluation.Metric = paste0(value.Type, ' - RMSE'),
    legend.text.Indicator   = 'Yes'
  )
  p_Bias <- fn.plot.restriction(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.Bias'),
    label.evaluation.Metric = paste0(value.Type, ' - Bias'),
    legend.text.Indicator   = 'Yes'
  )  
  p <- ggarrange(plotlist = list(p_RMSE, p_Bias), 
                 ncol = 1, common.legend = TRUE , legend = "bottom")
  png(paste0('03.Figures/01.MAR_vs_NMAR/', dimension, '.', var.Type, '.', imp.Method, 
             '.', gsub(' ', '\\.', value.Type), '.png'), 
      width = 800, height = 400, pointsize = 15)
  print(p)
  dev.off()
}


### PMM-Ordinal -----------------------
var.Type   <- 'ordinal'

for (value.Type in Options.value.Type) {
  
  if (value.Type == 'Imputed Values') {
    evaluation.Metric <- 'Imputed.Values'
  } else if (value.Type == 'Choice Probabilities') {
    evaluation.Metric <- 'Choice.Prob'
  } else if (value.Type == 'GLM Coefficients') {
    evaluation.Metric <- 'mdl.Y.coef'
  } else if (value.Type == 'GLM Fitted Values') {
    evaluation.Metric <- 'mdl.Y.fitted.All'
  } else {
    stop('value.Type not recognized.')
  }
  
  p_RMSE <- fn.plot.restriction(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.RMSE'),
    label.evaluation.Metric = paste0(value.Type, ' - RMSE'),
    legend.text.Indicator   = 'Yes'
  )
  p_Bias <- fn.plot.restriction(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.Bias'),
    label.evaluation.Metric = paste0(value.Type, ' - Bias'),
    legend.text.Indicator   = 'Yes'
  )  
  p <- ggarrange(plotlist = list(p_RMSE, p_Bias), 
                 ncol = 1, common.legend = TRUE , legend = "bottom")
  png(paste0('03.Figures/01.MAR_vs_NMAR/', dimension, '.', var.Type, '.', imp.Method, 
             '.', gsub(' ', '\\.', value.Type), '.png'), 
      width = 800, height = 400, pointsize = 15)
  print(p)
  dev.off()
}




# ========================================================!
## Miss Rate Y --------------------------------------------
dimension <- 'Miss.Rate.Y'

imp.Method <- 'ppd'

### PPD-Continuous --------------------
var.Type   <- 'continuous'

for (value.Type in Options.value.Type) {
  
  if (value.Type == 'Imputed Values') {
    evaluation.Metric <- 'Imputed.Values'
  } else if (value.Type == 'Choice Probabilities') {
    evaluation.Metric <- 'Choice.Prob'
  } else if (value.Type == 'GLM Coefficients') {
    evaluation.Metric <- 'mdl.Y.coef'
  } else if (value.Type == 'GLM Fitted Values') {
    evaluation.Metric <- 'mdl.Y.fitted.All'
  } else {
    stop('value.Type not recognized.')
  }
  
  p_RMSE <- fn.plot.miss.Rate.Y(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.RMSE'),
    label.evaluation.Metric = paste0(value.Type, ' - RMSE'),
    legend.text.Indicator   = 'Yes'
  )
  p_Bias <- fn.plot.miss.Rate.Y(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.Bias'),
    label.evaluation.Metric = paste0(value.Type, ' - Bias'),
    legend.text.Indicator   = 'Yes'
  )  
  p <- ggarrange(plotlist = list(p_RMSE, p_Bias), 
                 ncol = 1, common.legend = TRUE , legend = "bottom")
  png(paste0('03.Figures/02.Miss_Rate_Y/', dimension, '.', var.Type, '.', imp.Method, 
             '.', gsub(' ', '\\.', value.Type), '.png'), 
      width = 800, height = 400, pointsize = 15)
  print(p)
  dev.off()
}


### PPD-Ordinal -----------------------
var.Type   <- 'ordinal'

for (value.Type in Options.value.Type) {
  
  if (value.Type == 'Imputed Values') {
    evaluation.Metric <- 'Imputed.Values'
  } else if (value.Type == 'Choice Probabilities') {
    evaluation.Metric <- 'Choice.Prob'
  } else if (value.Type == 'GLM Coefficients') {
    evaluation.Metric <- 'mdl.Y.coef'
  } else if (value.Type == 'GLM Fitted Values') {
    evaluation.Metric <- 'mdl.Y.fitted.All'
  } else {
    stop('value.Type not recognized.')
  }
  
  p_RMSE <- fn.plot.miss.Rate.Y(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.RMSE'),
    label.evaluation.Metric = paste0(value.Type, ' - RMSE'),
    legend.text.Indicator   = 'Yes'
  )
  p_Bias <- fn.plot.miss.Rate.Y(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.Bias'),
    label.evaluation.Metric = paste0(value.Type, ' - Bias'),
    legend.text.Indicator   = 'Yes'
  )  
  p <- ggarrange(plotlist = list(p_RMSE, p_Bias), 
                 ncol = 1, common.legend = TRUE , legend = "bottom")
  png(paste0('03.Figures/02.Miss_Rate_Y/', dimension, '.', var.Type, '.', imp.Method, 
             '.', gsub(' ', '\\.', value.Type), '.png'), 
      width = 800, height = 400, pointsize = 15)
  print(p)
  dev.off()
}




imp.Method <- 'pmm'

### PMM-Continuous --------------------
var.Type   <- 'continuous'

for (value.Type in Options.value.Type) {
  
  if (value.Type == 'Imputed Values') {
    evaluation.Metric <- 'Imputed.Values'
  } else if (value.Type == 'Choice Probabilities') {
    evaluation.Metric <- 'Choice.Prob'
  } else if (value.Type == 'GLM Coefficients') {
    evaluation.Metric <- 'mdl.Y.coef'
  } else if (value.Type == 'GLM Fitted Values') {
    evaluation.Metric <- 'mdl.Y.fitted.All'
  } else {
    stop('value.Type not recognized.')
  }
  
  p_RMSE <- fn.plot.miss.Rate.Y(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.RMSE'),
    label.evaluation.Metric = paste0(value.Type, ' - RMSE'),
    legend.text.Indicator   = 'Yes'
  )
  p_Bias <- fn.plot.miss.Rate.Y(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.Bias'),
    label.evaluation.Metric = paste0(value.Type, ' - Bias'),
    legend.text.Indicator   = 'Yes'
  )  
  p <- ggarrange(plotlist = list(p_RMSE, p_Bias), 
                 ncol = 1, common.legend = TRUE , legend = "bottom")
  png(paste0('03.Figures/02.Miss_Rate_Y/', dimension, '.', var.Type, '.', imp.Method, 
             '.', gsub(' ', '\\.', value.Type), '.png'), 
      width = 800, height = 400, pointsize = 15)
  print(p)
  dev.off()
}


### PMM-Ordinal -----------------------
var.Type   <- 'ordinal'

for (value.Type in Options.value.Type) {
  
  if (value.Type == 'Imputed Values') {
    evaluation.Metric <- 'Imputed.Values'
  } else if (value.Type == 'Choice Probabilities') {
    evaluation.Metric <- 'Choice.Prob'
  } else if (value.Type == 'GLM Coefficients') {
    evaluation.Metric <- 'mdl.Y.coef'
  } else if (value.Type == 'GLM Fitted Values') {
    evaluation.Metric <- 'mdl.Y.fitted.All'
  } else {
    stop('value.Type not recognized.')
  }
  
  p_RMSE <- fn.plot.miss.Rate.Y(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.RMSE'),
    label.evaluation.Metric = paste0(value.Type, ' - RMSE'),
    legend.text.Indicator   = 'Yes'
  )
  p_Bias <- fn.plot.miss.Rate.Y(
    tmp.var.Type            = var.Type,
    tmp.imp.Method          = imp.Method,
    tmp.evaluation.Metric   = paste0(evaluation.Metric, '.Bias'),
    label.evaluation.Metric = paste0(value.Type, ' - Bias'),
    legend.text.Indicator   = 'Yes'
  )  
  p <- ggarrange(plotlist = list(p_RMSE, p_Bias), 
                 ncol = 1, common.legend = TRUE , legend = "bottom")
  png(paste0('03.Figures/02.Miss_Rate_Y/', dimension, '.', var.Type, '.', imp.Method, 
             '.', gsub(' ', '\\.', value.Type), '.png'), 
      width = 800, height = 400, pointsize = 15)
  print(p)
  dev.off()
}

