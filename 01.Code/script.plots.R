
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

load('01.Code/overall.Results.1000.RData')



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



fn.plot.Imputed.Values.restriction <- function(
    tmp.imp.Method, tmp.var.Type,
    tmp.evaluation.Metric, label.evaluation.Metric,
    tmp.metric.Type = 'RMSE'
) {
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
      value = mean(!!sym(paste0(tmp.evaluation.Metric, '.', tmp.metric.Type)))
    )
  # Plot
  p <- tmp %>% 
    ggplot() +
    geom_bar(aes(x = analysis.Method, y = value, fill = analysis.Method),
             width = 0.75, position = 'dodge', stat = 'identity') +
    xlab('Missing Mechanism') +
    ylab(paste0(label.evaluation.Metric, ' - ', tmp.metric.Type)) +
    labs(title = paste0(label.evaluation.Metric, ' - ', tmp.metric.Type)) +
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
  # Return plot
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



fn.plot.Imputed.Values.miss.Rate.Y <- function(
    tmp.imp.Method, tmp.var.Type,
    tmp.evaluation.Metric, label.evaluation.Metric,
    tmp.metric.Type = 'RMSE'
) {
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
      value = mean(!!sym(paste0(tmp.evaluation.Metric, '.', tmp.metric.Type)))
    )
  # Plot
  p <- tmp %>% 
    ggplot() +
    geom_bar(aes(x = analysis.Method, y = value, fill = analysis.Method),
             width = 0.75, position = 'dodge', stat = 'identity') +
    xlab('Miss Rate in Y') +
    ylab(paste0(label.evaluation.Metric, ' - ', tmp.metric.Type)) +
    labs(title = paste0(label.evaluation.Metric, ' - ', tmp.metric.Type)) +
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
    facet_wrap(~miss.Rate.Y) +
    theme(panel.spacing.x = unit(2, "lines"))
  # Return plot
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
Options.value.Type.Continuous <- c(
  'Imputed Values', 'GLM Coefficients', 'GLM Fitted Values'
)
Options.value.Type.Ordinal <- c(
  'Choice Probabilities', 'GLM Coefficients', 'GLM Fitted Values'
)



# ========================================================!
## MAR vs. NMAR -------------------------------------------

dimension <- 'MAR.vs.NMAR'



### RMSE --------------------------------------------------
metric.Type <- 'RMSE'



# ==============================================|
#### PPD:Continuous -----------------------------

var.Type <- 'continuous'
imp.Method <- 'ppd'

p.list <- list()
counter <- 0
for (value.Type in Options.value.Type.Continuous) {
  counter <- counter + 1
  # Prepare Metric
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
  # Plot
  p <- fn.plot.Imputed.Values.restriction(
    tmp.imp.Method          = imp.Method,
    tmp.var.Type            = var.Type,
    tmp.evaluation.Metric   = paste0(evaluation.Metric),
    label.evaluation.Metric = paste0(value.Type),
    tmp.metric.Type = metric.Type
  )
  p.list[[counter]] <- p
  counter <- counter+1
  p.list[[counter]] <- NULL
}

p <- ggarrange(plotlist = p.list, ncol = 1, heights = c(1, 0.1, 1, .1, 1), 
               common.legend = TRUE, legend = 'bottom')
png(paste0('03.Figures/01.MAR_vs_NMAR/', dimension, '.', var.Type, '.', 
           imp.Method, '.', metric.Type, '.png'), 
    width = 500, height = 1000, pointsize = 15)
print(p)
dev.off()




# ==============================================|
#### PPD:Ordinal --------------------------------

var.Type <- 'ordinal'
imp.Method <- 'ppd'

p.list <- list()
counter <- 0
for (value.Type in Options.value.Type.Ordinal) {
  counter <- counter + 1
  # Prepare Metric
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
  # Plot
  p <- fn.plot.Imputed.Values.restriction(
    tmp.imp.Method          = imp.Method,
    tmp.var.Type            = var.Type,
    tmp.evaluation.Metric   = paste0(evaluation.Metric),
    label.evaluation.Metric = paste0(value.Type),
    tmp.metric.Type = metric.Type
  )
  p.list[[counter]] <- p
  counter <- counter+1
  p.list[[counter]] <- NULL
}

p <- ggarrange(plotlist = p.list, ncol = 1, heights = c(1, 0.1, 1, .1, 1), 
               common.legend = TRUE, legend = 'bottom')
png(paste0('03.Figures/01.MAR_vs_NMAR/', dimension, '.', var.Type, '.', 
           imp.Method, '.', metric.Type, '.png'), 
    width = 500, height = 1000, pointsize = 15)
print(p)
dev.off()





# ==============================================|
#### PMM:Continuous -----------------------------

var.Type <- 'continuous'
imp.Method <- 'pmm'

p.list <- list()
counter <- 0
for (value.Type in Options.value.Type.Continuous) {
  counter <- counter + 1
  # Prepare Metric
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
  # Plot
  p <- fn.plot.Imputed.Values.restriction(
    tmp.imp.Method          = imp.Method,
    tmp.var.Type            = var.Type,
    tmp.evaluation.Metric   = paste0(evaluation.Metric),
    label.evaluation.Metric = paste0(value.Type),
    tmp.metric.Type = metric.Type
  )
  p.list[[counter]] <- p
  counter <- counter+1
  p.list[[counter]] <- NULL
}

p <- ggarrange(plotlist = p.list, ncol = 1, heights = c(1, 0.1, 1, .1, 1), 
               common.legend = TRUE, legend = 'bottom')
png(paste0('03.Figures/01.MAR_vs_NMAR/', dimension, '.', var.Type, '.', 
           imp.Method, '.', metric.Type, '.png'), 
    width = 500, height = 1000, pointsize = 15)
print(p)
dev.off()




# ==============================================|
#### PMM:Ordinal --------------------------------

var.Type <- 'ordinal'
imp.Method <- 'pmm'

p.list <- list()
counter <- 0
for (value.Type in Options.value.Type.Ordinal) {
  counter <- counter + 1
  # Prepare Metric
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
  # Plot
  p <- fn.plot.Imputed.Values.restriction(
    tmp.imp.Method          = imp.Method,
    tmp.var.Type            = var.Type,
    tmp.evaluation.Metric   = paste0(evaluation.Metric),
    label.evaluation.Metric = paste0(value.Type),
    tmp.metric.Type = metric.Type
  )
  p.list[[counter]] <- p
  counter <- counter+1
  p.list[[counter]] <- NULL
}

p <- ggarrange(plotlist = p.list, ncol = 1, heights = c(1, 0.1, 1, .1, 1), 
               common.legend = TRUE, legend = 'bottom')
png(paste0('03.Figures/01.MAR_vs_NMAR/', dimension, '.', var.Type, '.', 
           imp.Method, '.', metric.Type, '.png'), 
    width = 500, height = 1000, pointsize = 15)
print(p)
dev.off()






# ========================================================!
## Miss Rate Y --------------------------------------------

dimension <- 'Miss_Rate_Y'

### RMSE --------------------------------------------------
metric.Type <- 'RMSE'



# ==============================================|
#### PPD:Continuous -----------------------------

var.Type <- 'continuous'
imp.Method <- 'ppd'

p.list <- list()
counter <- 0
for (value.Type in Options.value.Type.Continuous) {
  counter <- counter + 1
  # Prepare Metric
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
  # Plot
  p <- fn.plot.Imputed.Values.miss.Rate.Y(
    tmp.imp.Method          = imp.Method,
    tmp.var.Type            = var.Type,
    tmp.evaluation.Metric   = paste0(evaluation.Metric),
    label.evaluation.Metric = paste0(value.Type),
    tmp.metric.Type = metric.Type
  )
  p.list[[counter]] <- p
  counter <- counter+1
  p.list[[counter]] <- NULL
}

p <- ggarrange(plotlist = p.list, ncol = 1, heights = c(1, 0.1, 1, .1, 1), 
               common.legend = TRUE, legend = 'bottom')
png(paste0('03.Figures/02.Miss_Rate_Y/', dimension, '.', var.Type, '.', 
           imp.Method, '.', metric.Type, '.png'), 
    width = 500, height = 1000, pointsize = 15)
print(p)
dev.off()




# ==============================================|
#### PPD:Continuous -----------------------------

var.Type <- 'ordinal'
imp.Method <- 'ppd'

p.list <- list()
counter <- 0
for (value.Type in Options.value.Type.Ordinal) {
  counter <- counter + 1
  # Prepare Metric
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
  # Plot
  p <- fn.plot.Imputed.Values.miss.Rate.Y(
    tmp.imp.Method          = imp.Method,
    tmp.var.Type            = var.Type,
    tmp.evaluation.Metric   = paste0(evaluation.Metric),
    label.evaluation.Metric = paste0(value.Type),
    tmp.metric.Type = metric.Type
  )
  p.list[[counter]] <- p
  counter <- counter+1
  p.list[[counter]] <- NULL
}

p <- ggarrange(plotlist = p.list, ncol = 1, heights = c(1, 0.1, 1, .1, 1), 
               common.legend = TRUE, legend = 'bottom')
png(paste0('03.Figures/02.Miss_Rate_Y/', dimension, '.', var.Type, '.', 
           imp.Method, '.', metric.Type, '.png'), 
    width = 500, height = 1000, pointsize = 15)
print(p)
dev.off()




# ==============================================|
#### PMM:Continuous -----------------------------

var.Type <- 'continuous'
imp.Method <- 'pmm'

p.list <- list()
counter <- 0
for (value.Type in Options.value.Type.Continuous) {
  counter <- counter + 1
  # Prepare Metric
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
  # Plot
  p <- fn.plot.Imputed.Values.miss.Rate.Y(
    tmp.imp.Method          = imp.Method,
    tmp.var.Type            = var.Type,
    tmp.evaluation.Metric   = paste0(evaluation.Metric),
    label.evaluation.Metric = paste0(value.Type),
    tmp.metric.Type = metric.Type
  )
  p.list[[counter]] <- p
  counter <- counter+1
  p.list[[counter]] <- NULL
}

p <- ggarrange(plotlist = p.list, ncol = 1, heights = c(1, 0.1, 1, .1, 1), 
               common.legend = TRUE, legend = 'bottom')
png(paste0('03.Figures/02.Miss_Rate_Y/', dimension, '.', var.Type, '.', 
           imp.Method, '.', metric.Type, '.png'), 
    width = 500, height = 1000, pointsize = 15)
print(p)
dev.off()




# ==============================================|
#### PMM:Continuous -----------------------------

var.Type <- 'ordinal'
imp.Method <- 'pmm'

p.list <- list()
counter <- 0
for (value.Type in Options.value.Type.Ordinal) {
  counter <- counter + 1
  # Prepare Metric
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
  # Plot
  p <- fn.plot.Imputed.Values.miss.Rate.Y(
    tmp.imp.Method          = imp.Method,
    tmp.var.Type            = var.Type,
    tmp.evaluation.Metric   = paste0(evaluation.Metric),
    label.evaluation.Metric = paste0(value.Type),
    tmp.metric.Type = metric.Type
  )
  p.list[[counter]] <- p
  counter <- counter+1
  p.list[[counter]] <- NULL
}

p <- ggarrange(plotlist = p.list, ncol = 1, heights = c(1, 0.1, 1, .1, 1), 
               common.legend = TRUE, legend = 'bottom')
png(paste0('03.Figures/02.Miss_Rate_Y/', dimension, '.', var.Type, '.', 
           imp.Method, '.', metric.Type, '.png'), 
    width = 500, height = 1000, pointsize = 15)
print(p)
dev.off()

