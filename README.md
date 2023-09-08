# stat536E_finalProject
Paper Review - Multiple Imputation for Continuous and Categorical Data - Comparing Joint Multivariate Normal and Conditional Approaches

# Stat 536E - Final Project

## A Study on multiple imputation techniques based on the *Multiple Imputation for Continuous and Ordinal Data: Comparing Joint Multivariate Normal and Conditional Approaches* paper by Jonathan Kropko, Ben Goodrich, Andrew Gelman and Jennifer Hill

Date: 04/20/2023

Abstract:

This report studies the performance of two popular approaches to multiple imputation (MI), namely joint multivariate normal (MVN) MI, in which the data is modeled as a jointly MVN distributed, and conditional MI, in which each variable is modeled conditionally on all of the other variables. As a consequence, implementing joint MVN MI requires an extra step of transforming discrete variables from continuous variables, which is often done via probabilistic models. In order to compare the relative performance, we simulate data and score each approach with two types of metrics: (1) the accuracy of the imputed values and (2) the accuracy of the coefcients and ftted values based on a model ftted to their completed datasets.

Building on top of the work of Kropko, Goodrich, Gelman, and Hill in [Kro+17], we extend the simulations on two dimensions: (1) the missingness mechanism to include missing not at random (MNAR) in additiona to missing at random (MAR) and (2) the miss rates in both the variable of interest and the explanatory variables. In consideration as continuous and ordinal (ordered-discrete) variables. We fnd that conditional MI tends to outperform in all metrics when the variable is continuous, whereas there does not appear a clearly dominating approach when the variable is ordinal.

