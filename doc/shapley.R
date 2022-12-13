## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(shapleyAIC)

## -----------------------------------------------------------------------------
Y <- swiss[,1]
X <- as.matrix(swiss[,2:6])

shaps <- shapley(Y,X)

shaps$AIC_shapleys
shaps$rsq_shapleys

## ----error=TRUE---------------------------------------------------------------
Y <- swiss[,1]
X <- swiss[,2:6]

shaps <- shapley(Y,X)

