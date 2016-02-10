library(genlasso)
library(assertthat)
library(testthat)
library(quadprog)
library(doMC)
library(foreach)

source("source_main.R") #main fused lasso stuff
source("source_plotter.R") #main plotting
source("source_recursive.R") #INCOMPLETE, recursively remove changepoints
source("source_threshold.R") #seeing if we can threshold post-fused lasso

set.seed(10)
DATE = Sys.Date()
