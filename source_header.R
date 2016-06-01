library(genlasso)
library(assertthat)
library(testthat)
library(quadprog)
library(doMC)
library(foreach)
library(glmnet)

source("source_main.R") #main fused lasso stuff
source("source_plotter.R") #main plotting
source("source_threshold.R") #seeing if we can threshold post-fused lasso
source("source_filter.R")
source("source_filter_localpermute.R")
source("source_filter_staircase.R")
source("source_puffer.R")
source("source_common.R")
source("source_plotter_haus.R")
source("source_intersim.R")
source("source_graph.R")
source("source_plotter_graph.R")

set.seed(10)
DATE = Sys.Date()
