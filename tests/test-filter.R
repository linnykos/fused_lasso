library(assertthat)
library(testthat)
source("~/ryan/fused.git/source_filter.R")

test_that("Test that the max of shuffling is right", {
  y = 1:1000
  max.dist = 10
  trials = 100

  for(i in 1:trials){
    set.seed(10*i)

    y2 = .shuffle.values(y, max.dist)
    assert_that(max(abs(y2-y)) <= max.dist)
  }
})
