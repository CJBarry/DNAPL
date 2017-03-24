library(DNAPL)
context("cstG")

cG <- cstG.DNmodel(.5, 20, 1, .5, .1, .1, .2, 1200, 1.2, 20, 3)

test_that("class is S4 DNAPLmodel", {
  expect_s4_class(cG, "DNAPLmodel")
})
