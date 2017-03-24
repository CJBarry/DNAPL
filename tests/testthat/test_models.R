library(DNAPL)
context("DNAPL models")

cG <- cstG.DNmodel(.5, 20, 1, .5, .1, .1, .2, 1200, 1.2, 20, 3)
cvG <- cnvG.DNmodel(.5, 20, 1, .5, .1, .1, .2, 1200, 1.2, 20, 3)
ddG <- DDpg.DNmodel(.5, 20, 1, .1, .1, .2, 1200, 1.2, 20, 3)

test_that("class is S4 DNAPLmodel", {
  expect_s4_class(cG, "DNAPLmodel")
  expect_s4_class(cvG, "DNAPLmodel")
  expect_s4_class(ddG, "DNAPLmodel")
})

test_that("cstG and cnvG model functions are set up correctly", {
  expect_silent(DNAPLmodel.check(cG))
  expect_silent(DNAPLmodel.check(cvG))
  expect_silent(DNAPLmodel.check(ddG))
})
