library(DNAPL)
context("DNAPL models")

cG <- cstG.DNmodel(.5, 20, 1, .5, .1, .1, .2,
                   CHCprops["TCE", "density"],
                   CHCprops["TCE", "solubility"], 20, 3)
cvG <- cnvG.DNmodel(.5, 20, 1, .5, .1, .1, .2,
                    CHCprops["TCE", "density"],
                    CHCprops["TCE", "solubility"], 20, 3)
ddG1 <- DDpg.DNmodel(.5, 20, 1, .1, .1, .2,
                     CHCprops["TCE", "density"],
                     CHCprops["TCE", "solubility"], 20, 3)
ddG2 <- DDpg.DNmodel(.5, 20, 1, .1, .1, .2,
                     CHCprops["TCE", "density"],
                     CHCprops["TCE", "solubility"], 20, 3, cnvG = TRUE)

test_that("class is S4 DNAPLmodel", {
  expect_s4_class(cG, "DNAPLmodel")
  expect_s4_class(cvG, "DNAPLmodel")
  expect_s4_class(ddG1, "DNAPLmodel")
  expect_s4_class(ddG2, "DNAPLmodel")
})

test_that("cstG and cnvG model functions are set up correctly", {
  expect_silent(DNAPLmodel.check(cG))
  expect_silent(DNAPLmodel.check(cvG))
  expect_silent(DNAPLmodel.check(ddG1))
  expect_silent(DNAPLmodel.check(ddG2))
})
