library(DNAPL)
context("datasets and utilities")

test_that("check usage data sets", {
  e <- environment()
  expect_silent(data("TCEusage", envir = e))
  expect_silent(data("PCEusage", envir = e))
  expect_silent(data("TCAusage", envir = e))
  expect_silent(data("TeCMusage", envir = e))

  expect_named(get("TCEusage", e), c("year", "cons"))
  expect_named(get("PCEusage", e), c("year", "cons"))
  expect_named(get("TCAusage", e), c("year", "cons"))
  expect_named(get("TeCMusage", e), c("year", "cons"))
})

test_that("check UK.to.site", {
  e <- environment()
  data("TCEusage", envir = e)
  tce <- get("TCEusage", e)
  pu <- sample(seq(9, 11, .1), 1L)

  # with loaded object
  expect_silent(su <- UK.to.site(tce, pu))
  expect_equal(max(su$cons), pu)

  # with name of data set
  expect_silent(su <- UK.to.site("TCEusage", pu))
  expect_equal(max(su$cons), pu)
})
