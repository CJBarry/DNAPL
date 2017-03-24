library(DNAPL)
context("DNST solver")

cG <- cstG.DNmodel(.5, 20, 1, .5, .1, .1, .2, 1200, 1.2, 20, 3)
mfdata <- RNetCDF::open.nc(system.file("rflow_mf_demo.nc",
                                       package = "Rflow"))

test_that("works with MODFLOW model", {
  fnm <- tempfile()
  expect_silent(dnst <- {
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = 0, end.t = 1500, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })
  expect_equal(dnst, readRDS(fnm))
})

test_that("works with explicit flow", {
  fnm <- tempfile()
  expect_silent(dnst <- {
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = 0, end.t = 1500, dt = 10,
         qh = rep(list(function(t) .1), cG@NLAY))
  })
  expect_equal(dnst, readRDS(fnm))

  expect_silent(dnst <- {
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = 0, end.t = 1500, dt = 10,
         qh = rep(.1, cG@NLAY))
  })
  expect_equal(dnst, readRDS(fnm))
})
