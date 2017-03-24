library(DNAPL)
context("DNST solver")

cG <- cstG.DNmodel(.5, 20, 1, .5, .1, .1, .2, 1200, 1.2, 20, 3)
cvG <- cnvG.DNmodel(.5, 20, 1, .5, .1, .1, .2, 1200, 1.2, 20, 3)
ddG <- DDpg.DNmodel(.5, 20, 1, .1, .1, .2, 1200, 1.2, 20, 3)
mfdata <- RNetCDF::open.nc(system.file("rflow_mf_demo.nc",
                                       package = "Rflow"))

test_that("works with MODFLOW model, cstG", {
  fnm <- tempfile()
  expect_silent(dnst <- {
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = 0, end.t = 1500, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })
  expect_equal(dnst, readRDS(fnm))
})

test_that("cnvG model", {
  fnm <- tempfile()
  expect_silent(dnst <- {
    DNST(fnm, "test", cvG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = 0, end.t = 1500, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })
  expect_equal(dnst, readRDS(fnm))
})

test_that("DDpg model", {
  fnm <- tempfile()
  expect_silent(dnst <- {
    DNST(fnm, "test", ddG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = 0, end.t = 1500, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })
  expect_equal(dnst, readRDS(fnm))
})

test_that("works with explicit flow", {
  # qh is list of functions
  fnm <- tempfile()
  expect_silent(dnst <- {
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = 0, end.t = 1500, dt = 10,
         qh = rep(list(function(t) .1), cG@NLAY))
  })
  expect_equal(dnst, readRDS(fnm))

  # qh is numeric vector of constant values
  expect_silent(dnst <- {
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = 0, end.t = 1500, dt = 10,
         qh = rep(.1, cG@NLAY))
  })
  expect_equal(dnst, readRDS(fnm))
})

test_that("catches timing and spillage errors", {
  mftrg <- range(Rflow::mftstime(mfdata, TRUE))
  mftmin <- mftrg[1L]; mftmax <- mftrg[2L]

  # DNAPL completely too early - error
  expect_error({
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = mftmin - 1000, end.t = mftmin, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })

  # DNAPL starts too early - warning
  expect_warning({
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = mftmin - 1000, end.t = mftmin + 1000, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })

  # DNAPL completely too late - error
  expect_error({
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = mftmax, end.t = mftmax + 1000, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })

  # DNAPL finishes too late - warning
  expect_warning({
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = mftmax - 1000, end.t = mftmax + 1000, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })

  # DNAPL starts too early and finishes too late - warning
  expect_warning({
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = mftmin - 1000, end.t = mftmax + 1000, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })

  # exactly matches MODFLOw time range - no problem
  expect_silent({
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         start.t = mftmin, end.t = mftmax, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })

  # spillage mistimed - warning
  expect_warning({
    DNST(fnm, "test", cG,
         data.frame(year = c(1800, 1805), cons = c(0, 10)),
         start.t = mftmin, end.t = mftmax, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })

  # negative spillage - error
  expect_error({
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, -10)),
         start.t = mftmin, end.t = mftmax, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })
  expect_error({
    DNST(fnm, "test", cG,
         data.frame(year = c(1900, 1905), cons = c(0, 10)),
         fsphist = data.frame(year = c(1900, 1905), f = c(0, -1)),
         start.t = mftmin, end.t = mftmax, dt = 10, x = 625, y = 825,
         mfdata = mfdata)
  })
})
