# DNAPL package - a convenient master function

#' DNAPL source term master
#'
#' A convenient master function for DNAPL source terms when a MODFLOW
#'  data set and the contaminant data frames in the package are being used.
#'
#' @param contnt
#' \code{"TCE"}, \code{"PCE"}, \code{"TCA"} or {"TeCM"}
#' @param DNmodel
#' function (as object or character string);
#' function to create object of class \link{DNAPLmodel}, such as
#'  \code{\link{cstG.DNmodel}}
#' @param Dm.pars
#' named or correctly ordered list;
#'  arguments for \code{DNmodel} (see \code{\link[base]{formals}}),
#'  excluding hL, NLAY, Cs and rho; \code{\link[base]{do.call}} is used
#' @inheritParams DNST
#' @inheritParams hL.setup
#' @inheritParams UK.to.site
#'
#' @return
#' \link{DNAPLSourceTerm} object, which is also written to
#'  \code{result.file}.  The results are returned invisibly as they will be
#'  too complex to be usefully printed to screen.  Assign the results to an
#'  object or read from \code{result.file} using
#'  \code{\link[base]{readRDS}}.
#'
#' @export
#'
DNST_MASTER <- function(result.file, description,
                        x, y = NULL, mfdata, contnt,
                        DNmodel, Dm.pars, pu,
                        start.t = td::td(1, 1, 1925), end.t, dt = 20,
                        fw = 1, fi = 1,
                        fsphist = data.frame(year = c(1974, 1990),
                                             f = c(1, 0)),
                        ndlpmfl = 1L, bph = 0, z0 = "base"){
  # set up layers
  hL <- hL.setup(mfdata, x, y, ndlpmfl, bph, z0)

  # read contaminant-related data
  e <- environment()
  data(list = paste0(contnt, "usage"), envir = e)
  nat.uhist <- get(paste0(contnt, "usage"), e)
  uhist <- UK.to.site(nat.uhist)
  data(list = "CHCprops", envir = e)
  Cs <- CHCprops[contnt, "solubility"]
  rho <- CHCprops[contnt, "density"]

  # set up DNAPL model
  DNmodel <- match.fun(DNmodel)
  Dm.pars <- c(Dm.pars, list(Cs = Cs, rho = rho, hL = hL))
  DNAPLmodel <- do.call(DNmodel, Dm.pars)

  # run solver
  DNST(result.file = result.file, description = description,
       DNAPLmodel = DNAPLmodel, uhist = uhist, fsphist = fsphist,
       fw = fw, fi = fi, start.t = start.t, end.t = end.t, dt = dt,
       x = x, y = y, z0 = z0, mfdata = mfdata)
}
