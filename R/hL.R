# DNAPL package - setting up model layers based on a MODFLOW model

#' Make DNAPL model layer heights appropriate to a MODFLOW data set
#'
#' @param mfdata
#' NetCDF; a MODFLOW data set (see \code{\link[Rflow]{GW.nc}})
#' @param x,y
#' for \code{\link[grDevices]{xy.coords}}; location of the source zone, in
#'  the same co-ordinate system as \code{mfdata}; one location only
#' @param ndlpmfl
#' integer \code{[1]} or \code{[mfNLAY]}, where \code{mfNLAY} is the number
#'  of layers in the MODFLOW model;
#' number of DNAPL model layers per MODFLOW layer (including any which are
#'  always dry at the location, although the values for these layers won't
#'  make any difference)
#' @param z0
#' numeric \code{[1]} or \code{"base"};
#' elevation of the bottom of the DNAPL source zone, with the same datum as
#'  \code{mfdata}; \code{"base"} instructs to read the bottom of the
#'  MODFLOW model at this location
#'
#' @return
#' numeric;
#' heights of layers of the DNAPL model, from top to bottom
#'
#' @import Rflow
#' @importFrom stats weighted.mean
#' @importFrom abind adrop
#' @importFrom RNetCDF var.get.nc
#' @importFrom RNetCDF dim.inq.nc
#' @importFrom grDevices xy.coords
#' @export
#'
#' @examples
#' mfdata <- RNetCDF::open.nc(system.file("rflow_mf_demo.nc",
#'                                        package = "Rflow"))
#'
#' hL.setup(mfdata, 625, 825, 1L)
#' hL.setup(mfdata, 625, 825, 3L)
#' hL.setup(mfdata, 625, 825, 3L, 10)
#'
hL.setup <- function(mfdata, x, y = NULL, ndlpmfl = 1L, z0 = "base"){
  nmfl <- dim.inq.nc(mfdata, "NLAY")$length
  ndlpmfl <- expand.vec(ndlpmfl, nmfl)
  if(length(ndlpmfl) != nmfl) stop({
    "hL.setup: incorrect length for ndlpmfl"
  })

  xy <- xy.coords(x, y)
  x <- xy$x; y <- xy$y
  stopifnot(length(x) == 1L); stopifnot(length(y) == 1L)
  mfC <- cellref.loc(x, gccs(mfdata, TRUE))
  mfR <- cellref.loc(y, grcs(mfdata, TRUE), TRUE)

  mfldivs <- var.get.nc(mfdata, "elev",
                        c(mfC, mfR, NA),
                        c(1L, 1L, NA))
  if(z0 == "base" || z0 < mfldivs[nmfl + 1L]) z0 <- mfldivs[nmfl + 1L]

  # find the average water table height
  hds <- adrop(var.get.nc(mfdata, "Head",
                          c(mfC, mfR, NA, NA), c(1L, 1L, NA, NA),
                          collapse = FALSE), c(TRUE, TRUE, FALSE, FALSE))
  wtab <- apply(hds, 2L, function(lhd) lhd[which(!is.na(lhd))[1L]])
  tslngs <- diff(mftstime(mfdata))
  avwt <- weighted.mean(wtab, tslngs)

  # find which layer contains the average water table
  ltop <- cellref.loc(avwt, rev(mfldivs), TRUE)

  # find which layer contains the base
  # - uses a very slightly elevated value to avoid NA
  adj <- min(-diff(mfldivs), na.rm = TRUE)/100
  lbot <- cellref.loc(z0 + adj, rev(mfldivs), TRUE)

  used.mfldivs <- c(avwt, mfldivs[mfldivs < avwt & mfldivs > z0 + adj], z0)

  c(Map(function(h, n) rep(h/n, n),
        diff(used.mfldivs), ndlpmfl[ltop:lbot]), recursive = TRUE)
}
