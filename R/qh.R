# DNAPL package - calculation of Darcy velocity

#' Horizontal Darcy Velocity
#'
#' Calculates the horizontal Darcy velocity, or specific discharge, at a
#'  specific location in a MODFLOW model.
#'
#' @details
#' Uses linear interpolation to calculate the velocity between cell faces.
#'  For the layer with the water table, the thickness (and therefore flow
#'  normal area) is transient.  NA is returned for dry layers.
#'
#' @param x,y
#' location, read by \code{xy.coords};
#' should be in absolute co-ordinates with the same origin and units as
#'  \code{mfdata}
#' @param L
#' integer \code{[]};
#' a continuous range is insisted upon and this is corrected to
#'  \code{min(L):max(L)} within function
#' @param mfdata
#' NetCDF object;
#' a MODFLOW dataset; see \code{\link[Rflow]{GW.nc}}
#'
#' @return
#' numeric array \code{[length(L), NTS]}, dropped, where \code{NTS} is the
#'  number of time steps in the MODFLOW results;
#' the Darcy velocity in each layer for each time step at the location x,y
#'
#' @import RNetCDF
#' @import Rflow
#' @importFrom grDevices xy.coords
#' @importFrom abind adrop
#' @importFrom stats qunif
#' @importFrom stats punif
#' @export
#'
#' @examples
#' library("RNetCDF")
#'
#' # get example MODFLOW data set
#' mfdata <- open.nc(system.file("rflow_mf_demo.nc", package = "Rflow"))
#' qh(625, 825, 1L, mfdata)
#'
qh <- function(x, y = NULL, L, mfdata){
  xy <- xy.coords(x, y)
  x <- xy$x; y <- xy$y

  # ensure L is a continuous integer sequence
  L <- min(L):max(L)
  nlay <- length(L)

  # initialise array
  nts <- dim.inq.nc(mfdata, "NTS")$length
  q <- array(NA_real_, c(nlay, nts),
             dimnames = list(paste0("L", L), NULL))

  # column and row reference
  ccs <- gccs(mfdata, TRUE)
  rcs <- grcs(mfdata, TRUE)
  C <- cellref.loc(x, ccs, FALSE)
  R <- cellref.loc(y, rcs, TRUE)

  # layer bottoms
  bots <- adrop(var.get.nc(mfdata, "elev",
                           c(C, R, L[1L] + 1L), c(1L, 1L, nlay),
                           collapse = FALSE), c(TRUE, TRUE, FALSE))

  # transient layer tops
  # - cell tops
  ctops <- adrop(var.get.nc(mfdata, "elev",
                            c(C, R, L[1L]), c(1L, 1L, nlay),
                            collapse = FALSE), c(TRUE, TRUE, FALSE))
  #
  # - transient water table
  wtab <- adrop(var.get.nc(mfdata, "Head",
                           c(C, R, L[1L], NA), c(1L, 1L, nlay, NA),
                           collapse = FALSE), c(TRUE, TRUE, FALSE, FALSE))
  #
  # - whichever is lower out of cell top or water table or, if the water
  #    table is below the cell bottom, then NA
  tops <- array(NA_real_, c(nlay, nts))
  for(ts in 1L:nts){
    tops[, ts] <- ifelse(wtab[, ts] > bots,
                         ifelse(wtab[, ts] > ctops, ctops, wtab[, ts]),
                         NA_real_)
  }

  # transient saturated heights
  # - explicit form would be:
  # thk <- array(NA_real_, c(nlay, nts))
  # for(ts in 1:nts) thk[, ts] <- tops[, ts] - bots
  #
  # - can be vectorised because of favourable dimension ordering
  thk <- tops - bots

  # cell face widths
  # - right face
  wx <- diff(rev(rcs[R + 1:0]))
  #
  # - front face
  wy <- diff(ccs[C + 0:1])

  # DV in x direction
  # - DV through left and right face
  if(C == 1L){
    flrf <- adrop(var.get.nc(mfdata, "FlowRightFace",
                             c(C, R, L[1L], NA), c(1L, 1L, nlay, NA),
                             collapse = FALSE), c(FALSE, TRUE, FALSE, FALSE))
    dvlrf <- array(NA_real_, c(2L, nlay, nts))
    dvlrf[1L,,] <- 0
    dvlrf[2L,,] <- flrf[1L,,]/(wx*thk)
  }else{
    flrf <- adrop(var.get.nc(mfdata, "FlowRightFace",
                             c(C - 1L, R, L[1L], NA), c(2L, 1L, nlay, NA),
                             collapse = FALSE), c(FALSE, TRUE, FALSE, FALSE))
    dvlrf <- array(NA_real_, c(2L, nlay, nts))
    dvlrf[1L,,] <- flrf[1L,,]/(wx*thk)
    dvlrf[2L,,] <- flrf[2L,,]/(wx*thk)
  }
  dvlrf[is.na(dvlrf)] <- 0
  #
  # - fractional position from left to right of cell
  fC <- punif(x, ccs[C], ccs[C + 1L])
  #
  # - linearly interpolate qx
  #  -- qunif cannot work backwards (max < min), so use this solution to
  #      detect when the limits should be swapped and 1 - fC used
  lgtr <- dvlrf[1L,,] > dvlrf[2L,,]
  qx <- qunif(ifelse(lgtr, 1 - fC, fC),
              ifelse(lgtr, dvlrf[2L,,], dvlrf[1L,,]),
              ifelse(lgtr, dvlrf[1L,,], dvlrf[2L,,]))

  # DV in -y direction
  # - DV through back and front face
  if(R == 1L){
    fbff <- adrop(var.get.nc(mfdata, "FlowFrontFace",
                             c(C, R, L[1L], NA), c(1L, 1L, nlay, NA),
                             collapse = FALSE),
                  c(TRUE, FALSE, FALSE, FALSE))
    dvbff <- array(NA_real_, c(2L, nlay, nts))
    dvbff[1L,,] <- 0
    dvbff[2L,,] <- fbff[1L,,]/(wy*thk)
  }else{
    fbff <- adrop(var.get.nc(mfdata, "FlowFrontFace",
                             c(C, R - 1L, L[1L], NA), c(1L, 2L, nlay, NA),
                             collapse = FALSE),
                  c(TRUE, FALSE, FALSE, FALSE))
    dvbff <- array(NA_real_, c(2L, nlay, nts))
    dvbff[1L,,] <- fbff[1L,,]/(wy*thk)
    dvbff[2L,,] <- fbff[2L,,]/(wy*thk)
  }
  dvbff[is.na(dvbff)] <- 0
  #
  # - fractional position from back to front of cell
  #  -- increases in opposite direction to y
  fR <- 1 - punif(y, rev(rcs)[R + 1L], rev(rcs)[R])
  #
  # - linearly interpolate qy
  #  -- qunif cannot work backwards (max < min), so use this solution to
  #      detect when the limits should be swapped and 1 - fR used
  bgtf <- dvbff[1L,,] > dvbff[2L,,]
  qy <- qunif(ifelse(bgtf, fR, 1 - fR),
              ifelse(bgtf, dvbff[2L,,], dvbff[1L,,]),
              ifelse(bgtf, dvbff[1L,,], dvbff[2L,,]))

  # Pythagoras
  result <- sqrt(qx^2L + qy^2L)
  if(nlay > 1L){
    structure(result, dimnames = list(paste0("L", L), NULL))
  }else result
}
