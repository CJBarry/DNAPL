# DNAPL package - the DNAPL source term solver

#' DNAPL source term model
#'
#' @param result.file
#' character string;
#' file to save results to
#' @param description
#' character string;
#' information about the model and what it represents; this is saved with
#'  the results for your record
#' @param DNAPLmodel
#' DNAPLmodel S4 object;
#' see \code{\link{DNAPLmodel}}
#' @param uhist
#' data.frame:\cr
#' \code{$year} (num)\cr
#' \code{$cons} (num): consumption rate\cr
#' the transient import rate of solvent to the site (usage history)
#' @param uhist.J_units
#' character string;
#' units of flux (or usage rate) in the \code{$cons} column of
#'  \code{uhist}: \code{"kg/day"}, (default), \code{"kg/year"} or
#'  \code{"t/year"} (metric tonnes per year)
#' @param fsphist
#' data.frame:\cr
#' \code{$year} (num)\cr
#' \code{$f} (num): fraction of imported solvent that is spilt, rather than
#'  disposed of safely\cr
#' the transient safe disposal fraction of solvent by the site
#' @param fw
#' numeric \code{[1]};
#' proportion of imported solvent that becomes waste (between 0 and 1, 1 is
#'  suggested)
#' @param fi
#' numeric \code{[1]};
#' proportion of spilt solvent that infiltrates to water table (between 0
#'  1, 1 is suggested)
#' @param x,y
#' x and y co-ordinates of the spill with the same origin and units as
#'  \code{mfdata}
#' @param mfdata
#' NetCDF object;
#' MODFLOW data in NetCDF format (see \code{\link{GW.nc, package = Rflow}}),
#'  from which the transient horizontal Darcy velocity through the source
#'  zone is determined
#' @param qh
#' 1-layer: function(t) or numeric \code{[1]};\cr
#' multi-layer: list of function(t), one per layer or numeric \code{[NLAY]};\cr
#' if \code{mfdata}, the horizontal Darcy velocity at the source zone may
#'  be given explicity as a function of time or a single constant value (or
#'  list of functions, one per layer, or numeric vector of values per layer
#'  if the model is multi-layer)
#'
#' @return
#'
#' @import Rflow
#' @import RNetCDF
#' @importFrom stats approxfun
#' @importFrom rlist list.save
#' @importFrom grDevices xy.coords
#' @export
#'
#' @examples
DNST <- function(result.file, description = "", DNAPLmodel,
                 uhist, uhist.J_units = "kg/day",
                 fsphist = data.frame(year = c(1974, 1990), f = c(1, 0)),
                 fw = 1, fi = 1,
                 x, y = NULL, mfdata, qh = NULL){

  # determine flux to water table as function of time ----
  #
  # --------------------------------------------------------------------- #
  # The output of this section is a function of time, representing the
  #  solvent flux to the top of the water column, based on the import rate
  #  (uhist), waste fraction (fw), unsafe spillage fraction (fsphist) and
  #  infiltration fraction (fi).
  #
  # Jin: function(t)
  # --------------------------------------------------------------------- #
  #
  # - convert units of data frames
  uhist$t <- vapply(uhist$year, Rflow:::td, numeric(1L), d = 1, m = 7)
  fsphist$t <- vapply(fsphist$year, Rflow:::td, numeric(1L), d = 1, m = 7)
  J.mult <- switch(uhist.J_units,
                   "kg/day" = 1,
                   "kg/year" = 1/365.25,
                   "t/year" = 1000/365.25,
                   stop("invalid units for uhist.J_units: accepts \"kg/day\", \"kg/year\" or \"t/year\""))
  uhist$cons <- uhist$cons*J.mult
  #
  # - make into functions of time
  uhist <- approxfun(uhist, yleft = 0, yright = 0)
  fsphist <- approxfun(sdhist, rule = 2L)
  #
  # - apply losses from exports, safe disposal and volatilisation (not
  #    infiltrating)
  Jin <- function(t) uhist(t)*fi*fw*fsphist(t)
  # --------------------------------------------------------------------- #


  # check the DNAPL model ----
  #
  # --------------------------------------------------------------------- #
  # Checks whether the DNAPL model has the correct features with the
  #  correct structure and is internally consistent.  Stops if the DNAPL
  #  model is not correct and issues a message indicating what is wrong.
  # --------------------------------------------------------------------- #
  #
  DNAPLmodel.check(DNAPLmodel)
  # --------------------------------------------------------------------- #


  # determine transient Darcy velocity through each layer of the source term ----
  #
  # --------------------------------------------------------------------- #
  #
  #
  # qh: list with length NLAY of functions of time
  # --------------------------------------------------------------------- #
  #
  qh <- if(missing(mfdata)){

  }else DNAPL.qh(...)
  # --------------------------------------------------------------------- #


  # pre-allocate arrays for the results ----
  #
  # --------------------------------------------------------------------- #
  #
  # --------------------------------------------------------------------- #
  #


  # run model ----
  #
  # --------------------------------------------------------------------- #
  #
  # --------------------------------------------------------------------- #
  #


  # post-process results ----
  #
  # --------------------------------------------------------------------- #
  #
  # --------------------------------------------------------------------- #
  #


  # save results ----
  #
  # --------------------------------------------------------------------- #
  #
  # --------------------------------------------------------------------- #
  #
}