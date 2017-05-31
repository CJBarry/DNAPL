# DNAPL package - pool dissolution models

# these models give the C_0 from pools, that is, the concentration that
#  would have been expected when mass was at its historic maximum

#' @rdname pdm
#' @title Flow-normal cross-sectional area PDM
#'
#' @param Cs
#' numeric [];
#' solubility
#' @param hp
#' numeric [];
#' pool height
#' @param hL
#' numeric [];
#' layer thickness
#'
#' @return
#' numeric [];
#' C0, the historic maximum of effluent concentration averaged through the
#'  layer
#'
C0_pool_XsectA <- function(Cs, hp, hL) Cs*hp/hL

#' @rdname pdm
#' @title Hunt et al. 1988 PDM
#'
#' @description
#' integral of \code{Cs*erfc(z/(2*sqrt(ax*wp0)))} between 0 and \code{hp},
#'  or else lots of smaller pools if \code{Np} isn't 1
#'
#' @inheritParams C0_pool_XsectA
#' @param wp0
#' numeric [];
#' historic maximum of pool width
#' @param aV
#' numeric [1];
#' vertical dispersivity, in units of length
#' @param Np
#' integer [1];
#' number of pools within the layer; more pools will result in a greater
#'  average effluent concentration as no point will be so far from a pool
#'
#' @import pracma
#' @return
#' numeric [];
#' C0, the historic maximum of effluent concentration averaged through the
#'  layer
#'
#' @note
#' It is assumed that dispersion only occurs above the pool, as pools rest
#'  on top of aquitards, so hydrodynamic dispersion below the pool should
#'  be fairly small, instead just molecular diffusion should occur into the
#'  underlying aquitard.  If you want to dispense with this assumption,
#'  then double the value for \code{Np}, as the two-way dispersion from one
#'  pool at the bottom of the layer will result in the same average
#'  concentration as the one-way dispersion from two pools, one at the
#'  bottom and one in the middle.
#'
C0_pool_Hunt88 <- function(Cs, wp0, hL, aV = wp0/1000, Np = 1L){
  Nsegs <- 100L
  zpts <- 10^seq(log10(hL/Np) - 10, log10(hL/Np), length.out = Nsegs + 1L)

  # pi*wp0/4 is the average pool length across its width, assuming circular
  Csamples <- Cs*Hunt88_C.Csprof(zpts, pi*wp0/4, aV)

  Np*sum((Csamples[-1L] + Csamples[-(Nsegs + 1L)])/2*diff(zpts))/hL
}
Hunt88_C.Csprof <- function(z, L, aV) erfc(z/(2*(aV*L)^.5))
