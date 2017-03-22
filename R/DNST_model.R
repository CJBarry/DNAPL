# DNAPL package - making DNAPL distribution and dissolution models

#' DNAPL distribution and dissolution model
#'
#' This is the structure of a DNAPL model used by \code{\link{DNST}} to
#'  specify the behaviour and dissolution of DNAPLs layer by layer.
#'
#' @slot NLAY integer;
#' number of layers in the DNAPL model
#' @slot hL numeric;
#' height of each DNAPL model layer
#' @slot params list;
#' a list of input parameters specific to the model
#' @slot domains vector of character strings;
#' names of the mass domains, one of which will probably be "plume"; other
#'  common examples include "NAPL" (lumped NAPL), "pool", "ganglia",
#'  "vapour", "fractureNAPL", ...
#' @slot domain1 character string;
#' the name of the domain to which the infiltrating DNAPL is sent initially;
#'  must feature in \code{domains}
#' @slot mdmax numeric matrix;
#' with rows representing layers and columns (named) representing domains,
#'  layer-by-layer values for the maximum mass in each domain; Inf is
#'  permitted and is often used for the plume domain
#' @slot mdredist list;
#' elements should be of class \code{\link{Mredistribution}}; number of
#'  elements is arbitrary as each element represents a mass redistribtion
#'  process
#' @slot spill.to data.frame;
#' one named row for each domain, except for domains with infinite capacity
#'  (as in \code{mdmax}) in every layer; this data frame redirects mass
#'  that exceeds the capacity for a domain in a given layer to a new domain
#'  or a new layer.  The data frame should have columns:\cr
#'  \code{$domain} chr/ factor: domain to which overspilt mass from this
#'   domain (in the \code{row.names}) should be directed; may be itself if
#'   \code{layer != 0}; for example excess mass in a "ganglia" domain may
#'   spill to the "pool" domain in the same layer, (in which case
#'   \code{layer} would be 0), or excess mass in a "NAPL" domain may spill
#'   to the "NAPL" domain in the underlying layer (in which case
#'   \code{layer} would be 1)
#'  \code{$layer} int: 1 for directing to the layer below; 0 for the same
#'   layer or (rarely) -1 for the layer above
#'  optional \code{$domain2} chr/ factor: a second choice domain
#'
#' @return
#' @export
#'
#' @examples
DNAPLmodel <- setClass("DNAPLmodel",
                       slots = c(NLAY = "integer",
                                 hL = "numeric",
                                 params = "list",
                                 domains = "character",
                                 domain1 = "character",
                                 mdmax = "matrix",
                                 mdredist = "list",
                                 spill.to = "data.frame"))

#' Mass Redistribution function
#'
#' @slot from character string;
#' the domain from which mass is taken
#' @slot to character string;
#' the domaint to which mass is sent
#' @slot flux function;
#' a function whose arguments are:\cr
#' \code{fromM}: mass of \code{from} domain\cr
#' \code{toM}: mass of \code{to} domain\cr
#' \code{LAY}: the layer of the DNAPL model\cr
#' \code{time}: in days since 30/12/1899\cr
#' and whose output represents the rate of mass transfer between the
#'  \code{from} and \code{to} domains (negative outputs are permitted and
#'  imply a flux the other way; in this way an equilibrium may be set up)
#'
#' @details
#' The \code{flux} function must have all four arguments even if in fact
#'  they are not all used.  The function may also use global variables that
#'  will exist in the \code{\link{DNST}} solver function.  Typically,
#'  \code{qh} will feature somewhere in some functions, as the dissolution
#'  mass flux, for example, will depend on the water speed through the
#'  source zone.  Some functions may also use the \code{M} array in order
#'  to include some dependency on the source zone history.
#'
#' @return
#' @export
#'
#' @examples
Mredistribution <- setClass("Mredistribution",
                            slots = c(from = "character", to = "character",
                                      flux = "function"),
                            contains = c("function"))

#' Constant power mass-flux relationship model
#'
#' @param wg
#' numeric \code{[1]} or \code{[NLAY]};
#' ganglia width
#' @param wpm
#' numeric \code{[1]} or \code{[NLAY]};
#' maximum pool width (may be \code{Inf} in the lowest layer, in which case
#'  no DNAPL is allowed to escape)
#' @param hp
#' numeric \code{[1]} or \code{[NLAY]};
#' pool height, or total height of pools in a layer
#' @param Gamma
#' numeric \code{[1]} or \code{[NLAY]};
#' empirical source depletion parameter, positive; small values (<1) imply
#'  a more persistent source term
#' @param Srn,Srw
#' numeric \code{[1]};
#' residual saturations of NAPL and water
#' @param phi
#' numeric \code{[1]} or \code{[NLAY]};
#' bulk porosity of each layer
#' @param rho
#' numeric \code{[1]};
#' solvent density (ensure consistent units, probably kg/m^3)
#' @param Cs
#' numeric \code{[1]};
#' solvent solubility in water (ensure consistent units, probably kg/m^3
#'  which is the same as g/l)
#' @param hL
#' numeric \code{[1]} or \code{[NLAY]};
#' the height of each layer; note that higher layers will contain more NAPL
#'  mass
#' @param NLAY
#' integer \code{[1]};
#' number of layers in the DNAPL model
#'
#' @return
#'
#' @import methods
#' @export
#'
#' @examples
cstG.DNmodel <- function(wg, wpm, hp, Gamma, Srn, Srw, phi, rho, Cs, hL,
                         NLAY = length(hL)){

  # expand vectors where necessary
  hL <- expand.vec(hL, NLAY)
  wg <- expand.vec(wg, NLAY)
  wpm <- expand.vec(wpm, NLAY)
  hp <- expand.vec(hp, NLAY)
  Gamma <- expand.vec(Gamma, NLAY)
  phi <- expand.vec(phi, NLAY)

  # mass capacity and area calculations
  # - mass per unit horizontal area in a pool
  #  -- assumes that the pool has saturation (1 - Srw) at the base,
  #      exponentially declining to Srn at the top (hp); Sbar is the mean
  #      saturation of the column
  l <- log(Srn/(1 - Srw))/hp
  Sbar <- (2*Srn)/(l*hp^2)*(hp*exp(l*hp) - (exp(l*hp) - 1)/l)
  mpua <- phi*rho*Sbar*hp
  #
  # - flow-normal ganglia and pool areas and maximum masses
  Aqg <- wg*hL
  mgmax <- hL*((pi*wg^2)/4)*phi*Srn*rho
  mpmax <- wpm^2*pi/4*mpua
  #
  # - maximum mass for each domain - no limit for plume
  mdmax <- cbind(NAPL = mgmax + mpmax, plume = Inf)

  # mass redistribution function
  # - dissolution
  # - uses qh and M which will be present in the solution function
  #    environment
  mdredist <- list(Mredistribution(from = "NAPL", to = "plume", flux = {
    function(fromM, toM, LAY, time){
      # value of m0 depends on saturation history
      m0 <- max(M[LAY,, "NAPL"])
      if(m0 == 0) return(0)

      # horizontal area of pool, hence width of pool assuming circular
      Abas <- m0/mpua
      wp <- sqrt(4*Abas/pi)

      # flow-normal area of pool
      Aqp <- (wp*hp)[LAY]

      Rc <- (1 - fromM/m0)^Gamma[LAY]
      C <- Cs*(1 - Rc)
      qh[[LAY]](time)*C*(Aqg[LAY] + Aqp)
    }
  })
  )

  # create DNAPL model S4 object
  DNAPLmodel(NLAY = as.integer(NLAY),
             hL = as.numeric(hL),
             params = mget(c("wg", "wpm", "hp", "Srn", "Srw",
                             "phi", "rho", "Cs", "Gamma",
                             "mpua")),
             domains = c("NAPL", "plume"),
             domain1 = c("NAPL"),
             mdmax = mdmax,
             mdredist = mdredist,
             spill.to = data.frame(row.names = "NAPL",
                                   domain = "NAPL",
                                   layer = 1L))
}
