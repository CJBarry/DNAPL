# DNAPL package - making DNAPL distribution and dissolution models

#' DNAPL distribution and dissolution model
#'
#' This is the structure of a DNAPL model used by \code{\link{DNST}} to
#'  specify the behaviour and dissolution of DNAPLs layer by layer.
#'
#' @details
#' Use the \code{\link{DNAPLmodel.check}} function for a thorough analysis
#'  of the DNAPL model, to make sure that it is set up correctly.
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
#'   \code{layer} would be 1)\cr
#'  \code{$layer} int: 1 for directing to the layer below; 0 for the same
#'   layer or (rarely) -1 for the layer above\cr
#'  optional \code{$domain2} chr/ factor: a second choice domain
#'
#' @export DNAPLmodel
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
#' A description of a process that transfers mass from one domain to
#'  another.  A common example is the dissolution of NAPL into groundwater.
#'  In this case, slot \code{from} may be \code{"NAPL"}, slot \code{to} may
#'  be \code{"plume"} and slot \code{flux} will be a function which returns
#'  the rate of mass transfer (dissolution) from the NAPL domain into the
#'  plume domain, as a function of the NAPL mass.
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
#' \code{env}: the calling function's environment\cr
#' other parameters used within the function that are not got from the DNST
#'  environment (see example)\cr
#' and whose output represents the rate of mass transfer between the
#'  \code{from} and \code{to} domains (negative outputs are permitted and
#'  imply a flux the other way; in this way an equilibrium may be set up)
#'
#' @details
#' The \code{flux} function must at least have the five arguments named
#'  above, as well as any arguments for the parameters within the function,
#'  even if in fact they are not all used.  The function may also use
#'  variables that will exist in the \code{\link{DNST}} solver function by
#'  using, for example, \code{get("qh", env)}. Typically, \code{qh} will
#'  feature somewhere in some functions, as the dissolution mass flux, for
#'  example, will depend on the water speed through the source zone.
#'  \code{qh} is given to \code{\link{DNST}} as a layer-by-layer function of
#'  time, so it is common to use \code{get("qh", env)[[LAY]](time)} in the
#'  \code{flux} function. Some functions may also use the \code{M} array in
#'  order to include some dependency on the source zone history.  When
#'  included as part of a \code{\link{DNAPLmodel}} (as intended), the
#'  function may use any of the parameters saved in the slot \code{params}
#'  as well. Remember to subset these parameters by layer
#'  (\code{par[[LAY]]}) when the parameters could vary by layer.
#'
#' When used in \code{\link{DNST}}, the solver function, no flux is allowed
#'  to reduce the mass of a domain below 0, so the solution is always
#'  stable.  Therefore, this needn't be worried about in the \code{flux}
#'  function.
#'
#' \code{Mredistribution} objects are created within the provided
#'  \code{\link{DNAPLmodel}} creation functions
#'  (\code{\link{cstG.DNmodel}}, \code{\link{cnvG.DNmodel}} and
#'  \code{\link{DDpg.DNmodel}} currently).  Defining your own
#'  \code{Mredistribution} is fairly advanced and will require some
#'  understanding of environments in \code{R}.
#'
#' @export Mredistribution
#'
#' @examples
#' # exponentially declining dissolution of a NAPL pool
#' \dontrun{
#' Ndiss <- Mredistribution(from = "pool",
#'                          to = "plume",
#'                          flux = function(fromM, toM, LAY, time, env,
#'                                          pool.width, pool.height, solubility){
#'                            # get historic maximum of mass from this layer
#'                            # - use get("M", env) to access the mass array from
#'                            #    the DNST master function which will call this
#'                            #    function
#'                            max.mass.so.far <- max(get("M", env)[LAY,, "pool"])
#'
#'                            # to use this redistribution model, pool.width and
#'                            #  pool.height will need to be defined in the global
#'                            #  environment or, better, included in the params
#'                            #  slot of the DNAPLmodel object which contains this
#'                            #  Mredistribution object
#'                            Xsect.area <- pool.width*pool.height
#'
#'                            # again, horizontal flow (Darcy velocity) is accessed
#'                            #  by get("qh", env)
#'                            # solubility must be defined in the same way as
#'                            #  pool.width and pool.height (see comment above)
#'                            solubility*Xsect.area*get("qh", env)[[LAY]](time)*fromM/max.mass.so.far
#'                          })
#' }
#'
Mredistribution <- setClass("Mredistribution",
                            slots = c(from = "character", to = "character",
                                      flux = "function"))

#' Check that a DNAPL model is set up correctly
#'
#' The set up of a DNAPL distribution and dissolution model is fairly
#'  detailed (see \code{\link{DNAPLmodel}}).  This function provides an
#'  easy way to check that a DNAPL model is set up correctly.  If it isn't,
#'  a vector of messages is returned that inform the user of all the things
#'  that need to be corrected.  A DNAPL model that passes this check should
#'  work with \code{\link{DNST}}.
#'
#' @note
#' Some more obscure errors will not be picked up by this function.  In the
#'  \code{spill.to} slot, it is possible to define a cycle of overspilling
#'  (e.g. if domain A spills to domain B and domain B spills to domain A,
#'  all in the same layer), which could result in an infinite loop during
#'  model calculation.  This error would not be detected by
#'  \code{DNAPLmodel.check}.
#'
#' @param d
#' DNAPLmodel object to be checked, as given by \code{\link{DNAPLmodel}}
#' @param stop.if.mistakes
#' logical \code{[1]};
#' should R stop executing if there are any mistakes (default \code{TRUE})
#'
#' @return
#' A character vector explaining all mistakes in the DNAPL model
#'  (\code{d}).  If there are any mistakes, the session will stop if
#'  \code{stop.if.mistakes == TRUE} and the mistakes will be printed as
#'  error messages.
#'
#' @importFrom methods is
#' @export
#'
DNAPLmodel.check <- function(d, stop.if.mistakes = TRUE){
  # for storing mistake messages
  m <- character(0L)

  # check class
  if(!isS4(d) || !is(d, "DNAPLmodel")) m <- c(m, {
    "the DNAPL model should be an S4 object with class 'DNAPLmodel': see help(DNAPLmodel)"
  })

  # check NLAY
  if(length(d@NLAY) != 1L) m <- c(m, {
    "The length of the NLAY slot should be 1: this slot gives the number of layers in the model."
  })

  # check no duplicated domains
  if(any(duplicated(d@domains))) m <- c(m, {
    "duplicated domain names detected"
  })

  # check domain 1
  d1 <- d@domain1
  tst <- length(d1) == 1L && d1 %in% d@domains
  if(!tst) m <- c(m, {
    "domain1 should be a single character string matching one of the named domains in the domains slot;\ndomain1 represents the mass domain to which input mass to layer 1 is directed and is likely to be a NAPL domain"
  })

  # check mdmax
  mdm <- d@mdmax
  expl <- "mdmax represents the maximum mass in each domain in each layer, with named columns representing mass domains and rows representing layers"
  #
  # - check columns
  subject <- colnames(mdm)
  tst <- ncol(mdm) == length(d@domains) &&
    all(d@domains %in% subject) &&
    all(subject %in% d@domains)
  if(!tst) m <- c(m, paste({
    paste0("There should be one named column (see help(colnames)) in the mdmax slot for each of the named domains in the domains slot.  Current colnames are: ", paste(subject, collapse = ", "), ", but should be: ", paste(d@domains, collapse = ", "), " (in any order).")
  }, expl, sep = "\n"))
  #
  # - check rows
  tst <- identical(nrow(mdm), d@NLAY)
  if(!tst) m <- c(m, paste({
    "The mdmax slot must have one row for each layer of the model (slot NLAY)."
  }, expl, sep = "\n"))

  # check spill.to
  st <- d@spill.to
  expl <- "The spill.to data frame governs what happens when the mass in each domain (in the row.names, see help(row.names)) exceeds the capacity in a given layer: which domain does the excess mass go into and does it go into the layer below?  See help(DNAPLmodel) for more details."
  #
  # - check rownames
  #  -- only needs to include finite-capacity domains (not the plume, for
  #      example)
  nec.rns <- d@domains[vapply(d@domains, function(dm){
    dm %in% colnames(d@mdmax) && !all(is.infinite(d@mdmax[, dm]))
  }, logical(1L))]
  tst <- all(nec.rns %in% row.names(st)) && !any(duplicated(row.names(st)))
  if(!tst) m <- c(m, paste({
    paste0("The slot spill.to must have one entry (row) for each mass domain, apart from domains which can contain infinite mass (such as the plume, normally).  The rows must be labelled with the domain it refers to using the row.names attribute (see help(row.names)).  The current row.names are: ", paste(row.names(st), collapse = ", "), ", but they should be: ", paste(nec.rns, collapse = ", "), " (in any order).")
  }, expl, sep = "\n"))
  #
  # - check structure
  tst <- `||`(identical(sapply(st, class),
                        c(domain = "factor", layer = "integer")),
              identical(sapply(st, class),
                        c(domain = "factor", layer = "integer",
                          domain2 = "factor")))
  if(!tst) m <- c(m, paste({
    "The columns in spill.to are not of the correct type.  They must be $domain: factor (input as character vector); $layer: integer (note, not double; use 0L, for example to denote an integer on input, rather than simply 0).  An optional $domain2 (factor from character vector) column is also permitted: see help(DNAPLmodel)."
  }, expl, sep = "\n"))
  #
  # - check logic
  #  -- not a complete check, just checks that a domain cannot spill to
  #      itself in the same layer, although spilling to the same domain in
  #      the layer below is quite common
  tst <- !any(prob <- (st$domain == row.names(st) & st$layer == 0L))
  if(!tst) m <- c(m, paste({
    paste0("The destination domain for overspills (domain column) can only be the same domain as the spilling domain (in the row.names) if it is to a different layer (layer column is -1L or 1L for that row).  Problem rows are: ", paste(row.names(st)[prob], collapse = ", "))
  }, expl, "\n"))

  # check mdredist
  mdr <- d@mdredist
  expl <- "The mdredist slot is a list of mass transfer processes between domains describing, for example, dissolution of each NAPL domain, evaporation, matrix diffusion, as necessary for the particular model.  Each element of the list should be a Mredistritution S4 object (see help(Mredistribution))."
  #
  # - check that mdredist is a list of exclusively Mredistribution objects
  tst <- all(sapply(mdr, class) == "Mredistribution")
  if(!tst) m <- c(m, paste({
    "Some or all elements of slot mdredist are not Mredistribution objects."
  }, expl, sep = "\n"))
  #
  # - check Mredistribution set ups
  if(tst){
    # no point in checking this if the previous test failed
    tst <- any(okay <- vapply(mdr, function(x){
      class(x@to) == "character" && length(x@to) == 1L &&
        class(x@from) == "character" && length(x@from) == 1L &&
        identical(names(formals(x@flux))[1:5],
                  c("fromM", "toM", "LAY", "time", "env"))
    }, logical(1L)))
    if(!tst) m <- c(m, paste({
      paste0("The Mredistribution objects in slot mdredist are not all set up correctly.  Items: (", paste(which(!okay), collapse = ", "), ") are incorrect.  See help(Mredistribution).")
    }, expl, sep = "\n"))
  }

  # finish up and return/ stop
  if(length(m) && stop.if.mistakes) stop(paste(m, collapse = "\n"))
  cat(m, sep = "\n")
  invisible(m)
}

#' @rdname cG
#' @name cst-cnvG
#'
#' @title Constant and Convering power mass-flux relationship model
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
#' @param Gamma,Gamma0
#' numeric \code{[1]} or \code{[NLAY]};
#' empirical source depletion parameter, positive; small values (<1) imply
#'  a more persistent source term\cr
#' for \code{cstG.DNmodel}, \Gamma is constant throughout the model,
#'  but for \code{cnvG.DNmodel}, \Gamma converges linearly from
#'  \code{Gamma0} (at peak mass) to 1 as the NAPL mass depletes.
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
#' a \link{DNAPLmodel} S4 object
#'
NULL

#' @rdname cG
#'
#' @import methods
#' @export
#'
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
  mpua <- pool_mpua(hp, phi, rho, Srn, Srw, "exponential")
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
    function(fromM, toM, LAY, time, env, mpua, hp, Gamma, Aqg, Cs){
      # value of m0 depends on saturation history
      m0 <- max(get("M", env)[LAY,, "NAPL"])
      if(m0 == 0) return(0)

      # horizontal area of pool, hence width of pool assuming circular
      Abas <- m0/mpua
      wp <- sqrt(4*Abas/pi)

      # flow-normal area of pool
      Aqp <- (wp*hp)[LAY]

      Rc <- (1 - fromM/m0)^Gamma[LAY]
      C <- Cs*(1 - Rc)
      get("qh", env)[[LAY]](time)*C*(Aqg[LAY] + Aqp)
    }
  }))

  # create DNAPL model S4 object
  DNAPLmodel(NLAY = as.integer(NLAY),
             hL = as.numeric(hL),
             params = mget(c("wg", "wpm", "hp", "Srn", "Srw",
                             "phi", "rho", "Cs", "Gamma",
                             "Aqg", "mpua")),
             domains = c("NAPL", "plume"),
             domain1 = c("NAPL"),
             mdmax = mdmax,
             mdredist = mdredist,
             spill.to = data.frame(row.names = "NAPL",
                                   domain = "NAPL",
                                   layer = 1L))
}

#' @rdname cG
#'
#' @import methods
#' @export
#'
cnvG.DNmodel <- function(wg, wpm, hp, Gamma0, Srn, Srw, phi, rho, Cs, hL,
                         NLAY = length(hL)){

  # expand vectors where necessary
  hL <- expand.vec(hL, NLAY)
  wg <- expand.vec(wg, NLAY)
  wpm <- expand.vec(wpm, NLAY)
  hp <- expand.vec(hp, NLAY)
  Gamma0 <- expand.vec(Gamma0, NLAY)
  phi <- expand.vec(phi, NLAY)

  # mass capacity and area calculations
  # - mass per unit horizontal area in a pool
  #  -- assumes that the pool has saturation (1 - Srw) at the base,
  #      exponentially declining to Srn at the top (hp); Sbar is the mean
  #      saturation of the column
  mpua <- pool_mpua(hp, phi, rho, Srn, Srw, "exponential")
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
    function(fromM, toM, LAY, time, env, mpua, hp, Gamma0, Aqg, Cs){

      # value of m0 depends on saturation history
      m0 <- max(get("M", env)[LAY,, "NAPL"])
      if(m0 == 0) return(0)

      # horizontal area of pool, hence width of pool assuming circular
      Abas <- m0/mpua
      wp <- sqrt(4*Abas/pi)

      # flow-normal area of pool
      Aqp <- (wp*hp)[LAY]

      # determine Gamma based on proportion of mass removed
      Rm <- 1 - fromM/m0
      Gamma <- Gamma0[LAY] + Rm*(1 - Gamma0[LAY])

      Rc <- (1 - fromM/m0)^Gamma
      C <- Cs*(1 - Rc)
      get("qh", env)[[LAY]](time)*C*(Aqg[LAY] + Aqp)
    }
  })
  )

  # create DNAPL model S4 object
  DNAPLmodel(NLAY = as.integer(NLAY),
             hL = as.numeric(hL),
             params = mget(c("wg", "wpm", "hp", "Srn", "Srw",
                             "phi", "rho", "Cs", "Gamma0",
                             "Aqg", "mpua")),
             domains = c("NAPL", "plume"),
             domain1 = c("NAPL"),
             mdmax = mdmax,
             mdredist = mdredist,
             spill.to = data.frame(row.names = "NAPL",
                                   domain = "NAPL",
                                   layer = 1L))
}

#' Dual-domain Pool and Ganglia model
#'
#' This model has separate domains for ganglia and pools.  Ganglia are
#'  flushed efficiently (effectively \code{Gamma = Inf}) and pools are
#'  flushed inefficiently, with a \code{poolGamma} depletion parameter that
#'  defaults to 1.
#'
#' @inheritParams cstG.DNmodel
#' @param poolGamma
#' numeric \code{[1]} or \code{[NLAY]};
#' the empirical depletion power to be applied to pools, analogous to
#'  \code{Gamma} in \code{\link{cstG.DNmodel}} which applies to the bulk
#'  NAPL
#'
#' @return
#' a \link{DNAPLmodel} S4 object
#'
#' @import methods
#' @export
#'
DDpg.DNmodel <- function(wg, wpm, hp, Srn, Srw, phi, rho, Cs, hL,
                         NLAY = length(hL), poolGamma = 1){

  # expand vectors where necessary
  hL <- expand.vec(hL, NLAY)
  wg <- expand.vec(wg, NLAY)
  wpm <- expand.vec(wpm, NLAY)
  hp <- expand.vec(hp, NLAY)
  phi <- expand.vec(phi, NLAY)
  poolGamma <- expand.vec(poolGamma, NLAY)

  # mass capacity and area calculations
  # - mass per unit horizontal area in a pool
  #  -- assumes that the pool has saturation (1 - Srw) at the base,
  #      exponentially declining to Srn at the top (hp); Sbar is the mean
  #      saturation of the column
  mpua <- pool_mpua(hp, phi, rho, Srn, Srw, "exponential")
  #
  # - flow-normal ganglia and pool areas and maximum masses
  Aqg <- wg*hL
  mgmax <- hL*((pi*wg^2)/4)*phi*Srn*rho
  mpmax <- wpm^2*pi/4*mpua
  #
  # - maximum mass for each domain - no limit for plume
  mdmax <- cbind(pool = mpmax, ganglia = mgmax, plume = Inf)

  # mass redistribution functions
  mdredist <- list(Mredistribution(from = "ganglia", to = "plume", flux = {
    function(fromM, toM, LAY, time, env, Aqg, Cs){
      C <- ifelse(fromM == 0, 0, Cs)
      get("qh", env)[[LAY]](time)*C*Aqg[LAY]
    }
  }),
  Mredistribution(from = "pool", to = "plume", flux = {
    function(fromM, toM, LAY, time, env, mpua, hp, poolGamma, Cs){
      # in this model, max. mass is based on saturation history, so there
      #  can be a limitless
      m0 <- max(get("M", env)[LAY,, "pool"], fromM)

      # if no mass has ever entered this cell
      if(m0 == 0) return(0)

      # horizontal area of pool, hence width of pool assuming circular
      Abas <- m0/mpua
      wp <- sqrt(4*Abas/pi)

      # determine pool Gamma based on proportion of mass removed
      Rm <- 1 - fromM/m0

      Rc <- (1 - fromM/m0)^poolGamma[LAY]
      C <- Cs*(1 - Rc)

      # flow-normal area of pool
      Aqp <- (wp*hp)[LAY]

      get("qh", env)[[LAY]](time)*C*Aqp
    }
  }))

  # create DNAPL model S4 object
  DNAPLmodel(NLAY = as.integer(NLAY),
             hL = as.numeric(hL),
             params = mget(c("wg", "wpm", "hp", "Srn", "Srw",
                             "phi", "rho", "Cs", "poolGamma",
                             "Aqg", "mpua")),
             domains = c("ganglia", "pool", "plume"),
             domain1 = c("ganglia"),
             mdmax = mdmax,
             mdredist = mdredist,
             spill.to = data.frame(row.names = c("ganglia", "pool"),
                                   domain = c("pool", "ganglia"),
                                   layer = c(0L, 1L)))
}

#' Mass per unit area of pool
#'
#' @inheritParams cstG.DNmodel
#' @param mode
#' character [1];
#' nature of saturation variation with height above pool base, one of
#'  \code{"uniform"}, \code{"linear"} or \code{"exponential"}; for
#'  \code{"uniform"}, it is taken that the pool has saturation
#'  \code{1 - Srw} and for the other two it is taken that the pool has
#'  saturation \code{1 - Srw} at the base and \code{Srn} at the top
#'
#' @return
#' numeric[]
#'
pool_mpua <- function(hp, phi, rho, Srn, Srw,
                      mode = c("uniform", "linear", "exponential")[1L]){
  Sbar <- switch(mode,
                 uniform = 1 - Srw,
                 linear = (1 - Srw + Srn)/2,
                 exponential = {
                   # see a3.p76
                   var <- log((1 - Srw)/Srn)

                   (1 - Srw)/(var)*(1 - exp(-var))
                 },
                 stop("DNAPL::pool_mpua: invalid mode; choose \"uniform\", \" linear\" or \"exponential\""))

  Sbar*phi*hp*rho
}
