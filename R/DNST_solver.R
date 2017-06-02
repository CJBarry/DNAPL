# DNAPL package - the DNAPL source term solver

#' DNAPL source term model
#'
#' @param result.file
#' character string;
#' file to save results to (should end with extension \code{".rds"})
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
#' @param start.t,end.t
#' numeric \code{[1]};
#' start and end time of the model, in days since 30/12/1899
#' @param dt
#' numeric \code{[1]};
#' time step length, in days, usually in the order of tens of days
#' @param x,y
#' x and y co-ordinates of the spill with the same origin and units as
#'  \code{mfdata};
#' used to find the grid cell from which to read MODFLOW results;
#'  alternatively, if \code{mfdata} is not used, x and y can simply be
#'  given as information to be stored with the model
#' @param z0
#' numeric \code{[1]} or \code{"base"};
#' elevation of the base of the DNAPL model, with the same datum as
#'  \code{mfdata}; \code{"base"} instructs to use the bottom elevation of
#'  the lowest layer at \code{x,y}, read from \code{mfdata}
#' @param mfdata
#' NetCDF object;
#' MODFLOW data in NetCDF format (see \code{\link[Rflow]{GW.nc}}),
#'  from which the transient horizontal Darcy velocity through the source
#'  zone is determined
#' @param qh
#' 1-layer: function(t) or numeric \code{[1]};\cr
#' multi-layer: list \code{[NLAY]} of function(t), one per layer or numeric
#'  \code{[NLAY]}, where \code{NLAY} is the number of layers in the
#'  \code{DNAPLmodel} (\code{DNAPLmodel@NLAY});\cr
#' if \code{mfdata} is missing, the horizontal Darcy velocity at the source
#'  zone may be given explicity as a function of time or a single constant
#'  value (or list of functions, one per layer, or numeric vector of values
#'  per layer if the model is multi-layer)
#' @param check.qh
#' logical \code{[1]};
#' whether to check the final set of qh functions for NA results (could be
#'  time-consuming for lots of time steps, but makes error catching easier)
#'
#' @return
#' A \code{\link{DNAPLSourceTerm}} S4 object
#'
#' @import Rflow
#' @import RNetCDF
#' @import methods
#' @import td
#' @importFrom stats approxfun
#' @importFrom rlist list.save
#' @importFrom grDevices xy.coords
#' @importFrom abind adrop
#' @export
#'
#' @examples
#' # to follow
#'
DNST <- function(result.file, description = "", DNAPLmodel,
                 uhist, uhist.J_units = "kg/day",
                 fsphist = data.frame(year = c(1974, 1990), f = c(1, 0)),
                 fw = 1, fi = 1, start.t = 9100, end.t, dt = 20,
                 x = NULL, y = NULL, z0 = "base", mfdata, qh = NULL,
                 check.qh = TRUE){

  # check the DNAPL model ----
  #
  # --------------------------------------------------------------------- #
  # Checks whether the DNAPL model has the correct features with the
  #  correct structure and is internally consistent.  Stops if the DNAPL
  #  model is not correct and issues a message indicating what is wrong.
  #
  # the code will stop if the check fails and the problem will be displayed
  #  in the error message
  # --------------------------------------------------------------------- #
  #
  # - ensure that spill.to$layer is integer
  mode(DNAPLmodel@spill.to$layer) <- "integer"
  #
  # - ensure that the columns of mdmax are in the correct order
  DNAPLmodel@mdmax <- DNAPLmodel@mdmax[, DNAPLmodel@domains, drop = FALSE]
  #
  DNAPLmodel.check(DNAPLmodel)
  # --------------------------------------------------------------------- #


  # check timings ----
  #
  # --------------------------------------------------------------------- #
  # Checks whether the MF timings encompass the DNAPL model timings
  #
  # If the DNAPL model period (start.t - end.t) starts before of finishes
  #  after the MF model period, then the start.t or end.t is adjusted as
  #  appropriate and a warning is given.  If the DNAPL model period is
  #  entirely outside the MF model period, then an error is given.
  # --------------------------------------------------------------------- #
  #
  # - ensure start.t, end.t and dt are single values
  if(length(start.t) != 1L){
    start.t <- start.t[1L]
    warning("DNST: start.t should be length-1 numeric")
  }
  if(length(end.t) != 1L){
    end.t <- end.t[1L]
    warning("DNST: end.t should be length-1 numeric")
  }
  if(length(dt) != 1L){
    dt <- dt[1L]
    warning("DNST: dt should be length-1 numeric")
  }
  #
  # - check correct sorting
  if(start.t >= end.t) stop("DNST: start.t must be less than end.t")
  if(dt > end.t - start.t) warning("DNST: dt longer than model period")
  #
  # - check period match with MODFLOW model
  if(!missing(mfdata)){
    MFtrg <- range(mftstime(mfdata, TRUE))

    if(start.t < MFtrg[1L]){
      if(end.t <= MFtrg[1L]){
        stop("DNST: DNAPL model period entirely before MODFLOW period")
      }else{
        start.t <- MFtrg[1L]
        warning("DNST: start.t brought forward to start of MODFLOW model")
      }
    }

    if(end.t > MFtrg[2L]){
      if(start.t >= MFtrg[2L]){
        stop("DNST: DNAPL model period entirely after MODFLOW period")
      }else{
        end.t <- MFtrg[2L]
        warning("DNST: end.t brought back to end of MODFLOW model")
      }
    }
  }
  # --------------------------------------------------------------------- #


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
  uhist$t <- vapply(uhist$year, td, numeric(1L), d = 1, m = 7)
  fsphist$t <- vapply(fsphist$year, td, numeric(1L), d = 1, m = 7)
  J.mult <- switch(uhist.J_units,
                   "kg/day" = 1,
                   "kg/year" = 1/365.25,
                   "t/year" = 1000/365.25,
                   stop("invalid units for uhist.J_units: accepts \"kg/day\", \"kg/year\" or \"t/year\""))
  uhist$cons <- uhist$cons*J.mult
  #
  # - make into functions of time
  uhist <- approxfun(uhist[, c("t", "cons")], yleft = 0, yright = 0,
                     ties = "ordered")
  fsphist <- approxfun(fsphist[, c("t", "f")], rule = 2L, ties = "ordered")
  #
  # - apply losses from exports, safe disposal and volatilisation (not
  #    infiltrating)
  Jin <- function(t) uhist(t)*fi*fw*fsphist(t)
  # --------------------------------------------------------------------- #


  # determine transient Darcy velocity through each layer of the source term ----
  #
  # --------------------------------------------------------------------- #
  # If qh is given explicitly and mfdata not given, then check whether the
  #  format is correct and convert to a consistent format (layer-by-layer
  #  function of time) if necessary.  Otherwise read from mfdata (NetCDF
  #  data set of MODFLOW results).  This requires the code to match up the
  #  DNAPL model layers to the MODFLOW model layers (may not be a simple
  #  corresondence).  This is done by comparing the elevations of the DNAPL
  #  model layer midpoints (hence needing z0) with the layer elevation
  #  divides in the MODFLOW model.  The Darcy velocity is then calculated
  #  using the external qh function.  The matrix (layer, time step) result
  #  is then converted to a list of functions using approxfun.
  #
  # qh: list with length NLAY of functions of time
  # --------------------------------------------------------------------- #
  #
  # - x,y co-ordinates and C,R reference
  if(is.null(x)){
    x <- y <- NA_real_
    if(!missing(mfdata)) stop({
      "DNST: no x or y values given, but mfdata is given"
    })
  }else if(!missing(mfdata)){
    xy <- xy.coords(x, y)
    x <- xy$x; y <- xy$y
    mfC <- cellref.loc(x, gccs(mfdata, TRUE))
    mfR <- cellref.loc(y, grcs(mfdata, TRUE), TRUE)
  }
  #
  # - determine which MODFLOW layer each DNAPL model layer is in (from
  #    midpoint of height)
  #  -- determine z0 if necessary
  if(z0 == "base"){
    z0 <- if(!missing(mfdata)){
      var.get.nc(mfdata, "elev",
                 c(mfC, mfR, dim.inq.nc(mfdata, "NLAY")$length + 1L),
                 c(1L, 1L, 1L))
    }else NA_real_
  }
  #
  if(!missing(mfdata)){
    #  -- elevations of the midpoints of each DNAPL model layer
    DNl.mp <- z0 + cumsum(DNAPLmodel@hL) - DNAPLmodel@hL/2
    #
    #  -- layer references
    mfL <- cellref.loc(DNl.mp,
                       c(rev(var.get.nc(mfdata, "elev",
                                        c(mfC, mfR, NA), c(1L, 1L, NA)))),
                       TRUE)
    #
  }
  # - calculate Darcy Velocity (or check it if given explicitly)
  qh <- if(missing(mfdata)){
    if(is.null(qh)) stop("no MODFLOW results given to mfdata and no value for qh given")

    if(is.function(qh)) qh <- list(qh)
    if(is.numeric(qh)) qh <- lapply(qh, function(q){function(t) q})

    if(length(qh) != DNAPLmodel@NLAY) stop("number of layers for qh do not match the number of layers in the DNAPL model (DNAPLmodel@NLAY)")
    qh
  }else{
    tmp <- qh(x, y, mfL, mfdata)

    if(is.null(dim(tmp))){
      # the qh result is a vector because there is only one MODFLOW layer
      #  used
      rep(list(approxfun(mftstime(mfdata, TRUE)[-1L], tmp, rule = 2:1)),
          DNAPLmodel@NLAY)
    }else{
      # the qh result is a matrix with a row for each unique MODFLOW layer
      #  used
      lapply(mfL, function(l){
        # uses the rownames attached to the matrix within qh
        approxfun(mftstime(mfdata, TRUE)[-1L], tmp[paste0("L", l),],
                  rule = 2L)
      })
    }
  }
  # --------------------------------------------------------------------- #


  # pre-allocate arrays for the results ----
  #
  # --------------------------------------------------------------------- #
  # tvals: time at end of each time step of DNAPL model
  # nts: number of time steps
  # nlay: number of layers in DNAPL model
  # M: 3D array with dim1 representing layers, dim2 representing time steps
  #  and dim3 (named) representing domains
  # Mtop, Mbot: mass escaped from top and bottom of model, by time step
  #  (converted to cumulative mass later)
  # dts: time step lengths
  # overspill: function to evaluate excess mass in domains
  # cascade: an expression for moving mass between domains and layers
  #  according to values calculated by the mdredist processes
  #
  # Also checks that none of the qh functions will make NA for any tval,
  #  if requested with check.qh
  # --------------------------------------------------------------------- #
  #
  # - determine time step end times
  tvals <- seq(start.t, end.t, dt)
  if(tvals[length(tvals)] < end.t - 1) tvals <- c(tvals, end.t)
  nts <- length(tvals)
  #
  nlay <- DNAPLmodel@NLAY
  #
  # - mass array
  M <- array(0,
             c(NLAY = nlay, NTS = nts, Ndom = length(DNAPLmodel@domains)),
             list(LAY = NULL, TS = NULL, dom = DNAPLmodel@domains))
  #
  # - mass lost from top and bottom
  Mtop <- Mbot <- double(nts)
  #
  # - time step lengths
  dts <- diff(tvals)
  #
  # - input flux vector, centrally weighted in time
  Jinv <- vapply((tvals[-1L] + tvals[-nts])/2, Jin, double(1L))
  #
  # - excess mass in domains
  overspill <- function(TS, mdmax){
    dif <- M[, TS,] - mdmax
    ifelse(dif > 0, dif, 0)
  }
  #
  # - cascade mass where domains have become overfull
  #  -- executed twice for each time step, so preassigned as an expression
  cascade <- expression(while(any((os <- {
    overspill(TS + 1L, DNAPLmodel@mdmax)
  }) > 0)){
    # os is a matrix with colnames corresponding to DNAPLmodel@domains and
    #  nlay rows
    lapply(DNAPLmodel@domains, function(d){
      osd <- os[, d]

      # if any overspill; otherwise skip
      if(any(osd > 0)){
        # take excess mass from overspilling domain
        M[, TS + 1L, d] <<- M[, TS + 1L, d] - osd

        # give excess mass to domain to which spill is directed, or cascade
        #  down into next cell
        use.d2 <- "domain2" %in% names(DNAPLmodel@spill.to)
        with(DNAPLmodel@spill.to[d,], {
          # convert from factor
          domain <- as.character(domain)
          if(use.d2) domain2 <- as.character(domain2)

          # domain2 needed if the following ever evaluates TRUE
          domains <- if(use.d2){
            ifelse(M[, TS, domain] >= DNAPLmodel@mdmax[, domain],
                   if(exists("domain2")) domain2 else NA_real_,
                   domain)
          }else rep(domain, nlay)

          if(identical(layer, 1L)){
            # case that mass is lost to layer below
            Mbot[TS + 1L] <<- Mbot[TS + 1L] + osd[nlay]
            lind <- seq(2L, by = 1L, length.out = nlay - 1L)
          }else if(identical(layer, -1L)){
            # (rare) case that mass is lost to layer above
            Mtop[TS + 1L] <<- Mtop[TS + 1L] + osd[1L]
            lind <- seq(1L, by = 1L, length.out = nlay - 1L)
          }else{
            # case that mass is lost within the layer
            lind <- seq_len(nlay)
          }

          for(i in seq_len(length(lind))){
            M[lind[i], TS + 1L, domains[lind[i]]] <<-
              M[lind[i], TS + 1L, domains[lind[i]]] +
              osd[lind[i] - layer]
          }
        })

      }
      NULL
    })
  })
  #
  # - this function environment
  fe <- environment()
  #
  # - check qh if requested
  if(check.qh){
    lapply(1:nlay, function(l){
      if(!is.function(qh[[l]])) stop({
        paste0("DNST: layer ", l, " qh is not a function")
      })
      if(any(is.na(vapply(tvals, qh[[l]], numeric(1L))))) stop({
        paste0("DNST: layer ", l, " qh returning 1 or more NA values within model time range")
      })
    })
  }
  #
  # - check spill input
  if(any(is.na(Jinv) | Jinv < -1e-10)) stop({
    "DNST: some invalid spill inputs (from uhist and fsphist); either negative or NA"
  })
  if(!any(Jinv > 1e-10)) warning("DNST: all spill input 0 (from uhist and fsphist)")
  # --------------------------------------------------------------------- #


  # run model ----
  #
  # --------------------------------------------------------------------- #
  # During each time step:
  # 1. bring mass forward from end of last time step
  # 2. add spill mass to top of model (from Jinv)
  # 3. cascade any excess mass
  # 4. calculate mass transfers from DNAPLmodel@mdredist (in a random
  #     order)
  # 5. cascade any excess mass
  # --------------------------------------------------------------------- #
  #
  for(TS in 1:(nts - 1L)){
    # bring mass forward
    M[, TS + 1L,] <- M[, TS,]

    # add new mass to top
    M[1L, TS + 1L, DNAPLmodel@domain1] <-
      M[1L, TS + 1L, DNAPLmodel@domain1] + Jinv[TS]*dts[TS]

    # cascade mass from overfull domains
    eval(cascade)

    # redistribute mass between domains (including plume) according to the
    #  mdredist functions
    # - the order in which each process is applied is randomised, so that
    #    there is no biassing in the long run
    lmdr <- length(DNAPLmodel@mdredist)
    lapply(DNAPLmodel@mdredist[sample(1:lmdr, lmdr)], function(process){
      # find which parameters are needed
      args <- names(formals(process@flux))
      needpars <- args[!args %in% c("fromM", "toM", "LAY", "time", "env")]

      # determine transfer
      transfer <- double(nlay)
      for(LAY in 1:nlay) transfer[LAY] <- {
        do.call(process@flux, c(list(fromM = M[LAY, TS + 1L, process@from],
                                     toM = M[LAY, TS + 1L, process@to],
                                     LAY = LAY,
                                     time = mean(tvals[TS + 0:1]),
                                     env = fe),
                                DNAPLmodel@params[needpars]))*dts[TS]
      }

      # ensure mass of any domain doesn't decrease below 0
      transfer <- ifelse(transfer < M[, TS + 1L, process@from],
                         transfer, M[, TS + 1L, process@from])
      transfer <- ifelse(-transfer < M[, TS + 1L, process@to],
                         transfer, M[, TS + 1L, process@to])

      # execute transfer
      M[, TS + 1L, process@from] <<-
        M[, TS + 1L, process@from] - transfer
      M[, TS + 1L, process@to] <<-
        M[, TS + 1L, process@to] + transfer
      NULL
    })

    # cascade again from overspills due to redistributed mass
    eval(cascade)
  }
  #
  # - convert bottom and top losses to cumulative losses
  Mbot <- cumsum(Mbot)
  Mtop <- cumsum(Mtop)
  # --------------------------------------------------------------------- #


  # post-process results ----
  #
  # --------------------------------------------------------------------- #
  # 1. Calculate the effluent flux to the plume domain, which serves as a
  #     source term to a contaminant transport model.
  # 2. Calculate mass imbalance by time step, as a check that the model has
  #     worked properly.
  # --------------------------------------------------------------------- #
  #
  # - calculate effluent flux to plume, if there is a plume domain
  J <- if("plume" %in% DNAPLmodel@domains){
    t(apply(adrop(M[,, "plume", drop = FALSE],
                  c(FALSE, FALSE, TRUE)), 1L,
            function(r) c(0, diff(r)/dts)))
  }else matrix(numeric(0L))
  #
  # - determine mass imbalance (should be 0 apart from machine imprecision)
  imbalance <- cumsum(c(0, Jinv*dts)) - apply(M, 2L, sum) - Mbot - Mtop
  # --------------------------------------------------------------------- #


  # save results ----
  #
  # --------------------------------------------------------------------- #
  #
  # --------------------------------------------------------------------- #
  #
  results <- DNAPLSourceTerm(Jspill = Jinv,
                             Jeffluent = J,
                             M = M,
                             Mbot = Mbot,
                             Mtop = Mtop,
                             time = tvals,
                             DNAPLmodel = DNAPLmodel,
                             imbalance = imbalance,
                             description = as.character(description),
                             xy = as.numeric(c(x, y)),
                             z0 = as.numeric(z0))
  #
  # - save to file
  saveRDS(results, result.file)
  #
  # - return invisibly
  invisible(results)
  # --------------------------------------------------------------------- #
}

#' DNAPL Source Term result
#'
#' Format for the results from \code{\link{DNST}}
#'
#' @slot Jspill numeric \code{[NTS]};
#' DNAPL spill flux by each time step
#' @slot Jeffluent numeric matrix \code{[NLAY, NTS]};
#' source term effluent flux to the plume, by layer (matrix row) and time
#'  step (matrix column); if there was no plume domain, this is an empty
#'  matrix
#' @slot M numeric array \code{[NLAY, NTS, Ndom]};
#' source zone mass by layer (dim1), time step (dim2) and mass domain
#'  (dim3)
#' @slot Mbot,Mtop numeric \code{[NTS]};
#' cumulative mass escaped from the bottom and top of the DNAPL model
#' @slot time numeric \code{[NTS]};
#' time values at the end of each time step
#' @slot DNAPLmodel DNAPLmodel;
#' The input DNAPL model (see \code{\link{DNAPLmodel}})
#' @slot imbalance numeric \code{[NTS]};
#' mass imbalance by time step; mass should be completely accounted for, so
#'  any \code{imbalance} asides from machine imprecision constitutes an
#'  error
#' @slot description character string;
#' user-supplied information specific to the model, for future reference
#' @slot xy numeric \code{[2]};
#' location of the spill
#' @slot z0 numeric \code{[1]};
#' elevation that the bottom of the source term model refers to
#'
#' @return
#' DNAPLSourceTerm object; the results from \code{\link{DNST}}
#'
#' @export
#'
DNAPLSourceTerm <- setClass("DNAPLSourceTerm",
                            slots = c(Jspill = "numeric",
                                      Jeffluent = "matrix",
                                      M = "array",
                                      Mbot = "numeric",
                                      Mtop = "numeric",
                                      time = "numeric",
                                      DNAPLmodel = "DNAPLmodel",
                                      imbalance = "numeric",
                                      description = "character",
                                      xy = "numeric",
                                      z0 = "numeric"))
