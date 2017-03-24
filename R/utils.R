# DNAPL package - utilities

#' Expand vectors
#'
#' Ensure that a vector is the correct length, recycling length-1 inputs if
#'  necessary
#'
#' @param x
#' vector \code{[1]} or \code{[N]}
#' @param N
#' target number of elements
#'
#' @return
#' vector \code{[N]}
#'
#' @export
#'
#' @examples
#' expand.vec(1, 5L)
#' expand.vec(1:5, 5L)
#'
#' # works with lists
#' expand.vec(list(1), 3L)
#'
#' # only recycles length-1 vectors
#' \dontrun{expand.vec(1:3, 5L)}
#' \dontrun{expand.vec(1:3, 6L)}
#'
expand.vec <- function(x, N){
  stopifnot(length(x) %in% c(1L, N))
  if(is.vector(x)){
    vec <- rep(x[1L], N)
    vec[] <- x
    vec
  }else{
    rep(list(x), N)
  }
}

#' Convert national usage trend to site usage rate through time
#'
#' @param nat.uhist
#' data.frame or one of \code{"TCEusage"}, \code{"PCEusage"},
#'  \code{"TCAusage"} or \code{"TeCMusage"};
#' data frame must have one column called \code{cons}, which will be read
#'  as the usage rate
#' @param pu
#' numeric \code{[1]};
#' the peak usage rate at the site ever - the national usage rate is
#'  normalised to give it a maximum value of \code{pu}; because \code{cons}
#'  is normalised, \code{pu} needn't be in the same units - it should be in
#'  the desired units of output and kg/day is suggested
#'
#' @return
#' data.frame with adapted \code{cons} column
#'
#' @export
#'
#' @examples
#' UK.to.site("TCEusage", 10)
#'
UK.to.site <- function(nat.uhist, pu){
  e <- environment()
  nat.uhist <- switch(class(nat.uhist),
                      character = {
                        data(list = nat.uhist, package = "DNAPL", envir = e)
                        get(nat.uhist, e)
                      },
                      data.frame = nat.uhist,
                      stop({
                        "UK.to.site: invalid input for the national usage history, nat.uhist"
                      }))

  if(!"cons" %in% names(nat.uhist))
    stop("UK.to.site: nat.uhist must have a 'cons' column, representing the nation-wide usage rate")

  nat.uhist$cons <- nat.uhist$cons/max(nat.uhist$cons)*pu
  nat.uhist
}
