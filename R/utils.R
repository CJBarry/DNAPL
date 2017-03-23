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

#' excess mass in domains
#'
#' @param TS time step number
#' @param mdmax from DNAPLmodel
#'
#' @return
#' matrix with colnames corresponding to DNAPLmodel@domains and nlay rows
#'
overspill <- function(TS, mdmax){
  dif <- M[, TS,] - mdmax
  ifelse(dif > 0, dif, 0)
}