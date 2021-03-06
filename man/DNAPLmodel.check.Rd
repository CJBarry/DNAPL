% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DNST_model.R
\name{DNAPLmodel.check}
\alias{DNAPLmodel.check}
\title{Check that a DNAPL model is set up correctly}
\usage{
DNAPLmodel.check(d, stop.if.mistakes = TRUE)
}
\arguments{
\item{d}{DNAPLmodel object to be checked, as given by \code{\link{DNAPLmodel}}}

\item{stop.if.mistakes}{logical \code{[1]};
should R stop executing if there are any mistakes (default \code{TRUE})}
}
\value{
A character vector explaining all mistakes in the DNAPL model
 (\code{d}).  If there are any mistakes, the session will stop if
 \code{stop.if.mistakes == TRUE} and the mistakes will be printed as
 error messages.
}
\description{
The set up of a DNAPL distribution and dissolution model is fairly
 detailed (see \code{\link{DNAPLmodel}}).  This function provides an
 easy way to check that a DNAPL model is set up correctly.  If it isn't,
 a vector of messages is returned that inform the user of all the things
 that need to be corrected.  A DNAPL model that passes this check should
 work with \code{\link{DNST}}.
}
\note{
Some more obscure errors will not be picked up by this function.  In the
 \code{spill.to} slot, it is possible to define a cycle of overspilling
 (e.g. if domain A spills to domain B and domain B spills to domain A,
 all in the same layer), which could result in an infinite loop during
 model calculation.  This error would not be detected by
 \code{DNAPLmodel.check}.
}

