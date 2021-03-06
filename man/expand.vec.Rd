% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{expand.vec}
\alias{expand.vec}
\title{Expand vectors}
\usage{
expand.vec(x, N)
}
\arguments{
\item{x}{vector \code{[1]} or \code{[N]}}

\item{N}{target number of elements}
}
\value{
vector \code{[N]}
}
\description{
Ensure that a vector is the correct length, recycling length-1 inputs if
 necessary
}
\examples{
expand.vec(1, 5L)
expand.vec(1:5, 5L)

# works with lists
expand.vec(list(1), 3L)

# only recycles length-1 vectors
\dontrun{expand.vec(1:3, 5L)}
\dontrun{expand.vec(1:3, 6L)}

}

