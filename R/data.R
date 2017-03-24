# DNAPL package - data documentation

#' @rdname usage
#' @name usage
#'
#' @title UK Chlorinated solvent usage
#'
#' @description
#' Estimated chlorinated solvent UK-wide consumption over time in the
#'  twentieth century.
#'
#' @format Data frames with 15 rows:
#' \describe{
#'   \item{year}{}
#'   \item{cons}{estimated total UK consumption, in kt/a (1000 metric
#'    tonnes per year)}
#' }
#'
#' @details
#' TCE: trichloroethene\cr
#' PCE: tetrachloroethene\cr
#' TCA: 1,1,1-trichloroethane\cr
#' TeCM: tetrachloromethane/ carbon tetrachloride
#'
#' Check \code{attr(TCEusage, "units")} to see units displayed.
#'
#' @references
#' Rivett, M.O. et al., 1990. Organic contamination of the Birmingham aquifer, U.K. Journal of Hydrology, 113(1–4), pp.307–323.
#'
#' @source
#' \url{http://linkinghub.elsevier.com/retrieve/pii/002216949090181V}
#'
"TCEusage"

#' @rdname usage
#'
"PCEusage"

#' @rdname usage
#'
"TCAusage"

#' @rdname usage
#'
"TeCMusage"

#' Chlorinated Hydrocarbon properties
#'
#' Properties of selected chlorinated solvents that are relevant to the
#'  DNAPL source term model, namely solubility and density.
#'
#' @format Data frame with 4 named rows:
#' \describe{
#'   \item{solubility}{solubility in water close to 20 degrees C, in
#'    grams/l}
#'   \item{density}{free phase density close to 20 degrees C, in kg/m3}
#' }
#'
#' @details
#' There is significant variation in the figures for TeCM solubility.
#'  Values are for temperatures between 20 and 25 degrees celcius.
#'
#' @source
#' \url{https://www.ncbi.nlm.nih.gov/pccompound}
#'
"CHCprops"
