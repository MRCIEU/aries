#' aries.time.points
#'
#' List of ARIES time-point identifiers
#'
#' @param path Base directory of ARIES.
#' @return A vector of valid time-point identifiers.
#' @export
aries.time.points <- function(path) {
    samples <- aries.samples(path)
    unique(samples$time_point)
}
