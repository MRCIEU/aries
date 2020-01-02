#' aries.feature.sets
#'
#' List of Illumina BeadChip formats available in ARIES
#'
#' @param path Base directory of ARIES.
#' @return Vector listing available format identifiers.
#' @export
aries.feature.sets <- function(path) {
    betas.path <- file.path(path, "betas")
    if (!dir.exists(betas.path))
        stop("Invalid path for ARIES data: ", path) 
    gds.filenames <- list.files(betas.path, pattern=".gds$", full.names=T)
    sub(".gds$", "", basename(gds.filenames))
}
