#' aries.samples
#'
#' Obtain sample information for all of ARIES.
#'
#' @param path Base directory of ARIES.
#' @return Data frame listing all samples.
#' @export
aries.samples <- function(path) {
    samplesheet.filename <- file.path(path, "samplesheet", "samplesheet.csv")
    if (!file.exists(samplesheet.filename))
        stop("Invalid path for ARIES data: ", path)
    read.csv(samplesheet.filename, stringsAsFactors=F)
}
