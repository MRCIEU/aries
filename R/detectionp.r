#' aries.detectionp
#'
#' Load detection p-values for a specified ARIES subset
#'
#' @param selection Output of \code{\link{aries.select}()}.
#' @param start Index of the first CpG in the profile (Default: 1).
#' Indexing is with respect to \code{selection$probe.names}.
#' @param count Number of CpG sites to load (Default: -1).
#' A value of \code{-1} indicates all available.
#' @return A matrix filled with detection p-values (rows=CpG sites, columns=samples).
#' @examples \dontrun{
#' ## Load the first 1000 CpG sites for all blood samples collected around age 15
#' aries <- aries.select(aries.dir, featureset="common", time.point="15up")
#' aries$detp <- aries.detectionp(aries, start=1, count=1000)
#' }
#' @export
aries.detectionp <- function(selection, start=1, count=-1) {
    if (!is.list(selection)
        && 
        !all(c("samples","probe.names","control.matrix","cell.counts") %in% names(selection)))
        stop("'selection' was not created by aries.select()")
    
    sample.names <- selection$samples$Sample_Name
    probe.names <- selection$probe.names
    if (start > length(probe.names) || start + count - 1 > length(probe.names))
        stop("'start/count' indexes past the end of the selected probes (", length(probe.names), ")")
    if (count > 0)
        probe.names <- probe.names[start:(start+count-1)]
    else
        probe.names <- probe.names[start:length(probe.names)]
    
    gds.filename <- file.path(selection$path, "detection_p_values", paste(selection$featureset, "gds", sep="."))

    retrieve.gds.matrix(gds.filename, probe.names, sample.names)
}
