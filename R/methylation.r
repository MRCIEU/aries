#' aries.methylation
#'
#' Load DNA methylation profiles for a specified ARIES subset
#'
#' @param selection Output of \code{\link{aries.select}()}.
#' @param start Index of the first CpG in the profile (Default: 1).
#' Indexing is with respect to \code{selection$probe.names}.
#' @param count Number of CpG sites to load (Default: -1).
#' A value of \code{-1} indicates all available.
#' @param probe.names Names of CpG sites to load (Default: NULL).
#' If not NULL, then `start/count` are ignored and
#' methylation levels are loaded for the specified CpG sites.
#' @param no.missing If TRUE, then load DNA methylation from a
#' file without undetected probe values missing (Default: FALSE).
#' @return A matrix filled with methylation estimates from
#' 0 to 1 (rows=CpG sites, columns=samples) corresponding to
#' the requested CpG sites.
#' @examples \dontrun{
#' ## Load the first 1000 CpG sites for all blood samples collected around age 15
#' aries <- aries.select(aries.dir, featureset="common", time.point="15up")
#' aries$meth <- aries.methylation(aries, start=1, count=1000)
#' aries$meth <- aries.methylation(aries, probe.names=c("cg19642007", "cg12038298", "cg09884480", "cg1291275"))
#' }
#' @export
aries.methylation <- function(selection, start=1, count=-1, probe.names=NULL, no.missing=FALSE) {
    if (!is.list(selection)
        && 
        !all(c("samples","probe.names","control.matrix","cell.counts") %in% names(selection)))
        stop("'selection' was not created by aries.select()")
    
    sample.names <- selection$samples$Sample_Name

    if (is.null(probe.names)) {    
        probe.names <- selection$probe.names
        if (start > length(probe.names))
            stop("'start' is too big")
        if (start + count - 1 > length(probe.names)) 
            count <- length(probe.names)-start+1
        if (count > 0)
            probe.names <- probe.names[start:(start+count-1)]
        else
            probe.names <- probe.names[start:length(probe.names)]
    } else {
        common.names <- intersect(probe.names, selection$probe.names)
        if (length(common.names) == 0)
            stop("none of 'probe.names' matches probes for this 'selection'")
        if (length(common.names) < length(probe.names))
            probe.names <- common.names     
    }
    
    gds.filename <- file.path(
        selection$path,
        ifelse(!no.missing, "betas", file.path("betas","no-missing")),
        paste(selection$featureset, "gds", sep="."))

    retrieve.gds.matrix(gds.filename, probe.names, sample.names)
}
