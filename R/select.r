#' aries.select
#'
#' Provide comprehensive sample and probe information for a subset of ARIES
#'
#' @param path Base directory of ARIES.
#' @param featureset Illumina BeadChip format (Default: 'common').
#' Run \code{\link{aries.feature.sets}()} to get a list of valid identifiers.
#' @param time.point Optional. Restrict to profiles for samples collected at specific
#' time points (Default: NULL).
#' Run \code{\link{aries.time.points}()} to get a list.
#' @param sample.names Optional. Restrict to profiles for a set of samples (Default: NULL).
#' Run \code{\link{aries.samples}()$Sample_Name} to get a list of sample identifiers.
#' @param probe.names Optional. Restrict to a set of Illumina probe identifiers (Default: NULL).
#' @param replicates Include replicate methylation profiles (Default: FALSE).
#' Ignored if \code{is.null(sample.names)==FALSE}.
#' @return A list containing ARIES information for the ARIES subset:
#' \itemize{
#' \item samples - Data frame providing sample information, one sample per row.
#' \item probe.names - Illumina probe identifiers.
#' \item control.matrix - Matrix providing control probe values (rows=samples, columns=control probes).
#' \item cell.counts - List of cell count estimates for various references. Each is a matrix (rows=samples, columns=cell types). 
#' }
#' @examples \dontrun{
#' ## Define subsets of ARIES restricted to samples collected around age 15
#' ## with Illumina 450K and EPIC BeadChip profiles.
#' aries.15.450 <- aries.select(aries.dir, featureset="450", time.point="15up")
#' aries.15.epic <- aries.select(aries.dir, featureset="epic", time.point="15up")
#' }
#' @export
aries.select <- function(path, featureset="common", time.point=NULL, sample.names=NULL, probe.names=NULL, replicates=F) {
    aries.info <- aries.info(path=path)

    if (!featureset %in% names(aries.info))
        stop("'", featureset, "' is not a valid feature set name. ",
             "Try one of the following: ", paste(names(aries.info), collapse=", "))
    
    aries <- aries.info[[featureset]]

    if (!is.null(sample.names))
        replicates <- T

    if (is.null(sample.names)) 
        is.selected.sample <- rep(T, nrow(aries$samples))
    else {
        is.selected.sample <- aries$samples$Sample_Name %in% sample.names
        if (!all(sample.names %in% aries$samples$Sample_Name)) {
            if (all(!is.selected.sample))
                stop("None of the sample names appear in this dataset.")
            warning("Some of the sample names do not appear in this dataset.")
        }
    }

    if (!is.null(time.point)) {
        if (!all(time.point %in% aries$samples$time_point))
            stop("Time point(s) ", paste(setdiff(time.point, aries$samples$time_point), collapse="/"),
                 " not available for feature set '", featureset, "'.")
        is.selected.sample <- is.selected.sample & aries$samples$time_point %in% time.point
    }
    
    if (!replicates)
        is.selected.sample <- is.selected.sample & !aries$samples$duplicate.rm

    if (is.null(probe.names))
        is.selected.probe <- rep(T,length(aries$probe.names))
    else {
        is.selected.probe <- aries$probe.names %in% probe.names
        if (!all(probe.names %in% aries$probe.names)) {
            if (all(!is.selected.probe))
                stop("None of the probe names can be found in this dataset.")
            warning("Some of the probe names cannot be found in this dataset.")
        }
    }

    list(path=path,
         featureset=featureset,
         samples=aries$samples[which(is.selected.sample),],
         probe.names=aries$probe.names[which(is.selected.probe)],
         control.matrix=aries$control.matrix[which(is.selected.sample),],
         cell.counts=lapply(aries$cell.counts, function(counts) counts[which(is.selected.sample),]))
}


