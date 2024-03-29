#' aries.info
#'
#' Provide comprehensive sample and probe information for all of ARIES
#'
#' @param path Base directory of ARIES.
#' @return A list containing ARIES information for each Illumina BeadChip format,
#' run \code{\link{aries.feature.sets}()} to get a list of valid identifiers.
#' Each of these is in turn a list of items describing the samples
#' with methylation profiles in that format:
#' \itemize{
#' \item samples - Data frame providing sample information, one sample per row.
#' \item probe.names - Illumina probe identifiers.
#' \item control.matrix - Matrix providing control probe values (rows=samples, columns=control probes).
#' \item cell.counts - List of cell count estimates for various references. Each is a matrix (rows=samples, columns=cell types). 
#' }
#' @examples \dontrun{
#' aries <- aries.info(aries.dir)
#' ## See CD4+ T cell estimate for the third sample
#' ## with an Illumina EPIC profile.
#' aries$epic$cell.counts[["blood-gse35069"]][3,"CD4T"]
#' }
#' @export
aries.info <- function(path) {
    samples <- aries.samples(path)

    gds.filenames <- list.files(file.path(path, "betas"), pattern=".gds$", full.names=T)
    names(gds.filenames) <- sub(".gds$", "", basename(gds.filenames))
    sapply(names(gds.filenames), function(featureset) {
        gds.filename <- gds.filenames[[featureset]]
        ids <- retrieve.gds.dims(gds.filename)
        probe.names <- ids[[1]]
        sample.names <- ids[[2]]        
        samples <- samples[match(sample.names, samples$Sample_Name),]
        
        control.matrix <- read.table(file.path(path, "control_matrix", paste(featureset, "txt", sep=".")),
                                     row.names=1, header=T, sep="\t",check.names=F)
        control.matrix <- as.matrix(control.matrix)

        cell.count.filenames <- list.files(file.path(path, "derived", "cellcounts"), ".txt$", full.names=T)
        counts <- sapply(cell.count.filenames, read.table, sep="\t", row.names=1, header=T, simplify=F)
        names(counts) <- sub(".txt", "", sapply(names(counts), basename))
        counts <- lapply(counts, as.matrix)
        counts <- lapply(counts, function(counts) counts[match(sample.names, rownames(counts)),])
        names(counts) <- gsub(" ", "-", names(counts))

        list(samples=samples,
             probe.names=probe.names,
             control.matrix=control.matrix,
             cell.counts=counts)
    }, simplify=F)
}




