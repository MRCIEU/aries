openfn.gds.safe <- function(filename) {
    tryCatch({
        openfn.gds(filename)
    }, error=function(e) {
        stop("\nFile has already been opened, please close and try again:\n ", filename,
             "\nIf in doubt, just close all GDS files using this command:\n showfile.gds(closeall=T)")
    })
}

## load a matrix from a GDS file
retrieve.gds.matrix <- function(gds.filename, sites, samples) {
    stopifnot(file.exists(gds.filename))
    gds.file <- openfn.gds(gds.filename)
    on.exit(closefn.gds(gds.file))        

    all.sites <- read.gdsn(index.gdsn(gds.file, "row.names"))
    if (is.null(sites)) sites <- all.sites
    else stopifnot(all(sites %in% all.sites))

    all.samples <- read.gdsn(index.gdsn(gds.file, "col.names"))
    if (is.null(samples)) samples <- all.samples
    else stopifnot(all(samples %in% all.samples))

    matrix.node <- index.gdsn(gds.file, "matrix")
    mat <- readex.gdsn(
        matrix.node,
        sel=list(all.sites %in% sites, all.samples %in% samples), 
        simplify="none")
    rownames(mat) <- all.sites[which(all.sites %in% sites)]
    colnames(mat) <- all.samples[which(all.samples %in% samples)]
    mat[sites,samples,drop=F]
}

## load the dimension names from the GDS matrix
retrieve.gds.dims <- function(gds.filename) {
    stopifnot(file.exists(gds.filename))
    gds.file <- openfn.gds(gds.filename)
    on.exit(closefn.gds(gds.file))
    list(read.gdsn(index.gdsn(gds.file, "row.names")),
         read.gdsn(index.gdsn(gds.file, "col.names")))                  
}
