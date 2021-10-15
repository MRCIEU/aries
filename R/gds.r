## load a matrix from a GDS file
retrieve.gds.matrix <- function(gds.filename, sites=NULL, samples=NULL) {
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

## create an new empty matrix in a GDS file
create.gds.matrix <- function(gds.file, storage, n.row, n.col, row.names=NULL, col.names=NULL, mat=NULL) {
    if (!is.null(mat)) {
        stopifnot(is.matrix(mat))
        stopifnot(nrow(mat) == n.row)
        stopifnot(ncol(mat) == n.col)
    }
    
    if (!is.null(row.names)) {
        stopifnot(n.row == length(row.names))
        add.gdsn(gds.file, "row.names", row.names)
    }
    if (!is.null(col.names)) {
        stopifnot(n.col == length(col.names))
        add.gdsn(gds.file, "col.names", col.names)
    }
    add.gdsn(gds.file,
             name="matrix",
             storage=storage,
             valdim=c(n.row, n.col),
             val=mat,
             replace=TRUE)
}

## save the columns in 'mat' to the GDS matrix
## starting at the given column
save.gds.columns <- function(node, mat, start) {
    if (!is.matrix(mat))
        mat <- matrix(mat, ncol=1)
    write.gdsn(node, mat, start=c(1,start), count=c(nrow(mat), ncol(mat)))
}

## save the rows in 'mat' to the GDS matrix
## starting at the given column
save.gds.rows <- function(node, mat, start) {
    if (!is.matrix(mat))
        mat <- matrix(mat, nrow=1)
    write.gdsn(node, mat, start=c(start,1), count=c(nrow(mat), ncol(mat)))
}

## copy one GDS matrix file to another file
## restricted to the samples identified in sample.names
copy.gds.matrix <- function(in.filename, out.filename, sample.names=NULL, block.size=1000) {
    stopifnot(file.exists(in.filename))
    gds.in <- openfn.gds(in.filename)
    on.exit(closefn.gds(gds.in))
    mat.in <- index.gdsn(gds.in, "matrix")
    col.names <- read.gdsn(index.gdsn(gds.in, "col.names"))
    row.names <- read.gdsn(index.gdsn(gds.in, "row.names"))
    ##row.names <- row.names[1:2500]

    if (is.null(sample.names))
        sample.names <- col.names
    else {
        if (!all(sample.names %in% col.names))
            stop("Some sample identifiers not in the file ", in.filename) 
    }
    col.idx <- match(sample.names, col.names)
    
    gds.out <- createfn.gds(out.filename)
    on.exit(closefn.gds(gds.out), add=T)
    mat.out <- create.gds.matrix(
        gds.out,
        "float64",
        n.row=length(row.names),
        n.col=length(sample.names),
        row.names=row.names,
        col.names=sample.names)

    start.idx <- 1
    while (start.idx < length(row.names)) {
        max.size <- length(row.names)-start.idx+1
        mat <- read.gdsn(
            mat.in,
            start=c(start.idx, 1),
            count=c(ifelse(max.size > block.size, block.size, max.size), -1),
            simplify="none")
        save.gds.rows(mat.out, mat[,col.idx], start=start.idx)
        start.idx <- start.idx + nrow(mat)
    }
    out.filename
}

## gdsfmt::showfile.gds(closeall=T)
