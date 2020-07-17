openfn.gds.safe <- function(filename) {
    tryCatch({
        openfn.gds(filename)
    }, error=function(e) {
        stop("\nFile has already been opened, please close and try again:\n ", filename,
             "\nIf in doubt, just close all GDS files using this command:\n showfile.gds(closeall=T)")
    })
}

