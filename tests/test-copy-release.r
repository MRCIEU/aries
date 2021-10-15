library(aries)

library(alspac)
alspac::setDataDir(alspac.dir)

exclusions <- alspac:::readExclusions()
for (group in names(exclusions)) {
    if (grepl("child", group)) 
        exclusions[[group]] <- paste0(exclusions[[group]], c("A","B"))
    if (grepl("mother", group))
        exclusions[[group]] <- paste0(exclusions[[group]], "M")
    if (grepl("partner", group))
        exclusions[[group]] <- paste0(exclusions[[group]], "F")
}
exclusions <- unlist(exclusions)

aries.copy.release(aries.dir, new.dir, exclusions) ## 3-4 hours

for (fs in c("common", "450", "epic")) {
    old.15up <- aries.select(aries.dir, featureset=fs, time.point="15up")
    new.15up <- aries.select(new.dir, featureset=fs, time.point="15up")
    old.15up$meth <- aries.methylation(old.15up)
    new.15up$meth <- aries.methylation(new.15up)

    cat("featureset = ", fs, "\n") 
    cat("identical meth = ", 
        identical(
            old.15up$meth[,colnames(new.15up$meth)],
            new.15up$meth),
        "\nidentical control.matrix = ",
        identical(
            old.15up$control.matrix[rownames(new.15up$control.matrix),],
            new.15up$control.matrix),
        "\n")
    for (type in names(old.15up$cell.counts)) {
        cat("identical ", type, " = ",
            identical(
                old.15up$cell.counts[[type]][colnames(new.15up$meth),],
                new.15up$cell.counts[[type]]),
            "\n")
    }
    cat("identical samples = ",
        identical(
            old.15up$samples[rownames(new.15up$samples),],
            new.15up$samples),
        "\nidentical probe.names = ",
        identical(old.15up$probe.names, new.15up$probe.names),
        "\n")
}


