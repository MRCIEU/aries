#' Saves an ARIES release following QC and normalization
#'
#' @param release.path Release directory path.
#' @param samplesheet Samplesheet for all samples.
#' @param qc.objects List of QC objects lists (one for each feature set, i.e. 'common', '450' and 'epic').
#' @param norm.objects List of normalised QC object lists corresponding to `qc.objects`.
#' @param cell.counts List of cell count estimate matrices (one for each cell count reference).
#' @param qc.summaries List of QC report summaries.
#' @param norm.summaries List of normalization report summaries.
#' @param qc.report.path Directory path of the QC reports.
#' @param norm.report.path Directory path of the normalization reports.
#' @param betas.path Directory path of the normalized betas GDS files, one for each element in `qc.objects`.
#' @param detp.path Directory path of the detection p-value GDS files, one for each element in `qc.objects`.
#' 
#' @export
aries.save.release <- function(release.path, samplesheet,
                               qc.objects, norm.objects,
                               cell.counts,
                               qc.summaries, norm.summaries,
                               qc.report.path, norm.report.path,
                               betas.path, detp.path) {

    ## check inputs
    stopifnot(is.list(qc.objects))
    stopifnot(is.list(norm.objects))
    stopifnot(is.list(qc.summaries))
    stopifnot(is.list(norm.summaries))
    stopifnot(identical(names(qc.objects), names(norm.objects)))
    stopifnot(identical(names(qc.objects), names(qc.summaries)))
    stopifnot(all(names(qc.summaries) %in% names(norm.summaries)))
    for (name in names(qc.objects))
        stopifnot(all(names(qc.objects[[name]]) %in% samplesheet$Sample_Name))
    for (name in names(norm.objects))
        stopifnot(identical(names(norm.objects[[name]]),
                            names(qc.objects[[name]])))
    for (reference in names(cell.counts))
        stopifnot(all(rownames(cell.counts[[reference]])
                      %in% samplesheet$Sample_Name))
    stopifnot(identical(list.files(betas.path, pattern=".gds$"), list.files(detp.path)))
    for (filename in list.files(betas.path, pattern=".gds$", full.names=T)) {
        name <- sub("\\.[^.]+$", "", basename(filename))
        stopifnot(identical(retrieve.gds.dims(filename)[[2]],
                            names(qc.objects[[name]])))
    }
    for (filename in list.files(detp.path, full.names=T)) {
        name <- sub("\\.[^.]+$", "", basename(filename))
        stopifnot(identical(retrieve.gds.dims(filename)[[2]],
                            names(qc.objects[[name]])))
    }

    ## create new release directory structure
    if (!dir.exists(release.path))
        dir.create(release.path, recursive=T)    
    dir.create(file.path(release.path, "samplesheet"))
    dir.create(file.path(release.path, "qc_objects"))
    dir.create(file.path(release.path, "control_matrix"))
    dir.create(file.path(release.path, "norm_objects"))
    dir.create(file.path(release.path, "derived"))
    dir.create(file.path(release.path, "derived", "reports"))
    dir.create(file.path(release.path, "derived", "cellcounts"))
    dir.create(file.path(release.path, "betas"))
    dir.create(file.path(release.path, "detection_p_values"))

    not.file.exists <- function(filename) {
        if (file.exists(filename)) {
            cat("Already exists, skipping: ", filename, "\n")
            FALSE
        }
        else
            TRUE
    }
    
    ## samplesheets
    cat(date(), "Saving samplesheets\n")
    save.path <- file.path(release.path, "samplesheet")
    filename <- file.path(save.path, "samplesheet.csv")
    if (not.file.exists(filename))
        write.csv(
            samplesheet, 
            file=filename,
            row.names=F)
    
    samplesheets <- sapply(names(qc.objects), function(name) {
        idx <- match(names(qc.objects[[name]]), samplesheet$Sample_Name)
        samplesheet[idx,]
    }, simplify=F)
    
    for (name in names(samplesheets)) {
        filename <- file.path(save.path, paste0("samplesheet-", name, ".csv"))
        if (not.file.exists(filename))
            write.csv(samplesheets[[name]],
                      file=filename,
                      row.names=F)
    }

    ## save qc objects
    cat(date(), "Saving QC objects\n")
    save.path <- file.path(release.path, "qc_objects")
    for (name in names(qc.objects)) {
        filename <- file.path(save.path, paste(name, "rds", sep="."))
        if (not.file.exists(filename))
            saveRDS(qc.objects[[name]],
                    file=filename)
    }

    ## save clean qc reports
    cat(date(), "Saving final QC reports\n")
    save.path <- file.path(release.path, "derived", "reports", "qc")
    R.utils::copyDirectory(qc.report.path, save.path)
    for (name in names(qc.summaries)) {
        filename <- file.path(save.path, paste(name, "rds", sep="."))
        if (not.file.exists(filename))
            saveRDS(qc.summaries[[name]], filename)
    }
    
    ## save normalization objects
    cat(date(), "Saving normalization objects\n")
    save.path <- file.path(release.path, "norm_objects")
    for (name in names(norm.objects)) {
        filename <- file.path(save.path, paste(name, "rds", sep="."))
        if (not.file.exists(filename))
            saveRDS(norm.objects[[name]],
                    file=filename)
    }

    ## save control matrices
    cat(date(), "Saving control matrices\n")
    save.path <- file.path(release.path, "control_matrix")
    for (name in names(qc.objects)) {
        cm <- sapply(qc.objects[[name]], function(obj) obj$controls)
        cm <- t(cm)
        filename <- file.path(save.path, paste(name, "txt", sep="."))
        if (not.file.exists(filename))
            write.table(cm, file=filename,
                        sep="\t", col.names=T, row.names=T)
    }

    ## save cell count estimates
    cat(date(), "Saving cell count estimates\n")
    save.path <- file.path(release.path, "derived", "cellcounts")
    for (reference in names(cell.counts)) {
        filename <- file.path(save.path, paste(reference, "txt", sep="."))
        if (not.file.exists(filename))
            write.table(cell.counts[[reference]],
                        file=filename,
                        sep="\t",quote=F,row.names=F,col.names=T)
    }

    ## save normalization reports
    cat(date(), "Saving normalization reports\n")
    save.path <- file.path(release.path, "derived", "reports", "normalization")
    R.utils::copyDirectory(norm.report.path, save.path)
    for (name in names(norm.summaries)) {
        filename <- file.path(save.path, paste(name, "rds", sep="."))
        if (not.file.exists(filename))
            saveRDS(norm.summaries[[name]], filename)
    }

    ## save normalized methylation data
    cat(date(), "Saving normalized methylation matrices\n")
    save.path <- file.path(release.path, "betas")
    for (name in names(samplesheets)) {
        gds.filename <- file.path(save.path, paste0(name, ".gds"))
        if (not.file.exists(gds.filename)) 
            R.utils::copyFile(file.path(betas.path, basename(gds.filename)), save.path)
    } 
    save.path <- file.path(release.path, "betas", "no-missing")
    if (file.exists(file.path(betas.path, "no-missing"))) {
        dir.create(save.path)
        for (name in names(samplesheets)) {
            gds.filename <- file.path(save.path, paste0(name, ".gds"))
            if (not.file.exists(gds.filename)) 
                R.utils::copyFile(file.path(betas.path, "no-missing", basename(gds.filename)), save.path)
        }
    }
    
    ## save detection p-values
    cat(date(), "Saving detection p-values\n")
    save.path <- file.path(release.path, "detection_p_values")
    for (name in names(qc.objects)) {
        gds.filename <- file.path(save.path, paste(name,"gds",sep="."))
        if (not.file.exists(gds.filename)) 
            R.utils::copyFile(file.path(detp.path, basename(gds.filename)), save.path)
    }

    release.path
}


#' Create a modified release after removing a set of individuals
#'
#' @param release.path Directory path of master release.
#' @param new.path Directory path of modified release.
#' @param remove.ids ALNQLET identifiers of individuals
#' who should be omitted from the new release.
#' 
#' @export
aries.copy.release <- function(release.path, new.path, remove.ids) {
    dir.create(new.path, recursive=T)
    
    directories <- unique(dirname(list.files(release.path, recursive=T)))
    directories <- directories[!grepl("reports", directories)]
    directories <- directories[!grepl("qc_objects", directories)]
    directories <- directories[!grepl("norm_objects", directories)]
    for (dir in directories)
        dir.create(file.path(new.path, dir), recursive=T)
    dir.create(reports.path <- file.path(new.path, "derived", "reports"))
    
    samplesheet <- aries.samples(release.path)
    is.included <- !samplesheet$alnqlet %in% as.character(remove.ids)
    samplesheet <- samplesheet[is.included,]
    write.csv(
        samplesheet, 
        file=file.path(new.path, "samplesheet", "samplesheet.csv"),
        row.names=F)
    
    info <- aries.info(release.path)
    for (name in names(info)) {
        sample.names <- with(
            info[[name]]$samples,
            Sample_Name[!alnqlet %in% as.character(remove.ids)])
        aries <- aries.select(
            release.path,
            featureset=name,
            sample.names=sample.names)

        gds.filename <- paste(name, "gds", sep=".")
        cat(date(), "Saving", name, "betas\n") 
        copy.gds.matrix(
            file.path(release.path, "betas", gds.filename),
            file.path(new.path, "betas", gds.filename),
            sample.names=sample.names)
        if (file.exists(file.path(release.path, "betas", "no-missing")))
            copy.gds.matrix(
                file.path(release.path, "betas", "no-missing", gds.filename),
                file.path(new.path, "betas", "no-missing", gds.filename),
                sample.names=sample.names)
        
        cat(date(), "Saving", name, "detection p-values\n") 
        copy.gds.matrix(
            file.path(release.path, "detection_p_values", gds.filename),
            file.path(new.path, "detection_p_values", gds.filename),
            sample.names=sample.names)

        cat(date(), "Saving", name, "samplesheet\n") 
        write.csv(
            aries$samples, 
            file=file.path(new.path, "samplesheet", paste0("samplesheet-", name, ".csv")),
            row.names=F)

        cat(date(), "Saving", name, "control matrix\n") 
        write.table(
            aries$control.matrix,
            file=file.path(new.path, "control_matrix", paste0(name,".txt")),
            sep="\t", col.names=T, row.names=T)
        
        if (name == "common") {
            for (reference in names(aries$cell.counts)) {
                cat(date(), "Saving", reference, "cell count estimates\n") 
                cc <- aries$cell.counts[[reference]]
                cc <- na.omit(cc)
                write.table(
                    data.frame(IID=rownames(cc), cc),
                    file=file.path(new.path, "derived", "cellcounts", paste0(reference, ".txt")),
                    sep="\t", quote=F,row.names=F,col.names=T)
            }
        }

        cat(date(), "Generating QC report\n")
        qc.objects <- readRDS(file.path(release.path, "qc_objects", paste0(name,".rds")))
        qc.parameters <- meffil.qc.parameters(
          beadnum.samples.threshold             = 0.1,  ## remove samples with >10% probes with bead count < 3
          detectionp.samples.threshold          = 0.1,  ## remove samples with >10% probes with detection p > 0.01 
          detectionp.cpgs.threshold             = 0.1,  ## remove cpgs with ...
          beadnum.cpgs.threshold                = 0.1,  ## remove cpgs with ...
          sex.outlier.sd                        = 5,    ## flag pred sex score outliers (>5 SDs)
          snp.concordance.threshold             = 0.95, ## flag snp probes with below 95% concordance genotypes
          sample.genotype.concordance.threshold = 0.8)  ## flag samples with <80% snp concordance with genotypes
        report.filename <- file.path(reports.path, "qc", paste0("qc-report-",name, ".html"))
        in.subset <- names(qc.objects) %in% aries$samples$Sample_Name
        qc.summary <- meffil.qc.summary(qc.objects[in.subset], parameters=qc.parameters)
        meffil.qc.report(qc.summary, output.file=report.filename)
        
        cat(date(), "Generating normalization report\n")
        norm.objects <- readRDS(file.path(release.path, "norm_objects", paste0(name,".rds")))
        meth.filename <- file.path(release.path, "betas", paste0(name,".gds"))
        report.filename <- file.path(
          reports.path,
          paste0("norm-report-", name, ".html"))
        in.subset <- names(norm.objects) %in% aries$samples$Sample_Name
        pcs <- meffil.methylation.pcs(
          meth.filename,
          probe.range=20000,
          samples=aries$samples$Sample_Name,
          autosomal=T)
        norm.parameters<-meffil.normalization.parameters(
          norm.objects,
          variables=intersect(c("Slide","BCD_plate","time_point"), colnames(aries$samples)),
          control.pcs=1:10,
          batch.pcs=1:10)
        norm.summary <- meffil.normalization.summary(norm.objects[in.subset], pcs,parameters=norm.parameters)
        meffil.normalization.report(norm.summary,output.file=report.filename)
    }
    
    invisible(new.path)
}



