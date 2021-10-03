## @knitr aries-init
library(aries)
aries.dir <- "/projects/MRC-IEU/research/data/alspac/epigenetic/methylation/aries/dev/release_candidate/data/release" ## "path/to/aries"

## @knitr alspac-init
## devtools::install_github("explodecomputer/alspac")
library(alspac)
alspac.dir <- "~/work/alspac/data" #"path/to/alspac"
alspac::setDataDir(alspac.dir)

## @knitr check-aries.feature.sets
stopifnot(identical(aries.feature.sets(aries.dir), c("450","common","epic")))

## @knitr check-aries.time.points
stopifnot(identical(aries.time.points(aries.dir), c("FOM","antenatal","cord","F7","15up","FOF","c43m","c61m","F9","F24")))

## @knitr apply-aries.select
ds.15up <- aries.select(aries.dir, featureset="common", time.point="15up")
ds.15up.450 <- aries.select(aries.dir, featureset="450", time.point="15up")
ds.15up.epic <- aries.select(aries.dir, featureset="epic", time.point="15up")

## @knitr combine-manually
controls <- intersect(colnames(ds.15up$control.matrix), colnames(ds.15up.epic$control.matrix))
ds.check <- list(cell.counts=sapply(names(ds.15up$cell.counts), function(ref)
                     rbind(ds.15up.epic$cell.counts[[ref]], ds.15up.450$cell.counts[[ref]]), simplify=F),
                 control.matrix=rbind(ds.15up.epic$control.matrix[,controls],
                                      ds.15up.450$control.matrix[,controls]),
                 samples=rbind(ds.15up.epic$samples, ds.15up.450$samples))

## @knitr identical.rows
identical.rows <- function(x,y) {
    common <- intersect(rownames(x), rownames(y))
    length(common) == nrow(x) && identical(x[common,], y[common,])
}

## @knitr check-aries.select-cell.counts
for (ref in names(ds.15up$cell.counts))
    stopifnot(identical.rows(ds.15up$cell.counts[[ref]], ds.check$cell.counts[[ref]]))

## @knitr check-aries.select-samples
stopifnot(identical.rows(ds.15up$samples, ds.check$samples))

## @knitr check-aries.select-control.matrix
idx <- match(rownames(ds.15up$control.matrix),
             rownames(ds.check$control.matrix))
r <- cor(t(scale(ds.15up$control.matrix[,controls])),
         t(scale(ds.check$control.matrix[idx,controls])))
rn <- apply(-r,1,rank)
stopifnot(sum(diag(rn) > 1)/nrow(rn) < 0.05)

## @knitr load-alspac
source("load-alspac.r")

## @knitr load-ewas
source("load-ewas.r")

## @knitr load-alcohol-models
## devtools::install_github("yousefi138/dnamalci")
library(dnamalci)

## @knitr load-age-models
## devtools::install_github("perishky/meffonym")
library(meffonym)

## @knitr add-models
reese <- read.csv("ewas/pte-reese-ehp-2017.csv",stringsAsFactors=F)
meffonym.add.model("pte-reese-ehp-2017",
                   reese$cpg, reese$effect,
                   "reese-ehp-2017 model for prenatal smoking")

meffonym.add.model("bmi-mendelson-plosmed-2017",
                   ewas[["bmi-mendelson-plosmed-2017"]]$cpg,
                   ewas[["bmi-mendelson-plosmed-2017"]]$effect,
                   "mendelson-plosmed-2017 model for bmi")

meffonym.add.model("smoking-joehanes-ccg-2016",
                   ewas[["smoking-joehanes-ccg-2016"]]$cpg,
                   ewas[["smoking-joehanes-ccg-2016"]]$effect,
                   "joehanes-ccg-2016 model for own smoking")

meffonym.add.model("sex-aries-15-unadj",
                   ewas[["sex-aries-15-unadj"]]$cpg,
                   ewas[["sex-aries-15-unadj"]]$effect,
                   "ARIES stats for PACE sex differences study")


## @knitr collect-cpg-sites
alcohol.sites <- names(dnamalci.get.model("dnamalc.144cpg")$coefficients)

model.names <- c(
    "hannum",
    "bmi-mendelson-plosmed-2017",
    "bohlin.1se",
    "smoking-joehanes-ccg-2016",
    "pte-reese-ehp-2017",
    "sex-aries-15-unadj")
model.sites <- sapply(model.names, function(model)
                      names(meffonym.get.model(model)$coefs),
                      simplify=F)

ewas.sites <- lapply(ewas, function(ewas) ewas$cpg)

sites <- unique(c(unlist(ewas.sites),
                  unlist(model.sites),
                  unlist(alcohol.sites)))


## @knitr load-methylation
time.points <- aries.time.points(aries.dir)
aries <- sapply(time.points, function(time.point) {
    ds <- aries.select(aries.dir,
                       time.point=time.point)
    ds$meth <- aries.methylation(ds)
    ds$meth <- ds$meth[which(rownames(ds$meth) %in% sites),]
    ds
}, simplify=F) ## 2 minutes

## @knitr site-coverage
kable(sapply(aries, function(x) {
    c(n=nrow(x$meth),
      sapply(
          c(model.sites, ewas.sites, list(alcohol=alcohol.sites)),
          function(sites) {
              floor(mean(sites %in% rownames(x$meth))*100)
          }))
    }))


## order aries time-point subsets by age for convenient outputs
for (i in 1:length(aries))
    aries[[i]]$median.age <- round(median(aries[[i]]$samples$age, na.rm=T),1)
aries$cord$median.age <- 0
aries <- aries[order(sapply(aries, function(ds) ds$median.age))]


## @knitr test-assocs
test.assoc <- function(aries, var, cpg, ...) {
    covs <- list(...)
    covs <- lapply(covs, as.data.frame)
    data <- data.frame(var=var, do.call(cbind, covs))
    meth <- aries$meth[cpg,]
    fit <- lm(meth ~ ., data=data)
    stats <- coef(summary(fit))
    colnames(stats) <- c("effect","se","t.stat","p.value")
    stats[grep("^var",rownames(stats))[1],]
}
test.assocs <- function(aries, var, cpgs, ...) {
    stats <- t(sapply(cpgs, function(cpg) test.assoc(aries, var, cpg, ...)))
    as.data.frame(stats)
}


## @knitr add-alspac-vars
time.points.children <- c("cord","c43m","c61m","F7","F9","15up","F24")

for (tp in time.points.children) {
    idx <- match(aries[[tp]]$samples$alnqlet, alspac$alnqlet)
    aries[[tp]]$samples$gestationalage <- alspac$gestational.age[idx]
    aries[[tp]]$samples$birthweight <- alspac$birthweight[idx]
    aries[[tp]]$samples$pte <- alspac$smoking.prenatal[idx]
    aries[[tp]]$samples$prenatal.alcohol <- alspac$prenatal.alcohol[idx]
}

idx <- match(aries[["F7"]]$samples$alnqlet, alspac$alnqlet)
aries[["F7"]]$samples$bmi <- alspac$bmi.child.7[idx]

idx <- match(aries[["F9"]]$samples$alnqlet, alspac$alnqlet)
aries[["F9"]]$samples$bmi <- alspac$bmi.child.7[idx]

idx <- match(aries[["15up"]]$samples$alnqlet, alspac$alnqlet)
aries[["15up"]]$samples$smoking <- alspac$smoking.child.17[idx]
aries[["15up"]]$samples$cotinine <- alspac$cotinine.child.17[idx]
aries[["15up"]]$samples$bmi <- alspac$bmi.child.17[idx]

idx <- match(aries[["F24"]]$samples$alnqlet, alspac$alnqlet)
aries[["F24"]]$samples$smoking <- alspac$smoking.child.24[idx]
aries[["F24"]]$samples$cotinine <- alspac$cotinine.child.17[idx]
aries[["F24"]]$samples$bmi <- alspac$bmi.child.24[idx]

idx <- match(aries[["antenatal"]]$samples$aln, alspac$aln)
aries[["antenatal"]]$samples$smoking <- alspac$smoking.prenatal[idx]
aries[["antenatal"]]$samples$alcohol <- alspac$prenatal.alcohol[idx]
aries[["antenatal"]]$samples$cotinine <- alspac$cotinine.prenatal[idx]
aries[["antenatal"]]$samples$bmi <- alspac$bmi.pre[idx]

idx <- match(aries[["FOM"]]$samples$aln, alspac$aln)
aries[["FOM"]]$samples$smoking <- alspac$smoking.middle[idx] == "current"
aries[["FOM"]]$samples$alcohol <- alspac$alcohol.middle[idx]
aries[["FOM"]]$samples$audit <- alspac$alcohol.audit[idx]
aries[["FOM"]]$samples$bmi <- alspac$bmi.middle[idx]

idx <- match(aries[["FOF"]]$samples$aln, alspac$aln)
aries[["FOF"]]$samples$smoking <- alspac$smoking.father[idx] == "current"
aries[["FOF"]]$samples$alcohol <- alspac$alcohol.father[idx]
aries[["FOF"]]$samples$bmi <- alspac$bmi.father[idx]


## @knitr test-effect-correlations

r.effects <- sapply(names(aries), function(time.point) {
    ds <- aries[[time.point]]
    sapply(names(ewas), function(ewasname) {
        varname <- sub("([^-]+).*", "\\1", ewasname)
        if (!varname %in% colnames(ds$samples)
            || length(unique(ds$samples[[varname]])) < 2)
            return(NA)        
        old.stats <- ewas[[ewasname]]
        old.stats <- old.stats[old.stats$cpg %in% rownames(ds$meth),]
        counts <- ds$cell.counts[["blood-gse35069-complete"]]
        new.stats <- test.assocs(
            ds,
            ds$samples[[varname]],
            old.stats$cpg,
            batch=as.factor(ds$samples$plate),
            counts=counts)
        cor(old.stats$effect, new.stats$effect, use="p")
    })
})


## @knitr test-model-correlations
r.models <- sapply(
    names(aries),
    function(time.point) {
        ds <- aries[[time.point]]
        counts <- ds$cell.counts[["blood-gse35069-complete"]]
        covs <- data.frame(batch=as.factor(ds$samples$plate),counts)
        covs <- model.matrix(~., covs)
        meth <- meffonym:::impute.mean(ds$meth)
        fit <- lm.fit(x=covs, t(meth))
        meth <- t(residuals(fit))
        scores <- list(
            age=meffonym.score(meth, "hannum")$score,
            bmi=meffonym.score(meth, "bmi-mendelson-plosmed-2017")$score,
            gestationalage=meffonym.score(meth, "bohlin.1se")$score,   
            cotinine=meffonym.score(meth,"smoking-joehanes-ccg-2016")$score,
            alcohol=dnamalci(meth, "dnamalc.144cpg")$score)
        scores$audit <- scores$alcohol
        sapply(
            names(scores),
            function(varname) {
                if (varname %in% colnames(ds$samples))
                    cor(scores[[varname]], ds$samples[[varname]], use="p")
                else
                    NA
            })
    })


## @knitr test-model-auc

library(pROC)

r.auc <- sapply(names(aries), function(time.point) {
    ds <- aries[[time.point]]
    counts <- ds$cell.counts[["blood-gse35069-complete"]]
    covs <- data.frame(batch=as.factor(ds$samples$plate),counts)
    covs <- model.matrix(~., covs)
    meth <- meffonym:::impute.mean(ds$meth)
    fit <- lm.fit(x=covs, t(meth))
    meth <- t(residuals(fit))
    scores <- list(
        pte=meffonym.score(meth, "pte-reese-ehp-2017")$score,
        sex=meffonym.score(meth, "sex-aries-15-unadj")$score,
        smoking=meffonym.score(meth, "smoking-joehanes-ccg-2016")$score)
    sapply(names(scores), function(varname) {
        if (varname %in% colnames(ds$samples)
            && length(na.omit(unique(ds$samples[[varname]]))) > 1)
            auc(ds$samples[[varname]], scores[[varname]])
        else
            NA
    })
})


## @knitr calculate-dnam-age

aries$all <- aries.select(aries.dir)
aries$all$meth <- aries.methylation(aries$all) ## 2 minutes
aries$all$meth <- aries$all$meth[which(rownames(aries$all$meth) %in% sites),]

cord.idx <- which(aries$all$samples$time_point == "cord")
alspac.idx <- match(aries$all$samples$alnqlet[cord.idx], alspac$alnqlet)
aries$all$samples$age[cord.idx] <- (alspac$gestational.age[alspac.idx]-40)/52

aries$all$samples$hannum <- meffonym.score(aries$all$meth, "hannum")$score
idx <- which(!is.na(aries$all$samples$age))
aries$all$samples$age.accel <- NA
aries$all$samples$age.accel[idx] <- residuals(
    lm(hannum ~ .,
       data=cbind(
           aries$all$samples[idx,c("hannum","age","plate")],
           aries$all$cell.counts[["blood-gse35069-complete"]][idx,])))

## @knitr dnam-age-correlations
time.points <- setdiff(names(aries), "all")
dnam.age.r <- cbind(
    "time point"=time.points,
    "median age"=sapply(aries[time.points],
        function(ds) ds$median.age),
    R=sapply(time.points, function(tp) {
        with(aries$all$samples, {
            idx <- which(time_point == tp)
            cor(age[idx], hannum[idx], use="p")
        })
    }))


## @knitr dnam-aa-by-sex
group <- with(aries$all$samples, {
    group <- rep(NA, length(time_point))
    group[which(time_point == "15up")] <- "age 15"
    group[which(time_point == "F24")] <- "age 24"
    group[which(time_point %in% c("FOM","FOF"))] <- "middle age"
    group
})

dnam.aa.by.sex <- t(sapply(c("age 15","age 24","middle age"), function(gp) {
    idx <- which(group == gp)
    with(aries$all$samples[idx,], {
        if (length(unique(na.omit(sex))) < 2) return(rep(NA,4))
        fit <- lm(age.accel ~ sex)
        stats <- coef(summary(fit))[2,]
        names(stats) <- c("difference","se","t-statistic","p.value")      
        stats
    })
}))






