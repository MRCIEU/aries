ewas <- c(## sex differences (ARIES stats for PACE study)
          "sex-aries-0-unadj", 
          "sex-aries-7-unadj",
          "sex-aries-15-unadj",
          ## Smoking (Joehanes et al. CCG 2016)
          "smoking-joehanes-ccg-2016",
          ## Prenatal smoking (ARIES stats in Joubert et al. EHP 2012)
          "pte-joubert-ehp-2012-aries",
          ## Birthweight (Simpkin et al. HMG 2015)
          "birthweight-simpkin-hmg-2015",
          ## Gestational age (Simpkin et al. HMG 2015)
          "gestationalage-simpkin-hmg-2015",
          ## Adult BMI (Mendelson et al. PloSMed 2017)
          "bmi-mendelson-plosmed-2017")

read.mendelson <- function(filename) {
    require(readxl)
    ret <- read_xlsx(filename, skip=2)
    ret <- as.data.frame(ret)
    cols <- colnames(ret)
    names(cols) <- cols
    cols[c("Probename","β...9","p-value...11")] <- c("cpg","effect","p.value")
    colnames(ret) <- cols
    ret
}

ewas <- sapply(ewas, function(name) {
    filename <- file.path("ewas", name)
    if (grepl("mendelson", name))
        ret <- read.mendelson(paste0(filename, ".xlsx"))
    else if (file.exists(paste0(filename, ".txt")))
        ret <- read.table(paste0(filename, ".txt"), sep="\t",header=T)
    else
        ret <- read.csv(paste0(filename, ".csv"), comment.char="#")
    if ("t.stat" %in% colnames(ret) && !"effect" %in% colnames(ret))
        ret$effect <- ret$t.stat    
    n <- min(100, nrow(ret))
    ret[order(ret$p.value)[1:n],]
}, simplify=F)


