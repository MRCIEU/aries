#' The steps below are needed to regenerate
#' the documentation
#' included with the package.

#' install.packages("devtools")
#' devtools::install_github("klutometis/roxygen")
library(devtools)
library(roxygen2)

document("aries")

system("R CMD Rd2pdf aries")
system("mkdir aries/docs")
system("mv aries.pdf aries/docs")

system("R CMD INSTALL aries")
reload(inst("aries"))


