dir.create("test-output")

sink("test-output/copy-release.Rout")
source("test-copy-release.r")
sink()
