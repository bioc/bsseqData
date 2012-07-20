library(bsseq)
library(tools) ## for resaveRdaFiles below

sampleNames <- list.files("../umtab")
dirs <- file.path("../umtab", sampleNames, "umtab")

umData <- read.umtab(dirs, sampleNames)

BS.cancer.ex <- umData$BSdata
seqlevels(BS.cancer.ex) <- paste0("chr", seqlevels(BS.cancer.ex))
pData(BS.cancer.ex)$Type <- rep(c("cancer", "normal"), times = c(3,3))
pData(BS.cancer.ex)$Pair <- rep(c("pair1", "pair2", "pair3"), 2)
validObject(BS.cancer.ex)
BS.cancer.ex <- chrSelectBSseq(BS.cancer.ex,
                               seqnames = c("chr21", "chr22"), order = TRUE)
save(BS.cancer.ex, file = "BS.cancer.ex.rda")
resaveRdaFiles("BS.cancer.ex.rda") ## additional compression

## Assumes you have 6 cores available, and enough RAM
BS.cancer.ex.fit <- BSmooth(BS.cancer.ex, ns = 70, h = 1000,
                            mc.cores = 6, verbose = TRUE)
save(BS.cancer.ex.fit, file = "BS.cancer.ex.fit.rda")
resaveRdaFiles("BS.cancer.ex.fit.rda") ## additional compression

BS.cov <- getCoverage(BS.cancer.ex.fit)
keepLoci.ex <- which(rowSums(BS.cov[, c("C1", "C2", "C3")] >= 2) >= 2 &
                     rowSums(BS.cov[, c("N1", "N2", "N3")] >= 2) >= 2)
save(keepLoci.ex, file = "keepLoci.ex.rda")
resaveRdaFiles("keepLoci.ex.rda")

BS.cancer.ex.tstat <- BSmooth.tstat(BS.cancer.ex.fit[keepLoci.ex,], 
                                    group1 = c("C1", "C2", "C3"),
                                    group2 = c("N1", "N2", "N3"), 
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
save(BS.cancer.ex.tstat, file = "BS.cancer.ex.tstat.rda")
resaveRdaFiles("BS.cancer.ex.tstat.rda")

