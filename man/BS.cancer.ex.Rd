\name{BS.cancer.ex}
\alias{BS.cancer.ex}
\docType{data}
\title{
  WGBS for colon cancer, chr 21 and 22
}
\description{
  Whole-genome bisulfite sequencing data (WGBS) for colon cancer on
  chromosome 21 and 22.  3 patients were sequenced and the data contains
  matched colon cancer and normal colon.
}
\usage{data(BS.cancer.ex)}
\format{
  The data is stored as an object of class \code{"BSseq"}.
}
\details{
  The file \file{scripts/create_BS.cancer.R} (see example for location)
  is a script that generates all data objects in this package from the
  raw alignment output, contained in the directory \file{umtab} (see
  example for location).  The raw alignment output is the output from
  the BSmooth alignment suite, using an old (legacy) format.

  This dataset \code{\link{BS.cancer.ex.fit}} is the same basic data,
  but it also contains smoothed methylation values.
}
\seealso{
  \code{\link{BS.cancer.ex.fit}}, \code{\link{BS.cancer.ex.tstat}}
  (t-statistics for this dataset) and \code{\link{keepLoci.ex}} for
  related datasets and the \code{"\link[bsseq:BSseq-class]{BSseq}"}
  class.  Also see the vignette(s) in the \pkg{bsseq} package.
}
\references{
Hansen, K. D. et al. (2011) \emph{Increased methylation variation in
  epigenetic domains across cancer types}.  Nature Genetics 43, 768-775.
}
\examples{
data(BS.cancer.ex)
BS.cancer.ex

script <- system.file("scripts", "create_BS.cancer.R",
                      package = "bsseqData")
script
readLines(script)
system.file("umtab", package = "bsseqData")
}
\keyword{datasets}
