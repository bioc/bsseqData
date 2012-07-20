\name{BS.cancer.ex.fit}
\alias{BS.cancer.ex.fit}
\docType{data}
\title{
  WGBS for colon cancer, chr 21 and 22, including smoothed methylation values
}
\description{
  Whole-genome bisulfite sequencing data (WGBS) for colon cancer on
  chromosome 21 and 22.  3 patients were sequenced and the data contains
  matched colon cancer and normal colon.  This dataset includes smoothed
  methylation values.
}
\usage{data(BS.cancer.ex.fit)}
\format{
  The data is stored as an object of class \code{"BSseq"}.
}
\details{
  The file \file{scripts/create_BS.cancer.R} (see example for location)
  is a script that generates all data objects in this package from the
  raw alignment output, contained in the directory \file{umtab} (see
  example for location).  The raw alignment output is the output from
  the BSmooth alignment suite, using an old (legacy) format.

  This dataset is exactly like \code{\link{BS.cancer.ex}} except it also
  contains smoothed methylation values.
}
\seealso{
  \code{\link{BS.cancer.ex}}, \code{\link{BS.cancer.ex.tstat}}
  (t-statistics for this dataset) and \code{\link{keepLoci.ex}} for
  related datasets as well as the
  \code{"\link[bsseq:BSseq-class]{BSseq}"} class and the
  \code{\link[bsseq:BSmooth]{BSmooth}} function.  Also see the
  vignette(s) in the \pkg{bsseq} package.
}
\references{
Hansen, K. D. et al. (2011) \emph{Increased methylation variation in
  epigenetic domains across cancer types}.  Nature Genetics 43, 768-775.
}
\examples{
data(BS.cancer.ex.fit)
BS.cancer.ex.fit

script <- system.file("scripts", "create_BS.cancer.R",
                      package = "bsseqData")
script
readLines(script)
system.file("umtab", package = "bsseqData")
}
\keyword{datasets}
