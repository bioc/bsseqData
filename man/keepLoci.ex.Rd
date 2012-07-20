\name{keepLoci.ex}
\alias{keepLoci.ex}
\docType{data}
\title{
  Which methylation loci were included in an analysis of BS.cancer.ex.
}
\description{
  This object describes which methylation loci were kept, when
  t-statistics were generated from \code{BS.cancer.fit.ex} using the
  function \code{BSmooth.tstat}.
}
\usage{data(keepLoci.ex)}
\format{
  A vector of indices into \code{BS.cancer.fit.ex}.
}
\details{
  See below how this object was created and used.
}
\seealso{
  \code{\link{BS.cancer.ex.fit}} (this is the data the subsetting index
  works on) and \code{\link{BS.cancer.ex.tstat}} and the
  \code{\link[bsseq:BSmooth.tstat]{BSmooth.tstat}} function .  Also see
  the vignette(s) in the \pkg{bsseq} package.
}
\examples{
data(keepLoci.ex)
## This script shows how the object was created
script <- system.file("scripts", "create_BS.cancer.R",
                      package = "bsseqData")
script
readLines(script)
}
\keyword{datasets}
