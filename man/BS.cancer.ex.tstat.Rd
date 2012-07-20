\name{BS.cancer.ex.tstat}
\alias{BS.cancer.ex.tstat}
\docType{data}
\title{
  T-statistics for WGBS data for colon cancer, chr 21 and 22
}
\description{
  T-statistics produced by the \code{BSmooth.tstat} function, run on the
  \code{BS.cancer.ex.fit} object subsetted by \code{keepLoci.ex}.
}
\usage{data(BS.cancer.ex.tstat)}
\format{
  The data is stored as an object of class \code{"BSseqTstat"}.
}
\details{
  See below for the script creating this object.
}
\seealso{
    \code{\link{BS.cancer.ex.fit}} (data used to produce the t-statistics)
  and \code{\link{keepLoci.ex}} (used for subsetting) as well as the
  \code{"\link[bsseq:BSseqTstat-class]{BSseqTstat}"} class and
  \code{\link[bsseq:BSmooth.tstat]{BSmooth.tstat}}.  Also see the
  vignette(s) in the \pkg{bsseq} package.
}
\references{
Hansen, K. D. et al. (2011) \emph{Increased methylation variation in
  epigenetic domains across cancer types}.  Nature Genetics 43, 768-775.
}
\examples{
data(BS.cancer.ex.tstat)
BS.cancer.ex.tstat
## This script shows how the object was created
script <- system.file("scripts", "create_BS.cancer.R",
                      package = "bsseqData")
script
readLines(script)
}
\keyword{datasets}
