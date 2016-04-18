\name{xsAnnotate_stats}
\alias{xsAnnotate_stats}
\title{
xsAnnotate statistics.
}
\description{
Give some statistics from an xsAnnotate object annotated by CAMERA.
}
\usage{
xsAnnotate_stats(object)
}
\arguments{
  \item{object}{\code{\link{xsAnnotate}} to analyze.}
}
\value{
A matrix giving the following statistics: \cr
'# Feature' \cr
'# Feature groups' \cr
'# Features groups with several features' \cr
'# Features in feature groups with several features' \cr
'# Features annotated as isotopes' \cr
'# Feature with adducts/fragments annotation' \cr
'# Feature groups with fragment/adduct annotation' \cr
'# Annotation groups' \cr
'Mean number of annotation groups per feature group with annotations' \cr
'Median number of annotation groups per feature group with annotations'
}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
\examples{
 library(CAMERA)
 library(chemhelper)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
 an   <- findIsotopes(an)  # optional but recommended.
 an   <- groupCorr(an) # optional but very recommended step
 an   <- findAdducts(an,polarity="positive")
 xsAnnotate_stats(an)
}
