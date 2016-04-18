\name{mzmine2xcmsSet}
\alias{mzmine2xcmsSet}
\title{
Convert an mzMine peaklist to an xcmsSet.
}
\description{
Convert an mzMine peaklist to an xcmsSet. This allows to use R based tools such as CAMERA.
}
\usage{
mzmine2xcmsSet(peaklist_mzmine, filepath = NULL, get.scantimes = FALSE)
}
\arguments{
  \item{peaklist_mzmine}{The peaktable as exported from mzMine}
  \item{filepath}{The file path to the raw files (to find the files).}
  \item{get.scantimes}{Logical indicating if the scantimes should be added to the xcmsSet. It is slow and not always needed depending on the down-stream analysis.}
}
\value{An \code{\link{xcmsSet-class}}.}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}

