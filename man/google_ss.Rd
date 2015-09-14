\name{google_ss}
\alias{google_ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download data from google spreadsheet
}
\description{
Download data from google spreadsheet.
}
\usage{
google_ss(gid = NA, key = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{gid}{Gid of the spreadsheet. Need to be public}
\item{key}{Key of the specific sheet/tab. Integer starting from 0.}
}

\value{data.frame with spreadsheet data.}
\references{
http://www.r-bloggers.com/r-function-to-read-data-from-google-docs-spreadsheets/

}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
