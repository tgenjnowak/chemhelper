\name{group_correct}
\alias{group_correct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch-correct dataset feature-wise.
}
\description{
Corrects dataset column-wise (feature-wise) such that groups of rows (samples) have the same mean. For each column the mean of each row group (class) is calculated and adjusted to the global row mean.
This is intended to correct for batch effects in LC-MS based metabolomics studies.
}
\usage{
group_correct(dataset,f,zero.is.na=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{The matrix to correct.}
  \item{f}{Vector defining the row (sample) grouping.}
  \item{zero.is.na}{Treat zero values as missing data for the calculation of means.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Matrix or data.frame depending in the input object.}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
