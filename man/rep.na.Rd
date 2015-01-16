\name{rep.na}
\alias{rep.na}
\alias{rand.to.min}
\alias{half.min}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Replacement of NA values
}
\description{
Column-wise replacement of NA values based on non-NA data.
}
\usage{
rep.na(mat,zero.as.na=TRUE,rep.FUN=half.min)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{Matrix or data.frame where each column will be treated seperately.}
  \item{zero.as.na}{If TRUE zero values are considered NA values.}
  \item{rep.FUN}{Function to calculate replacement data. The function is invoked once per missing \code{\link{value.half.min}} will replace by half the minimum value. \code{\link{rand.to.min}} will replace by a random value between zero and the minimum non-NA value.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A matrix with NA values replaced.}
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
