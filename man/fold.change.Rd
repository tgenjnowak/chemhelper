\name{fold.change}
\alias{fold.change}
\alias{ratio.w.invert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fold change
}
\description{
Calculates fold change for all columns/variables in a matrix based on grouping vector.
}
\usage{
fold.change(MAT,f,aggr_FUN=mean,combi_FUN="/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MAT}{Matrix or data.frame where each column is considered a variable and each row are the samples.}
  \item{f}{Sample grouping vector. Fold change will be calculated between all group conbinations.}
  \item{aggr_FUN}{function to aggregate replicates. Usually \code{\link{mean}} or \code{\link{median}}.}
  \item{combi_FUN}{Function to calculate the change between aggregated group values. "/" is default (just the ratio). If \code{\link{ratio.w.invert}} is specified the fold change will be positive if the the first group in the comparision notation is higher and negative if the second group in the comparision notation is higher.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A matrix giving the fold changes for each factor combination.}
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
