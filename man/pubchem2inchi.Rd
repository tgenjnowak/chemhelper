\name{pubchem2inchi}
\alias{pubchem2inchi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert PubChem CID to InChIs
}
\description{
Convert PubChem CID to InChIs. 
Calls obabel locally. Hence obabel need to be available in path.
}
\usage{
pubchem2inchi(cid,skip=NULL,silent=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cid}{Numeric vector of PubChem CIDs.}
  \item{skip}{Vector of indexes to skip in cid. Output index will still match input indexes.}
  \item{silent}{Suppress printing of progress info.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A character vector of InChIs.}
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
