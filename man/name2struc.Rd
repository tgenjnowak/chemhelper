\name{name2struc}
\alias{name2struc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lookup compound name
}
\description{
Lookup compound name and get matching structures using chemical translation service. A table is presented with all hits including PubChem CIDs and InChIs.
}
\usage{
name2struc(input_names,    input_pubchem=as.numeric(matrix(data=NA,nrow=length(input_names)))     )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_names}{Character vector of compound names to query.}
  \item{input_pubchem}{List of PubChems where they are known. This is used to skip query of input_names when pubchem CIP already know. This should be change to a skip argument.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A character matrix containing columns "org_row" (row number in original input), "input_name", "output_name", "pubchem_CID" and "inchi"      }
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
