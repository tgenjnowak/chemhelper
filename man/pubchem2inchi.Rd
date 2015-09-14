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
\value{A character vector of InChIs.}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
