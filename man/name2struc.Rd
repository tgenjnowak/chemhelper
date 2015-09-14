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
name2struc(input_names,
           input_pubchem=as.numeric(matrix(data=NA,nrow=length(input_names)))
           )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_names}{Character vector of compound names to query.}
  \item{input_pubchem}{List of PubChems where they are known. This is used to skip query of input_names when pubchem CIP already know. This should be change to a skip argument.}
}
\value{A character matrix containing columns "org_row" (row number in original input), "input_name", "output_name", "pubchem_CID" and "inchi"      }
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
