\name{db.comp.assign}
\alias{db.comp.assign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotation of retention time / m/z data
}
\description{
Match a list of retention times and mz values with a database containing known rt and m/z values of compounds.
}
\usage{
db.comp.assign(mz,rt,comp_name_db,mz_db,rt_db,mzabs=0.01,ppm=15,ret_tol=30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mz}{Vector of mass values to assign compounds to.}
  \item{rt}{Vector of retention times to assign compounds to.}
  \item{comp_name_db}{Character vector of compound names.}
  \item{mz_db}{Numeric vector of m/z values corresponding to known compounds.}
  \item{rt_db}{Numeric vector of retention times corresponding to known compounds.}
  \item{ppm,mzabs}{Mass tolerance defined both by ppm and absolute deviation. Both limits must be observed a database match.}
  \item{ret_tol}{Retention time tolerance for a database match.}
}

\value{A character vector of same length as mz with compound names where a match was found. Compound names are annotated with "non-unique hit" if multiple input values match the database entry.}

\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
