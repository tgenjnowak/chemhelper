\name{clipboard}
\alias{clipboard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copy to clipboard on linux
}
\description{
Copied to clipboard on linux. This works only on linux with xclip installed.
}
\usage{
clipboard(x, sep="\t", row.names=FALSE, col.names=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{The Object to be copied to the clipboard. Preferably a matrix or data frame. If not, it is attempted to coerce x to a data frame.}
\item{sep}{Seperator tu use for copied content.}
\item{row.names}{Copy row names.}
\item{col.names}{Copy column names.}
}

\value{None}

\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
