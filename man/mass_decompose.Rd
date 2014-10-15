\name{mass_decompose}
\alias{mass_decompose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mass decomposition
}
\description{
Calculate the elementary compositions from an exact mass. This function is s wrapper for Rdisop's decomposeMass function. The output is a table with diagnostic information for each possible formula.
}
\usage{
mass_decompose(mass,intensities=NULL,ion="neutral",elements=c('C','H','N','O','P','S'),ppm=30, filter.DBE=TRUE,filter.nitrogen=TRUE,simplify=TRUE, minElements="C0", maxElements="C999999") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mass}{A single exact mass (or m/z value).}
  \item{intensities}{A vector of intensities of the the monoisotopic peak and subsequent isotopes. If NULL the mass is decomposed without considering the isotopic pattern.}
  \item{ion}{A character string or a character vector with definitions of the type of ion. Multiple species hypotheses can be specified. Possible values are Neutral, Pos, Neg, [M+H]+, [M-H]-, [M+Na]+, [M+K]+, [M-2H+Na]-, [M-2H+K]-, [M+Cl]-, [M-H+HCOONa]-, [M-H+HCOOH]-.}
  \item{elements}{Character vector of allowed chemical elements, defaults to c('C','H','N','O','P','S').}
  \item{ppm}{Allowed deviation of hypotheses from given mass.}
  \item{filter.DBE}{Only allow positive integer DBE values.}
  \item{filter.nitrogen}{Only allow formulas that adhere to the nitrogen rule.}
  \item{simplify}{Output is always a list if multiple values of "ion" is supplied. If TRUE output is a data.frame if a single ion is supplied.}
  \item{minElements, maxElements}{Molecular formulas, which contain lower and upper boundaries of allowed formula respectively.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list containing a data.frame with the following columns: Formula, Nitrogen rule, DBE, Calc. m/z, ppm, Rdisop score. A single data.frame if simplify=TRUE.}
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
