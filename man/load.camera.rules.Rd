\name{load.camera.rules}
\alias{load.camera.rules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load adduct rules for CAMERA.
}
\description{
This function loads rules for use with CAMERA adduct annotation.
}
\usage{
load.camera.rules(mode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mode}{The ionization mode. Can be "pos", "neg" or "ei".}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{matrix defining adduct rules.}
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
\dontrun{
## Excel files with the rules can be found at:
system.file("extdata","CAMERA_rules_pos.xlsx", package="chemhelper")
system.file("extdata","CAMERA_rules_neg.xlsx", package="chemhelper")
system.file("extdata","CAMERA_rules_ei.xlsx", package="chemhelper")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
