\name{diff.mat2letter.not}
\alias{diff.mat2letter.not}
\alias{edges2diff.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Letter notation for overlapping confidence intervals.
}
\description{
Letter notation for overlapping confidence intervals familiar in multiple comparisons testing.
edges2diff.mat turns upper and lower bounds of confidence intervals into a logical matrix indicating which entries overlap. diff.mat2letter.not creates the letter notation.
}
\usage{
edges2diff.mat(lower,upper)
diff.mat2letter.not(g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower}{Lower bound of confidence intervals.}
  \item{upper}{Upper bound of confidence intervals.}
  \item{g}{Turns a logical vector of overlaps into letter notation.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Character vector with the letter notation.}
\references{
%% ~put references to the literature/web site here ~
http://menugget.blogspot.it/2014/05/automated-determination-of-distribution.html
http://stackoverflow.com/questions/27293770/compact-letter-display-from-logical-matrix/27296818
}
\author{
Original code by "Marc in the box" and modified by Jan Stanstrup, \email{stanstrup@gmail.com}.
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

test_data <- data.frame(        mean = c(1.48, 1.59, 1.81,1.94),
                            CI_lower = c(1.29,1.38,1.54, 1.62),
                            CI_upper = c(1.56,1.84, 2.3, 2.59)
                        )


g <- edges2diff.mat(test_data$CI_lower, test_data$CI_upper)
letter_notation <- diff.mat2letter.not(g)
letter_notation

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
