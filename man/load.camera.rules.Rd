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
\value{matrix defining adduct rules.}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
\examples{
\dontrun{
## Excel files with the rules can be found at:
system.file("extdata","CAMERA_rules_pos.xlsx", package="chemhelper")
system.file("extdata","CAMERA_rules_neg.xlsx", package="chemhelper")
system.file("extdata","CAMERA_rules_ei.xlsx", package="chemhelper")
}

}
