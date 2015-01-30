\name{heat.clust}
\alias{heat.clust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare dendrograms for gplots' heatmap.2.
}
\description{
This function will prepare dendrograms for the heatmap.2. The type of scaling can be adjusted and is performed before dendrogram calculations (as opposed to native heatmap.2), reordering can be turned on/off and distance and clustering functions can be customized.
}
\usage{
heat.clust(x, scaledim = "column", zlim=c(-3,3), reorder=TRUE,
              scalefun = scale,
              distfun = function(x) as.dist(1-cor(t(x))),
              hclustfun = function(x) hclust(x, method="complete")
          )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric matrix of the values to be plotted.}
  \item{scaledim}{character indicating if the values should be centered and scaled in either the row direction or the column direction, or none.}
  \item{zlim}{Reassign the extremes within the scaled data: zlim=c(-3,3).}
  \item{reorder}{Switch on/off dendrogram reordering.}
  \item{scalefun}{Function to do the data scaling.}
  \item{distfun}{Function used to compute the distance (dissimilarity) between both rows and columns.}
  \item{hclustfun}{Function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms. }
  
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with the scaled data ($data), row ($Rowv) and column ($Colv) dendrograms.}
\references{
%% ~put references to the literature/web site here ~
http://stackoverflow.com/questions/17924828/differences-in-heatmap-clustering-defaults-in-r-heatplot-versus-heatmap-2
}
\author{
Original function by Thomas W. Leja.
Extended by Jan Stanstrup, \email{stanstrup@gmail.com}.
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\dontrun{
#scalefun <- BioMark:::scalefun("auto")
scalefun <- scale

distfun <- function(x) as.dist(1-cor(t(x)))
#distfun <- function(x) dist(x,method="canberra")
hclustfun <- function(x) hclust(x, method="complete")

z <- heat.clust(x,scaledim="column",zlim=c(-3,3),reorder=T,
                distfun=distfun, 
                hclustfun=hclustfun,
                scalefun=scalefun)


library(gplots)
p2 <- heatmap.2(z$data, trace="none", scale="none", 
                Rowv=z$Rowv, 
                Colv=z$Colv,
                reorder=F,
                symbreak=T,
                col=bluered(50),
                na.color="black"
                )
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
