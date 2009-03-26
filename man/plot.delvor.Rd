\name{plot.delvor}
\alias{plot.delvor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the Voronoi diagram and Delaunay traingulation}
\description{
    This function returns a plot of the Voronoi diagram and Delaunay traingulation. 
}
\usage{
\method{plot}{delvor}(x, add = FALSE, wlines = c("both", "del", "vor"),
	wpoints = TRUE, number = FALSE, col = NULL, 
	xlim = NULL, ylim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{An object of class \code{"delvor"} as constructed by the function delvor. }
\item{add}{Logical, if TRUE add to a current plot.} 
\item{wlines}{"Which lines?". I.e. should the Delaunay triangulation be plotted (wlines='del'), should the Voronoi diagram be plotted (wlines='vor'), or should both be plotted (wlines='both', the default)?}
 \item{wpoints}{Logical, indicates if sample points should be plotted.}
\item{number}{Logical, defaulting to FALSE; if TRUE then the points plotted will be labelled with their index numbers.}  
  \item{col}{The colour numbers for plotting the data points, Delaunay triangulation, Voronoi diagram, and the point numbers, in that order; defaults to c(1,1,1,1). If fewer than four numbers are given, they are recycled. (If more than four numbers are given, the redundant ones are ignored.)}
\item{xlim}{The limits on the x-axis.}  
\item{ylim}{The limits on the y-axis.}  
\item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).}
}
\seealso{\code{\link{delvor}}.}
\examples{
# Random sample in the unit square
x<-matrix(runif(100),nc=2)
# Delaunay triangulation and Voronoi diagram
delvor.obj<-delvor(x)
# Plot Voronoi diagram and Delaunay triangulation 
plot(delvor.obj)
}
\keyword{nonparametric}