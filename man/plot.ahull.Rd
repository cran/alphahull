\name{plot.ahull}
\alias{plot.ahull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the alpha-convex hull}
\description{
   This function returns a plot of the \eqn{\alpha}-convex hull. If desired, it also adds the Delaunay triangulation, Voronoi diagram and \eqn{\alpha}-shape of the sample. 
}
\usage{
\method{plot}{ahull}(x, add = FALSE, do.shape = FALSE, 
	wlines = c("none", "both", "del", "vor"), wpoints = TRUE, 
	number = FALSE, col = NULL, xlim = NULL, 
	ylim = NULL, lwd = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{Object of class \code{"ahull"}.}
\item{add}{Logical, if TRUE add to a current plot.} 
\item{do.shape}{Logical, indicates if the \eqn{\alpha}-shape should also be plotted.}
\item{wlines}{"Which lines?". I.e. should the Delaunay triangulation be plotted (wlines='del'), should the Voronoi diagram be plotted (wlines='vor'), should both be plotted (wlines='both'), or none (wlines='none', the default)?}
 \item{wpoints}{Logical, indicates if sample points should be plotted.}
\item{number}{Logical, defaulting to FALSE; if TRUE then the points plotted will be labelled with their index numbers.}  
 \item{col}{The colour numbers for plotting the \eqn{\alpha}-convex hull, \eqn{\alpha}-shape, data points, Delaunay triangulation, Voronoi diagram, and the point numbers, in that order; defaults to c(1,1,1,1,1,1). If fewer than six numbers are given, they are recycled. (If more than six numbers are given, the redundant ones are ignored.)}
\item{xlim}{The limits on the x-axis.}  
\item{ylim}{The limits on the y-axis.}
  \item{lwd}{The line widths for plotting the tesselations, the \eqn{\alpha}-shape, and the \eqn{\alpha}-convex hull, in that order; defaults to c(1,1,2).}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).}
}

\seealso{\code{\link{ahull}}, \code{\link{ashape}}.}
\examples{
\dontrun{
# Random sample in the unit square
x <- matrix(runif(100), nc = 2)
# Value of alpha
alpha <- 0.2
# alpha-convex hull
ahull.obj <- ahull(x, alpha = alpha)
# Plot including the alpha-convex hull in pink, alpha-shape in blue, 
# sample points in black, voronoi diagram in green 
# and Delaunay triangulation in red 
plot(ahull.obj, do.shape = TRUE, wlines = "both", col = c(6, 4, 1, 2, 3))

# Random sample  from a uniform distribution on a Koch snowflake 
# with initial side length 1 and 3 iterations
x <- rkoch(2000, side = 1, niter = 3)
# Value of alpha
alpha <- 0.05
# Alpha-convex hull
ahull.obj <- ahull(x, alpha = alpha)
plot(ahull.obj)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{nonparametric}