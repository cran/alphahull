\name{delvor}
\alias{delvor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Delaunay triangulation and Voronoi diagram}
\description{
   This function returns a matrix with information about the Delaunay triangulation and Voronoi diagram of a given sample. 
}
\usage{
delvor(x, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x, y}{The \code{x} and \code{y} arguments provide the \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
}
\details{
An attempt is made to interpret the arguments x and y in a way suitable for computing the Delaunay triangulation and Voronoi diagram . Any reasonable way of defining the coordinates is acceptable, see \code{\link{xy.coords}}. 

The function \code{\link[interp]{tri.mesh}} from package \pkg{interp} calculates the Delaunay triangulation of at least three non-collinear points. Using the Delaunay triangulation, the function \code{delvor} calculates the correspondig Voronoi diagram. For each edge of the Delaunay triangulation there is a segment in the Voronoi diagram, given by the union of the circumcenters of the two neighbour triangles that share the edge. For those triangles with edges on the convex hull, the corresponding line in the Voronoi diagram is a semi-infinite segment, whose boundless extreme is calculated by the function \code{\link{dummycoor}}. The function \code{delvor} returns the sample, the output object of class \code{"triSht"} from the function \code{\link[interp]{tri.mesh}} and a matrix \code{mesh} with all the necessary information of the Delaunay triangulation and Voronoi diagram. Thus, for each edge of the Delaunay triangulation the output matrix contains the indexes and coordinates of the sample points that form the edge, the indexes and coordinates of the extremes of the corresponding segment in the Voronoi diagram, and an indicator that takes the value 1 for those extremes of the Voronoi diagram that represent a boundless extreme.
}
\value{
   A list with the following components:
\item{mesh}{A \eqn{n.edges}-row matrix, where \eqn{n.edges} is the total number of different edges of the Delaunay triangulation.}
  \item{x}{A 2-column matrix with the coordinates of the sample points.}
    \item{tri.obj}{Object of class \code{"tri"}. See \code{\link[interp]{tri.mesh}} in package \pkg{interp}.}
}
\references{Renka, R. J. (1996). Algorithm 751: TRIPACK: a constrained two-dimensional Delaunay triangulation package, \emph{ACM Trans. Math. Softw.}, 22(1), pp.1-8.}
\seealso{\code{\link{plot.delvor}}.}
\examples{
\dontrun{
# Random sample in the unit square
x <- matrix(runif(20), nc = 2)
# Delaunay triangulation and Voronoi diagram calculation
delvor.obj <- delvor(x)
}
}
\keyword{nonparametric}