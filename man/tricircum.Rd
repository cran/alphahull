\name{tricircum}
\alias{tricircum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Returns the circumcentres of a Delaunay triangulation}
\description{
 This function converts a triangulation data structure to a triangle list with vertexes, opposite triangles and arcs indexes. It also returns the coordinates of the circumcentre of each triangle. 
}
\usage{
tricircum(tri.obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tri.obj}{Object of class \code{"tri"}, see \code{\link[tripack]{tri.mesh}} from package \pkg{tripack}.}
}
\details{
The function \code{\link{tricircum}} calls the Fortran function \code{delaunaycircum} and it is internally used in the function \code{\link{delvor}}. 
}
\value{
\item{tri.info}{A \eqn{n.tri}-row matrix, where \eqn{n.tri} is the total number of different triangles of the Delaunay triangulation. For each row \code{i}, \code{tri.info[i,]} contains the vertex nodal indexes (first three columns), the indexes of the  triangles and arcs which are opposite to the triangle \code{i} (columns 4-9) and the coordinates of the circumcentre of the triangle \code{i} (columns 10-11).}  
}
\references{Renka, R. J. (1996). Algorithm 751: TRIPACK: a constrained two-dimensional Delaunay triangulation package, \emph{ACM Trans. Math. Softw.}, 22(1), pp.1-8.}
\examples{
# Random sample in the unit square
x <- runif(20)
y <- runif(20)
del <- tricircum(tri.mesh(x, y))
}
\keyword{nonparametric}