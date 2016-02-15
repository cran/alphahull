\name{ashape}
\alias{ashape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{alpha-shape calculation}
\description{
This function calculates the \eqn{\alpha}-shape of a given sample for \eqn{\alpha>0}.
}
\usage{
ashape(x, y = NULL, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{The \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
  \item{alpha}{Value of \eqn{\alpha}.}
}
\details{

An attempt is made to interpret the arguments x and y in a way suitable for computing the \eqn{\alpha}-shape, see \code{\link{xy.coords}}. 

The \eqn{\alpha}-shape is defined for any finite number of points. However, since the algorithm is based on the Delaunay triangulation, at least three non-collinear points are required. 

If \code{y} is NULL and \code{x} is an object of class \code{"delvor"}, then the \eqn{\alpha}-shape is computed without invoking again the function \code{\link{delvor}} (it reduces the computational cost). 

The function \code{\link{ashape}} returns (among other values) the matrix \code{edges}. The structure of \code{edges} is that of matrix \code{mesh} returned by the function \code{\link{delvor}}. Note that the \eqn{\alpha}-shape is a subgraph of the Delaunay triangulation and, therefore, \code{edges} is a submatrix of \code{mesh}.
}
\value{
 A list with the following components:
\item{edges}{A \emph{n.seg}-row matrix with the coordinates and indexes of the edges of the Delaunay triangulation that form the \eqn{\alpha}-shape. The number of rows \emph{n.seg} coincides with the number of segments of the \eqn{\alpha}-shape. The matrix also includes information of the Voronoi extremes corresponding to each segment.}
\item{length}{Length of the \eqn{\alpha}-shape.}  
\item{alpha}{Value of \eqn{\alpha}.}
\item{alpha.extremes}{Vector with the indexes of the sample points that are \eqn{\alpha}-extremes. See Edelsbrunnner \emph{et al.} (1983).}
\item{delvor.obj}{Object of class \code{"delvor"} returned by the \code{\link{delvor}} function.}
\item{x}{A 2-column matrix with the coordinates of the set of points.}
}

\references{Edelsbrunner, H., Kirkpatrick, D.G. and Seidel, R. (1983). On the shape of a set of points in the plane. \emph{IEEE Transactions on Information Theory}, 29(4), pp.551-559.}
\seealso{\code{\link{plot.ashape}}, \code{\link{delvor}}.}
\examples{
\dontrun{
# Uniform sample of size n=300 in the annulus B(c,0.5)\B(c,0.25), 
# with c=(0.5,0.5). 
n <- 300
theta<-runif(n,0,2*pi)
r<-sqrt(runif(n,0.25^2,0.5^2))
x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
# Value of alpha
alpha <- 0.1
# alpha-shape 
ashape.obj <- ashape(x, alpha = alpha)
# If we change the value of alpha there is no need to compute
# again the Delaunay triangulation and Voronoi Diagram
alpha <- 0.4
ashape.obj.new <- ashape(ashape.obj$delvor.obj, alpha = alpha)

# Random sample  from a uniform distribution on a Koch snowflake 
# with initial side length 1 and 3 iterations
x <- rkoch(2000, side = 1, niter = 3)
# Value of alpha
alpha <- 0.05
# alpha-shape 
ashape.obj <- ashape(x, alpha = alpha)
}
}
\keyword{nonparametric}