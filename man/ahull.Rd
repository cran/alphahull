\name{ahull}
\alias{ahull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{alpha-convex hull calculation}
\description{
This function calculates the \eqn{\alpha}-convex hull of a given sample of points in the plane for \eqn{\alpha>0}.
}
\usage{
ahull(x, y = NULL, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{The \code{x} and \code{y} arguments provide the \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
  \item{alpha}{Value of \eqn{\alpha}.}
}
\details{

An attempt is made to interpret the arguments x and y in a way suitable for computing the \eqn{\alpha}-convex hull. Any reasonable way of defining the coordinates is acceptable, see \code{\link{xy.coords}}. 

The \eqn{\alpha}-convex hull is defined for any finite number of points. However, since the algorithm is based on the Delaunay triangulation, at least three non-collinear points are required. 

If \code{y} is NULL and \code{x} is an object of class \code{"delvor"}, then the \eqn{\alpha}-convex hull is computed with no need to invoke again the function \code{\link{delvor}} (it reduces the computational cost). 

The complement of the \eqn{\alpha}-convex hull can be written as the union of \eqn{O(n)} open balls and halfplanes, see \code{\link{complement}}. 
The boundary of the \eqn{\alpha}-convex hull is formed by arcs of open balls of radius \eqn{\alpha} (besides possible isolated sample points). The arcs are determined by the intersections of some of the balls that define the complement of the \eqn{\alpha}-convex hull. The extremes of an arc are given by \eqn{c+rA_\theta v} and \eqn{c+rA_{-\theta}v} where \eqn{c} and \eqn{r} represent the center and radius of the arc, repectively, and \eqn{A_\theta v} represents the clockwise rotation of angle \eqn{\theta} of the unitary vector \eqn{v}. Joining the end points of adjacent arcs we can define polygons that help us to determine the area of the estimator , see \code{\link{areaahull}}.
}
\value{
 A list with the following components:
  \item{arcs}{For each arc in the boundary of the \eqn{\alpha}-convex hull, the columns of the matrix \code{arcs} store the center \eqn{c} and radius \eqn{r} of the arc, the unitary vector \eqn{v}, the angle \eqn{\theta} that define the arc and the indices of the end points, see Details. The coordinates of the end points of the arcs are stored in \code{xahull}. For isolated points in the boundary of the \eqn{\alpha}-convex hull, columns 3 to 6 of the matrix \code{arcs} are equal to zero. }
  \item{xahull}{A 2-column matrix with the coordinates of the original set of points besides possible new end points of the arcs in the boundary of the \eqn{\alpha}-convex hull.}
\item{length}{Length of the boundary of the \eqn{\alpha}-convex hull, see \code{\link{lengthahull}}.}  
\item{complement}{Output matrix from \code{\link{complement}}.}
\item{alpha}{Value of \eqn{\alpha}.}
\item{ashape.obj}{Object of class \code{"ashape"} returned by the function \code{\link{ashape}}.}
}
\references{Edelsbrunner, H., Kirkpatrick, D.G. and Seidel, R. (1983). On the shape of a set of points in the plane. \emph{IEEE Transactions on Information Theory}, 29(4), pp.551-559.

Rodriguez-Casal, R. (2007). Set estimation under convexity type assumptions. \emph{Annales de l'I.H.P.- Probabilites \& Statistiques}, 43, pp.763-774.

Pateiro-Lopez, B. (2008). \emph{Set estimation under convexity type restrictions}. Phd. Thesis. Universidad de Santiago de Compostela. ISBN 978-84-9887-084-8.}
\seealso{\code{\link{plot.ahull}}.}
\examples{
# Random sample in the unit square
x <- matrix(runif(100), nc = 2)
# Value of alpha
alpha <- 0.2
# Alpha-convex hull
ahull.obj <- ahull(x, alpha = alpha)
plot(ahull.obj)

# Uniform sample of size n=300 in the annulus B(c,0.5)\B(c,0.25), 
# with c=(0.5,0.5). 
n <- 300
theta<-runif(n,0,2*pi)
r<-sqrt(runif(n,0.25^2,0.5^2))
x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
# Value of alpha
alpha <- 0.1
# Alpha-convex hull
ahull.obj <- ahull(x, alpha = alpha)
# The arcs defining the boundary of the alpha-convex hull are ordered
plot(x)
for (i in 1:dim(ahull.obj$arcs)[1]){
arc(ahull.obj$arcs[i,1:2],ahull.obj$arcs[i,3],ahull.obj$arcs[i,4:5],
ahull.obj$arcs[i,6],col=2)
Sys.sleep(0.5)
}

# Random sample  from a uniform distribution on a Koch snowflake 
# with initial side length 1 and 3 iterations
x <- rkoch(2000, side = 1, niter = 3)
# Value of alpha
alpha <- 0.05
# Alpha-convex hull
ahull.obj <- ahull(x, alpha = alpha)
plot(ahull.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
