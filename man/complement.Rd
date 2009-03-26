\name{complement}
\alias{complement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Complement of the alpha-convex hull}
\description{
   This function calculates the complement of the \eqn{\alpha}-convex hull of a given sample for \eqn{\alpha>0}.
}
\usage{
complement(x, y = NULL, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x, y}{The \code{x} and \code{y} arguments provide the \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
  \item{alpha}{Value of \eqn{\alpha}.}
}
\details{
An attempt is made to interpret the arguments x and y in a way suitable for computing the \eqn{\alpha}-shape. Any reasonable way of defining the coordinates is acceptable, see \code{\link{xy.coords}}. 

If \code{y} is NULL and \code{x} is an object of class \code{"delvor"}, then the complement of the \eqn{\alpha}-convex hull is computed with no need to invoke again the function \code{\link{delvor}} (it reduces the computational cost). 

   The complement of the \eqn{\alpha}-convex hull is calculated as a union of open balls and halfplanes that do not contain any point of the sample. See Edelsbrunnner \emph{et al.} (1983) for a basic description of the algorithm. The construction of the complement is based on the Delaunay triangulation and Voronoi diagram of the sample, provided by the function \code{\link{delvor}}. The function \code{\link{complement}} returns a matrix \code{compl}. For each row \code{i},  \code{compl[i,]} contains the information relative to an open ball or halfplane of the complement. The first three columns are assigned to the characterization of the ball or halfplane \code{i}. The information relative to the edge of the Delaunay triangulation that generates the ball or halfplane \code{i} is contained in \code{compl[i,4:16]}. Thus, if the row \code{i} refers to an open ball, \code{compl[i,1:3]} contains the center and radius of the ball. Furthermore,  \code{compl[i,17:18]} and \code{compl[i,19]} refer to the unitary vector \eqn{v} and the angle \eqn{\theta} that characterize the arc that joins the two sample points that define the ball \code{i}. If the row \code{i} refers to a halfplane, \code{compl[i,1:3]} determines its equation. For the halfplane \eqn{y>a+bx}, \code{compl[i,1:3]=(a,b,-1)}. In the same way,   for the halfplane \eqn{y<a+bx}, \code{compl[i,1:3]=(a,b,-2)}, for the halfplane \eqn{x>a}, \code{compl[i,1:3]=(a,0,-3)} and for the halfplane \eqn{x<a}, \code{compl[i,1:3]=(a,0,-4)}.  
}
\value{
  \item{compl}{Output matrix. For each row \code{i},  \code{compl[i,]} contains the information relative to an open ball or halfplane of the complement of the \eqn{\alpha}-convex hull, see Details.}
 }
\references{Edelsbrunner, H., Kirkpatrick, D.G. and Seidel, R. (1983) \emph{On the shape of a set of points in the plane}. IEEE Transactions on Information Theory, Vol IT-29, No. 4.}
\seealso{\code{\link{delvor}}, \code{\link{ahull}}.}
\examples{
# Random sample in the unit square
x<-matrix(runif(100),nc=2)
# Value of alpha
alpha<-0.2
# Complement of the alpha-convex hull
compl<-complement(x,alpha=alpha)
}
\keyword{nonparametric}