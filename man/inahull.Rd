\name{inahull}
\alias{inahull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine whether a point belongs to the alpha-convex hull}
\description{
   This function determines whether a given point \eqn{p} belongs to the \eqn{\alpha}-convex hull of a sample.
}
\usage{
inahull(ahull.obj, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ahull.obj}{Object of class \code{"ahull"} returned by the funcion \code{\link{ahull}}.}
  \item{p}{ Numeric vector with two components describing a point in the plane.}
}
\details{
   The complement of the \eqn{\alpha}-convex hull of a sample is calculated by \code{\link{complement}}. The function \code{\link{inahull}} checks whether the point \eqn{p} belongs to any of the open balls or halfplanes that define the complement. 
}
\value{
  \item{in.ahull}{A logical value specifying whether the point \eqn{p} belongs to the \eqn{\alpha}-convex hull.}
}
\seealso{\code{\link{ahull}}, \code{\link{complement}}.}
\examples{
\dontrun{
# Random sample in the unit square
x <- matrix(runif(100), nc = 2)
# Value of alpha
alpha <- 0.2
# alpha-convex hull
ahull.obj <- ahull(x, alpha = alpha)
# Check if the point (0.5, 0.5) belongs to the alpha-convex hull
inahull(ahull.obj, p = c(0.5, 0.5))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}