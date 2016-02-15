\name{areaahull}
\alias{areaahull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Area of the alpha-convex hull
}
\description{
This function calculates the area of the \eqn{\alpha}-convex hull of a sample of points.
}
\usage{
areaahull(x, timeout = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"ahull"}.}
  \item{timeout}{A numeric specifying the maximum number of seconds the expression is allowed to run before being interrupted by the timeout.}
}
\value{
\item{area}{Area of the \eqn{\alpha}-convex hull. If the area cannot be computed, the output will be NA with a warning.}
}

\seealso{\code{\link{ahull}}.}
\examples{
\dontrun{
# Random sample in the unit square
x <- matrix(runif(500), nc = 2)
# Value of alpha
alpha <- 1
# alpha-convex hull
ahull.obj <- ahull(x, alpha = alpha)
# Area of the alpha-convex hull
areaahull(ahull.obj)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}