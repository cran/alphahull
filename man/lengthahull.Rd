\name{lengthahull}
\alias{lengthahull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Length of the boundary of the alpha-convex hull}
\description{
   This function calculates the length of the boundary of the \eqn{\alpha}-convex hull of a given sample.
}
\usage{
lengthahull(ahull.arcs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ahull.arcs}{Output matrix of arcs returned by \code{\link{ahull}}.}
}
\details{
The function \code{lengthahull} is internally called by the function \code{\link{ahull}}. 
}
\value{
  \item{length}{Length of the boundary of the \eqn{\alpha}-convex hull.}
}
\seealso{\code{\link{ahull}}.}
\examples{
\dontrun{
# Random sample in the unit square
x <- matrix(runif(100), nc = 2)
# Value of alpha
alpha <- 0.2
# alpha-convex hull
ahull.obj <- ahull(x, alpha = alpha)
# Length of the alpha-convex hull
ahull.obj$length
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}