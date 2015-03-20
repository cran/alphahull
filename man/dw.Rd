\name{dw}
\alias{dw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Devroye-Wise estimator}
\description{
This function calculates the Devroye-Wise estimator of a given sample of points in the plane for \eqn{\epsilon>0}.
}
\usage{
dw(x, y = NULL, eps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{The \code{x} and \code{y} arguments provide the \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
  \item{eps}{Value of \eqn{\epsilon}.}
}
\details{

An attempt is made to interpret the arguments x and y in a way suitable for computing the Devroye-Wise estimator. Any reasonable way of defining the coordinates is acceptable, see \code{\link{xy.coords}}. 
}
\value{
Given a sample of points in the plane, the estimator is defined as union of balls of radius \eqn{\epsilon} with centers in the sample points. For each arc in the boundary of the Devroye-Wise estimator, the columns of the output matrix store the center \eqn{c} and radius \eqn{r} of the arc, the unitary vector \eqn{v}, the angle \eqn{\theta} that define the arc and the indices of the end points.
}
\references{Devroye, L. and Wise, G. (1980) Detection of abnormal behaviour via nonparametric estimation of the support. \emph{SIAM J. Appl. Math.} 3, pp. 480-488.
}
\examples{
\dontrun{
# Uniform sample of size n = 200 in the annulus B(c, 0.5)\B(c, 0.25), 
# with c = (0.5, 0.5).
n <- 200
theta <- runif(n, 0, 2*pi)
r <- sqrt(runif(n, 0.25^2, 0.5^2))
x <- cbind(0.5 + r*cos(theta), 0.5 + r*sin(theta))
eps <- 0.05
dw.obj <- dw(x, eps = eps)
plot(x)
for(i in 1:dim(dw.obj)[1]){arc(dw.obj[i, 1:2], eps, dw.obj[i, 4:5], dw.obj[i, 6])}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
