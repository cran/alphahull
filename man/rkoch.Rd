\name{rkoch}
\alias{rkoch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random generation on a Koch snowflake curve}
\description{
This function generates ramdom points from a uniform distribution on a Koch snowflake.
}
\usage{
rkoch(n, side = 3, niter = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of observations.
}
  \item{side}{
Side length of the initial equilateral triangle of the Koch snowflake curve.
}
  \item{niter}{Number of iterations in the development of the snowflake curve.
}
}

\value{
A 2-column matrix with the coordinates of generated points.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{koch}}.}
\examples{
\dontrun{
unifkoch <- rkoch(2000, side = 1, niter = 3)
plot(unifkoch, asp = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
