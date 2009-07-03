\name{koch}
\alias{koch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct a Kock snowflake curve
}
\description{
This function uses recursion to construct a Kock snowflake curve.
}
\usage{
koch(side = 3, niter = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{Side length of the initial equilateral triangle.
}
  \item{niter}{Number of iterations in the development of the snowflake curve.
}
}
\details{
The Koch snowflake is a fractal curve described by the Swedish mathematician Helge von Koch in 1904. It is built by starting with an equilateral triangle, removing the inner third of each side, building another equilateral triangle at the location where the side was removed, and then repeating the process.
}
\value{
 \item{vertices}{A 2-column matrix with the coordinates of the snowflake vertices.}
}
\references{
von Koch, H. (1904). Sur une courbe continue sans tangente, obtenue par une construction geometrique elementaire. \emph{Arkiv for Matematik}, 1, pp.681-704.
}
\seealso{\code{\link{rkoch}}.}
\examples{
# The first four iterations of a Koch snowflake 
# with side length of the initial equilateral triangle equal to 3.
vertices <- koch(side = 2, niter = 4)
plot(vertices[, 1], vertices[, 2], type = "l", asp = TRUE, 
main = "Koch snowflake", xlab = "", ylab = "", col = 4)
polygon(vertices[, 1], vertices[, 2] , col = 4)
}
\keyword{nonparametric}
