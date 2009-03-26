\name{rotation}
\alias{rotation}
%- Also NEED an '\alias' for EACH other topic documented here.
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clockwise rotation}
\description{
   This function calculates the clockwise rotation of angle \eqn{\theta} of a given vector \eqn{v} in the plane.
}
\usage{
rotation(v, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{v}{Vector \eqn{v} in the plane.}
  \item{theta}{Angle \eqn{\theta} (in radians).}
}
\value{
 \item{v.rot}{Vector after rotation.}
 }
\examples{
# Rotation of angle pi/4 of the vector (0,1)
rotation(v=c(0,1),theta=pi/4)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}