\name{anglesArc}
\alias{anglesArc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Angles of the extremes of an arc}
\description{
  Given a vector \eqn{v} and an angle \eqn{\theta}, \code{anglesArc} returns the angles that \eqn{A_\theta v} and \eqn{A_{-\theta} v} form with the axis \emph{OX}, where \eqn{A_\theta v} represents the clockwise rotation of angle \eqn{\theta} of the vector \eqn{v}.
}
\usage{
anglesArc(v, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{Vector \eqn{v} in the plane.}
  \item{theta}{Angle \eqn{\theta} (in radians).}
}
\details{
   The angle that forms the vector \eqn{v} with the axis \emph{OX} takes its value in \eqn{[0,2\pi)}.
}

\value{
  \item{angs}{Numeric vector with two components.}
}
\examples{
# Let v=c(0,1) and theta=pi/4
# Consider the arc such that v is the internal angle bisector that 
# divides the angle 2*theta into two equal angles
# The angles that the arc forms with the OX axis are pi/4 and 3*pi/4
v<-c(0,1)
theta<-pi/4
anglesArc(v,theta)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}

