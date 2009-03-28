\name{arc}
\alias{arc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add an arc to a plot}
\description{
   This function adds the arc of \eqn{B(c,r)} between the angles that \eqn{A_\theta v} and \eqn{A_{-\theta} v} form with the axis \emph{OX}, where \eqn{A_\theta v} represents the clockwise rotation of angle \eqn{\theta} of the vector \eqn{v}.
}
\usage{
arc(c, r, v, theta, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c}{ Center \eqn{c} of the arc. }
  \item{r}{ Radius \eqn{r} of the arc. }
  \item{v}{ Vector \eqn{v} in the plane. }
  \item{theta}{ Angle \eqn{\theta} (in radians).}
  \item{\dots}{ Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).}
}
\seealso{\code{\link{plot.ahull}}.}
\examples{
# Plot of the circumference of radius 1
theta<-seq(0,2*pi,length=100)
r<-1
plot(r*cos(theta),r*sin(theta),type="l")
# Add in red the arc between pi/4 and 3*pi/4
arc(c(0,0),1,c(0,1),pi/4,col=2,lwd=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}