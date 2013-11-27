\name{inter}
\alias{inter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Intersection of two circumferences}
\description{
   This function calculates the intersection of two circumferences, given their centers and radius \eqn{c1,r1} and \eqn{c2,r2}, respectively.
}
\usage{
inter(c11, c12, r1, c21, c22, r2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c11}{\emph{X}-coordinate of the center \eqn{c1}.}
  \item{c12}{\emph{Y}-coordinate of the center \eqn{c1}.}
  \item{r1}{ Radius \eqn{r1}.}
  \item{c21}{\emph{X}-coordinate of the center \eqn{c2}.}
  \item{c22}{\emph{Y}-coordinate of the center \eqn{c2}.}
  \item{r2}{ Radius \eqn{r2}.}
}
\details{
The function \code{inter} is internally called by the function \code{\link{ahull}}. 
}
\value{
  A list with the following components:
  \item{n.cut}{Number of intersection points (0,1,2, or Inf).}
  \item{v1}{If there are two intersection points, \code{v1} is the numeric vector whose components are the coordinates of the unitary vector that has its origin in \eqn{c1} and it's perpendicular to the chord that joins the intersection points of the two circumferences. Otherwise, \code{v1=(0,0)}}
 \item{theta1}{Angle that forms \code{v1} with the radius that joins the center \eqn{c1} with an intersection point. }
\item{v2}{If there are two intersection points, \code{v2} is the numeric vector whose components are the coordinates of the unitary vector that has its origin in \eqn{c2} and it's perpendicular to the chord that joins the intersection points of the two circumferences. Otherwise, \code{v2=(0,0)}}
 \item{theta2}{Angle that forms \code{v2} with the radius that joins the center \eqn{c2} with an intersection point.}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
