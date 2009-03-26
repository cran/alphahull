\name{dummycoor}
\alias{dummycoor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Semi-infinite edge of the Voronoi diagram}
\description{
   This function determines fictitious coordinates for the boundless extreme of a semi-infinite edge of the Voronoi diagram.
}
\usage{
dummycoor(tri.obj, l1, l2, m, away)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tri.obj}{Object of class \code{"tri"}. See \code{\link[tripack]{tri.mesh}} in package \pkg{tripack}.}
  \item{l1}{ Index of the sample point correponding to one vertex of a triangle of Delaunay that lies on the convex hull, see Details.}
  \item{l2}{ Index of the sample point correponding to other vertex of a triangle of Delaunay that lies on the convex hull, see Details.}
  \item{m}{ Index of the circumcenter of the triangle of Delaunay with one edge on the convex hull.}
  \item{away}{ Constant that determines how far away the fictitious boundless extreme is located.}
}
\details{
   When a triangle of the Delaunay triangulation has one of its edges (given by the segment that joins the sample points with indexes \code{l1} and \code{l2}) on the convex hull, the corresponding segment of the Voronoi diagram is semi-infinite. The finite extreme coincides with the circumcenter of the triangle and the direction of the line is given by the perpendicular bisector of the edge that lies on the convex hull. 
}
\value{
  \item{dum}{Fictitious coordinates of the boundless extreme.}
}
\seealso{ \code{\link{delvor}}.}
\keyword{nonparametric}