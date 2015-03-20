\name{ahull_track}
\alias{ahull_track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{alpha-convex hull calculation of tracking data}
\description{
This function approximates the \eqn{\alpha}-convex hull of tracking data and returns a list of geom_path objects of the boundary.
}
\usage{
ahull_track(x, y = NULL, alpha, nps = 10000, sc = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{The \code{x} and \code{y} arguments provide the \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
  \item{alpha}{Value of \eqn{\alpha}.}
  \item{nps}{Number of points to generate in each segment connecting two locations, see Details}
  \item{sc}{Scale factor.}
}
\details{

An attempt is made to interpret the arguments x and y in a way suitable for computing the \eqn{\alpha}-convex hull. Any reasonable way of defining the coordinates is acceptable, see \code{\link{xy.coords}}. 

Increase \code{nps} if the trajectory is not contained in the computed estimator.
}
\value{
A list of geom_path objects defining the boundary of the \eqn{\alpha}-convex
}
\references{Cholaquidis, A., Fraiman, R., Lugosi, G. and Pateiro-Lopez, B. (2014) Set estimation from reflected Brownian motion. \emph{arXiv:1411.0433}.

Wikelski, M., and Kays, R. (2014). Movebank: archive, analysis and sharing of animal movement data. World Wide Web electronic publication.
}
\examples{
\dontrun{
library(move)
library(ggmap)
# Data from Movebank
# Study Name: Dunn Ranch Bison Tracking Project
# Principal Investigator: Stephen Blake, Randy Arndt, Doug Ladd
# Max Planck Institute for Ornithology Radolfzell Germany
study <- "Dunn Ranch Bison Tracking Project" 
cainfo <- system.file("CurlSSL", "cacert.pem", package = "RCurl")
options(RCurlOptions = list(verbose = FALSE, capath = cainfo, ssl.verifypeer = FALSE))
# Login to movebank (first create the login object) 
curl <- movebankLogin(username = "xxx", password = "zzz") 
# Downloads study stored in Movebank
track <- getMovebankData(study = study, login = curl) 
dat <- track@data[track@data[, "deployment_id"] == 13848432,]
# Map of animal locations 
bbox <- ggmap::make_bbox(dat[,"location_long"], dat[,"location_lat"], f = 0.3) 
map_loc <- get_map(location = bbox, source = "google", maptype = 'satellite')  
map <- ggmap(map_loc, extent = 'panel', maprange=FALSE) 
p <- map + geom_path(data = dat, aes(x = location_long, y = location_lat), col=2, size=0.3)
p
ah_gp <- ahull_track(x = dat[, c("location_long", "location_lat")], alpha = 0.005)
p + ah_gp
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
