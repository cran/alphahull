\name{dw_track}
\alias{dw_track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{RBM-sausage calculation of tracking data}
\description{
This function approximates the RBM-sausage of tracking data and returns a list of geom_path objects of the boundary.
}
\usage{
dw_track(x, y = NULL, eps, nps = 20000, sc = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{The \code{x} and \code{y} arguments provide the \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
  \item{eps}{Value of \eqn{\epsilon}.}
  \item{nps}{Number of points to generate in each segment connecting two locations, see Details.}
  \item{sc}{Scale factor.}
}
\details{

An attempt is made to interpret the arguments x and y in a way suitable for computing the RBM-sausage. Any reasonable way of defining the coordinates is acceptable, see \code{\link{xy.coords}}. 
Given a trajectory in the plane, the estimator is defined as the set of points whose distance to the trajectory is less than or equal to \eqn{\epsilon} (this estimator is analogous to the one of Devroye and Wise (1980) for i.i.d. data). Increase \code{nps} if the trajectory is not contained in the computed estimator.
}
\value{
A list of geom_path objects defining the boundary of the estimator
}
\references{Cholaquidis, A., Fraiman, R., Lugosi, G. and Pateiro-Lopez, B. (2014) Set estimation from reflected Brownian motion. \emph{arXiv:1411.0433}.

Devroye, L. and Wise, G. (1980) Detection of abnormal behaviour via nonparametric estimation of the support. \emph{SIAM J. Appl. Math.} 3, pp. 480-488.

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
ah_dw <- dw_track(x = dat[, c("location_long", "location_lat")], eps = 0.001)
p + ah_dw
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
