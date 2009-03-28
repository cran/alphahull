\name{plot.ashape}
\alias{plot.ashape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the alpha-shape}
\description{
    This function returns a plot of the \eqn{\alpha}-shape. If desired, it also adds the Delaunay triangulation and Voronoi diagram of the sample. 

}
\usage{
\method{plot}{ashape}(x, add = FALSE, wlines = c("none", "both", "del", "vor"),
	wpoints = TRUE, number = FALSE, col = NULL, 
	xlim = NULL, ylim = NULL, lwd = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{Object of class \code{"ashape"}.}
\item{add}{Logical, if TRUE add to a current plot.} 
 \item{wlines}{"Which lines?". I.e. should the Delaunay triangulation be plotted (wlines='del'), should the Voronoi diagram be plotted (wlines='vor'), should both be plotted (wlines='both'), or none (wlines='none', the default)?}
 \item{wpoints}{Logical, indicates if sample points should be plotted.}
\item{number}{Logical, defaulting to FALSE; if TRUE then the points plotted will be labelled with their index numbers.}  
 \item{col}{The colour numbers for plotting the \eqn{\alpha}-shape, data points, Delaunay triangulation, Voronoi diagram, and the point numbers, in that order; defaults to c(1,1,1,1,1). If fewer than five numbers are given, they are recycled. (If more than five numbers are given, the redundant ones are ignored.)}
\item{xlim}{The limits on the x-axis.}  
\item{ylim}{The limits on the y-axis.}
  \item{lwd}{The line widths for plotting the tesselations and the \eqn{\alpha}-shape; defaults to c(1,2).}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).}
}
\seealso{\code{\link{ashape}}.}
\examples{
# Uniform sample of size n=300 in the disc B(c,0.5)\B(c,0.25) 
# with c=(0.5,0.5). 
n<-300
m<-0
data<-matrix(0,n,2)
while(m<n){
x<-runif(1)
y<-runif(1)
if(((x-0.5)^2+(y-0.5)^2<=(0.5)^2 )&((x-0.5)^2+(y-0.5)^2>=(0.25)^2)){
m<-m+1
data[m,]<-c(x,y)
}
}
# Value of alpha
alpha<-0.1
# alpha-shape 
ashape.obj<-ashape(data,alpha=alpha)
# Plot alpha-shape in red and Delaunay triangulation in black
plot(ashape.obj)
}
\keyword{nonparametric}