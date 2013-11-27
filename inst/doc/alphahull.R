### R code from vignette source 'alphahull.rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: alphahull.rnw:78-79
###################################################
library(alphahull)


###################################################
### code chunk number 2: alphahull.rnw:119-134
###################################################
# Uniform sample of size n=300 in the disc B(c,0.5)\B(c,0.25) 
# with c=(0.5,0.5). 
n<-150
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
x<-data
alpha<-0.1


###################################################
### code chunk number 3: alphahull.rnw:140-143
###################################################
par(mfrow=c(1,1))
print(plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F))
print(plot(ahull(x,alpha=alpha),col=c(6,1,1),xlab="",ylab="",add=T))


###################################################
### code chunk number 4: alphahull.rnw:158-161
###################################################
par(mfrow=c(1,1))
print(plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F))
print(plot(ashape(x,alpha=alpha),col=c(4,1),xlab="",ylab="",add=T))


###################################################
### code chunk number 5: alphahull.rnw:170-184
###################################################
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
x<-data


###################################################
### code chunk number 6: alphahull.rnw:189-197
###################################################
par(mfrow=c(1,3))
alpha1=0.02
alpha2=0.25
alpha3=1
plot(ashape(x,alpha=alpha1),col=c(4,1),xlab="",ylab="",main=expression(paste(alpha," = 0.02 ")))
plot(ashape(x,alpha=alpha2),col=c(4,1),xlab="",ylab="",main=expression(paste(alpha," = 0.25")))
plot(ashape(x,alpha=alpha3),col=c(4,1),xlab="",ylab="",main=expression(paste(alpha," = 1")))
par(mfrow=c(1,1))


###################################################
### code chunk number 7: alphahull.rnw:212-216
###################################################
#x<-matrix(runif(20),nc=2)
x1<-c(0.5915,0.6230,0.9689,0.8248,0.9392,0.8156,0.2050,0.9757,0.0957,0.4139)
y1<-c(0.472,0.619,0.304,0.197,0.716,0.575,0.507,0.574,0.996,0.893)
x<-cbind(x1,y1)


###################################################
### code chunk number 8: alphahull.rnw:222-225
###################################################
par(mfrow=c(1,1))
print(plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F))
print(plot(delvor(x),col=1:3,xlab="",ylab="",add=T))


###################################################
### code chunk number 9: alphahull.rnw:240-244
###################################################
xdeldir<-c(0.7016204,0.8409771,0.1298654,0.6324295,0.9755833,0.9558375,0.3511134,0.9891238,0.4572095,0.545416)
ydeldir<-c(0.1698826,0.3654726,0.5892516,0.6108104,0.8855026,0.9910547,0.2433913,0.1406111,0.4103610,0.4450315)
xydeldir<-cbind(xdeldir,ydeldir)
delvor.obj<-delvor(xydeldir)


###################################################
### code chunk number 10: alphahull.rnw:248-250
###################################################
xtri<-xdeldir
ytri<-ydeldir


###################################################
### code chunk number 11: alphahull.rnw:255-262
###################################################
par(mfrow=c(1,2))
print(plot(voronoi.mosaic(xtri,ytri),do.points=FALSE,main="",sub="",col=3,xlim=c(-0.5,1.5),ylim=c(-0.5,1.5)))
axis(1)
axis(2)
print(plot(tri.mesh(xtri,ytri),col=2,add=T))
points(xtri,ytri,col=1)
print(plot(delvor(xtri,ytri),col=c(1,2,3),xlim=c(-0.5,1.5),ylim=c(-0.5,1.5),xlab="",ylab=""))


###################################################
### code chunk number 12: alphahull.rnw:270-274
###################################################
x <- c(0.905, 0.606, 0.458, 0.988, 0.744)
y <- c(0.763, 0.937, 0.095, 0.259, 0.731) 
dv <- delvor(x,y)
dv


###################################################
### code chunk number 13: alphahull.rnw:280-281
###################################################
plot(dv, main = "Delaunay triangulation and Voronoi diagram",col = 1:3, xlab = "x-coordinate", ylab = "y-coordinate",xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5), number = TRUE)


###################################################
### code chunk number 14: alphahull.rnw:288-290
###################################################
par(mfrow=c(1,1))
print(plot(dv, main = "Delaunay triangulation and Voronoi diagram",col = 1:3, xlab = "x-coordinate", ylab = "y-coordinate",xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5), number = TRUE))


###################################################
### code chunk number 15: alphahull.rnw:322-325
###################################################
x<-matrix(runif(40),ncol=2)
alpha <- 0.2
alphashape <- ashape(x, alpha = alpha)


###################################################
### code chunk number 16: alphahull.rnw:328-332
###################################################
names(alphashape)
alphashape$alpha.extremes
alphashape$edges[, 1:2]
alphashape$length


###################################################
### code chunk number 17: alphahull.rnw:337-338
###################################################
plot(alphashape, col = c(4, 1), xlab = "x-coordinate", ylab = "y-coordinate",number= TRUE, main = expression(paste(alpha, "-shape")))


###################################################
### code chunk number 18: alphahull.rnw:344-346
###################################################
par(mfrow=c(1,1))
print(plot(alphashape, col = c(4, 1), xlab = "x-coordinate", ylab = "y-coordinate",number= TRUE, main = expression(paste(alpha, "-shape"))))


###################################################
### code chunk number 19: alphahull.rnw:354-355
###################################################
plot(alphashape, wlines = "del", col = c(4, 1, 2), xlab = "x-coordinate", ylab = "y-coordinate")


###################################################
### code chunk number 20: alphahull.rnw:361-363
###################################################
par(mfrow=c(1,1))
print(plot(alphashape, wlines = "del", col = c(4, 1, 2), xlab = "x-coordinate",ylab = "y-coordinate",main = expression(paste(alpha, "-shape and Delaunay triangulation"))))


###################################################
### code chunk number 21: alphahull.rnw:417-443
###################################################
par(mfrow=c(1,1))
plot(0,type="n",axes=FALSE,xlim=c(0,0.5),ylim=c(0,0.5),xlab="",ylab="")
r<-0.5
t<-0
segments(0,0,r*cos(t),r*sin(t),col=4,lty=2)
points(r*cos(t),r*sin(t),pch=19,col=4)

t<-pi/6
arrows(0,0,0.3*cos(t),0.3*sin(t))
v<-c(0.3*cos(t),0.3*sin(t))
v<-v/sqrt(sum(v^2))
arc(c(0,0),0.5,v,pi/6,col=4,lwd=2)
t<-pi/3
segments(0,0,r*cos(t),r*sin(t),col=4,lty=2)
points(r*cos(t),r*sin(t),pch=19,col=4)

t<-pi/12
v<-c(0.3*cos(t),0.3*sin(t))
v<-v/sqrt(sum(v^2))
arc(c(0,0),0.3,v,pi/12,col=1)


text(0.15,0.11,expression(italic(v)),cex=1.5)
text(0.31,0.07,expression(italic(theta)),cex=1.5)
text(-0.015,0,expression(italic(c)),cex=1.5)
text(0.15,0.3,expression(italic(r)),cex=1.5)


###################################################
### code chunk number 22: alphahull.rnw:452-458
###################################################
n <- 200
theta<-runif(n,0,2*pi)
r<-sqrt(runif(n,0.25^2,0.5^2))
x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
alpha <- 0.15
alphahull <- ahull(x, alpha = alpha)


###################################################
### code chunk number 23: alphahull.rnw:461-465
###################################################
names(alphahull)
alphahull$complement[1:5,1:3]
alphahull$arcs[1:5,]
alphahull$length


###################################################
### code chunk number 24: alphahull.rnw:471-472
###################################################
plot(alphahull, col = c(6, rep(1, 5)), xlab = "x-coordinate",ylab = "y-coordinate", main = expression(paste(alpha, "-hull")))


###################################################
### code chunk number 25: alphahull.rnw:478-480
###################################################
par(mfrow=c(1,1))
print(plot(alphahull,col=c(6,rep(1,5)), xlab = "x-coordinate",ylab = "y-coordinate", main = expression(paste(alpha, "-hull"))))


###################################################
### code chunk number 26: alphahull.rnw:488-493
###################################################
plot(alphahull, col = c(6, rep(1, 5)), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-hull")))
warcs<- which(alphahull$arcs[,3]>0)
for (i in warcs) {
arc(alphahull$arcs[i, 1:2], alphahull$arcs[i,3], c(0,1), pi, col = "gray", lty = 2)
}


###################################################
### code chunk number 27: alphahull.rnw:498-499
###################################################
plot(alphahull, do.shape = TRUE, col = c(6, 4, rep(1, 4)), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-hull and ",alpha, "-shape")))


###################################################
### code chunk number 28: alphahull.rnw:505-513
###################################################
par(mfrow=c(1,2))
print(plot(alphahull, col = c(6, rep(1, 5)), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-hull"))))
warcs<- which(alphahull$arcs[,3]>0)
for (i in warcs) {
arc(alphahull$arcs[i, 1:2], alphahull$arcs[i,3], c(0,1), 2 * pi, col = "gray", lty = 2)
}
print(plot(alphahull, do.shape = TRUE, col = c(6, 4, rep(1, 4)),xlab = "x-coordinate", ylab = "y-coordinate",main = expression(paste(alpha, "-hull and ",alpha, "-shape")))
)


###################################################
### code chunk number 29: alphahull.rnw:521-523
###################################################
inahull(alphahull,c(0.5,0.5))
inahull(alphahull,c(0.6,0.2))


###################################################
### code chunk number 30: alphahull.rnw:528-529
###################################################
areaahull(alphahull)


###################################################
### code chunk number 31: alphahull.rnw:541-547
###################################################
par(mfrow=c(2,2))
for (k in 1:4){
snow <- koch(side = 1, niter = k)
plot(snow[, 1], snow[, 2], type = "l", xlab = "", ylab = "",axes=F,asp=TRUE)
polygon(snow[, 1], snow[, 2] , col = 4)
}


###################################################
### code chunk number 32: alphahull.rnw:553-559
###################################################
par(mfrow=c(2,2))
for (i in 1:4){
snow <- koch(side = 1, niter = i)
plot(snow[, 1], snow[, 2], type = "l", xlab = "", ylab = "",axes=F,asp=TRUE)
polygon(snow[, 1], snow[, 2] , col = 4)
}


###################################################
### code chunk number 33: alphahull.rnw:567-570
###################################################
unifkoch <- rkoch(2000, side = 1, niter = 3)
alpha=0.05
alphahull <- ahull(unifkoch, alpha = alpha)


###################################################
### code chunk number 34: alphahull.rnw:575-577
###################################################
alphahull$length
alphahull$ashape$length


###################################################
### code chunk number 35: alphahull.rnw:583-586
###################################################
par(mfrow=c(1,2))
plot(alphahull$ashape, col = c(4, 1), xlab = "x-coordinate", ylab = "y-coordinate",main = expression(paste(alpha, "-shape")),asp=TRUE)
plot(alphahull, col = c(6, rep(1, 5)), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-hull")),asp=TRUE)


###################################################
### code chunk number 36: alphahull.rnw:595-607
###################################################
unifkoch <- rkoch(2000, side = 1, niter = 3)
alpha1=0.15
alphashape1 <- ashape(unifkoch, alpha = alpha1)
alpha2=0.05
alphashape2 <- ashape(alphashape1$delvor.obj, alpha = alpha2)
alpha3=0.02
alphashape3 <- ashape(alphashape1$delvor.obj, alpha = alpha3)

par(mfrow=c(1,3))
plot(alphashape1, col = c(4, 1), xlab = "", ylab = "",asp=TRUE)
plot(alphashape2, col = c(4, 1), xlab = "", ylab = "",asp=TRUE)
plot(alphashape3, col = c(4, 1), xlab = "", ylab = "",asp=TRUE)


