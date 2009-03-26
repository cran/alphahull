`dummycoor` <-
function(tri.obj,l1,l2,m,away){
v<-l2-l1
v<-c(v[2],-v[1])
norm<-sum(v^2)
if (norm>0){
v<-v/norm
}
mp<-(l1+l2)/2
# We displace pm in the direction of v and we check if the new point is in the convex hull to determine the direction in which we have to set the infinite edge
eps<-10e-6
test<-mp+eps*v 
# We determine if m is in the convex hull
inconv<-in.convex.hull(tri.obj,test[1],test[2])
if (inconv){
dum<-mp-away*v
}
else{
dum<-mp+away*v
}
return(dum)
}

