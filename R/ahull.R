`ahull` <-
function(x,y=NULL,alpha){

ashape.obj<-ashape(x,y,alpha)
compl<-complement(ashape.obj$delvor.obj,alpha=alpha)

# We consider the edges of the alpha-shape and also the alpha-extremes
ashape.edges<-matrix(ashape.obj$edges[,c("ind1","ind2")],ncol=2,byrow=FALSE)
noforget<-ashape.obj$alpha.extremes
num<-length(noforget)
# For each alpha-extreme we define the ball of radius 0 centered in the alpha-extreme to guarantee that isolated sample points are included in the alpha-hull
mat.noforget<-cbind(ashape.obj$x[noforget,],rep(0,num),rep(0,num),rep(0,num),rep(0,num))
# We determine the possible arcs of the alpha-hull: arcs.old
# If balls are coming for both sides of the edge the corresponding arcs can not define the boundary of the alpha-hull
ind2<-integer()
j<-0
nshape<-length(ashape.edges)*0.5
arcs<-matrix(0,nrow=nshape,ncol=6)

if (nshape>0){
for (i in 1:nshape){
ind<-which(ashape.edges[i,1]==compl[,"ind1"]&ashape.edges[i,2]==compl[,"ind2"])
if (length(ind)>0){
if(!((1<=sum((compl[ind,"ind"]==1)))&(sum((compl[ind,"ind"]==1))<length(ind)))){
which<-which(compl[ind,"r"]== min(compl[ind[compl[ind,"r"]>0],"r"]))
j<-j+1
arcs[j,]<-c(compl[ind[which],1],compl[ind[which],2],compl[ind[which],3],compl[ind[which],"v.x"],compl[ind[which],"v.y"],compl[ind[which],"theta"])
}
ind2<-c(ind2,ind)
}
}
}


arcs.old<-arcs[arcs[,3]>0,]
colnames(arcs.old)<-c("c1","c2","r","v.x","v.y","theta")

# We define the definitive arcs

arcs<-arcs.old
n.arc<-dim(arcs)[1]
watch<-1
nowatch<-0
j<-1

if (n.arc>0){
while (watch<=n.arc){
	ind.arc<-1:n.arc
	while(j<=n.arc){
      	if (j!=watch&(is.na(match(j,nowatch)))){
# The function inter calculates the intersection between two circles, given their radius and centers
			intersection<-inter(arcs[watch,1],arcs[watch,2],arcs[watch,3],arcs[j,1],arcs[j,2],arcs[j,3])
# If the balls intersect
			if (intersection$n.cut==2){
# Unitary vector defining the arc
				v.arc<-c(arcs[watch,"v.x"],arcs[watch,"v.y"])
# ang.OX is the angle of v.arc with the OX axis (defined in [0,2pi))
				if (v.arc[2]>=0){ang.OX<-acos(v.arc[1])}else {ang.OX<-2*pi-acos(v.arc[1])}
# We rotate v.arc (clockwise) an angle ang.OX, then
				v.arc.rot<-rotation(v.arc,ang.OX)
# Unitary vector defining the intersection
				v.int<-intersection$v1
# We rotate v.int (clockwise) an angle ang.OX, then
				v.int.rot<-rotation(v.int,ang.OX)
# We determine the relative position of v.int.rot in respect to v.arc.rot 
				if (v.int.rot[2]>=0){
					ang.v.int.rot.OX<-acos(v.int.rot[1])
					angles<-c(-arcs[watch,"theta"],arcs[watch,"theta"],ang.v.int.rot.OX-intersection$theta1,ang.v.int.rot.OX+intersection$theta1)
					names(angles)<-c("theta1","theta2","beta1","beta2")
					order<-names(sort(angles))

					theta1<- ang.OX-arcs[watch,"theta"]
					theta2<- ang.OX+arcs[watch,"theta"]
					beta1<- ang.v.int.rot.OX+ang.OX-intersection$theta1
					beta2<- ang.v.int.rot.OX+ang.OX+intersection$theta1

				if (all(order==c("theta1","theta2","beta1","beta2"))|all(order==c("beta1","beta2","theta1","theta2"))){
# The original arc remains complete
				}
				else if (all(order==c("theta1","beta1","theta2","beta2"))){
# theta1<beta1<theta2<beta2
# A piece of the original arc is removed 
					ang.middle<-(angles["beta1"]-angles["theta1"])/2
					v.new<-rotation(c(1,0),arcs[watch,"theta"]-ang.middle-ang.OX)
					arcs[watch,]<-c(arcs[watch,1],arcs[watch,2],arcs[watch,3],v.new[1],v.new[2],ang.middle)
				}
				else if (all(order==c("beta1","theta1","theta2","beta2"))){
# beta1<theta1<theta2<beta2
# The complete arc is removed
					nowatch<-c(nowatch,watch)
				}
				else  if (all(order==c("theta1","beta1","beta2","theta2"))){
# theta1<beta1<beta2<theta2
# The original arc turns into two arcs  (we add a row to the arcs matrix)
					ang.middle<-(angles["beta1"]-angles["theta1"])/2
					v.new<-rotation(c(1,0),arcs[watch,"theta"]-ang.middle-ang.OX)
					ang.middle2<-(angles["theta2"]-angles["beta2"])/2
					v.new2<-rotation(c(1,0),-arcs[watch,"theta"]+ang.middle2-ang.OX)

					arcs<-rbind(arcs,c(arcs[watch,1],arcs[watch,2],arcs[watch,3],v.new2[1],v.new2[2],ang.middle2))
					arcs[watch,]<-c(arcs[watch,1],arcs[watch,2],arcs[watch,3],v.new[1],v.new[2],ang.middle)
					n.arc<-n.arc+1
				}
				else if (all(order==c("beta1","theta1","beta2","theta2"))){
# beta1<theta1<beta2<theta2
# A piece of the original arc is removed
					ang.middle2<-(angles["theta2"]-angles["beta2"])/2
					v.new2<-rotation(c(1,0),-arcs[watch,"theta"]+ang.middle2-ang.OX)
					arcs[watch,]<-c(arcs[watch,1],arcs[watch,2],arcs[watch,3],v.new2[1],v.new2[2],ang.middle2)
				}
				} # close if v.int
###################################################
				else{
					ang.v.int.rot.OX<-acos(v.int.rot[1])
					angles<-c(-arcs[watch,"theta"],arcs[watch,"theta"],-ang.v.int.rot.OX-intersection$theta1,-ang.v.int.rot.OX+intersection$theta1)
					names(angles)<-c("theta1","theta2","beta1","beta2")
					order<-names(sort(angles))

					theta1<- ang.OX-arcs[watch,"theta"]
					theta2<- ang.OX+arcs[watch,"theta"]
					beta1<- -ang.v.int.rot.OX+ang.OX-intersection$theta1
					beta2<- -ang.v.int.rot.OX+ang.OX+intersection$theta1

				if (all(order==c("theta1","theta2","beta1","beta2"))|all(order==c("beta1","beta2","theta1","theta2"))){
# The original arc remains complete
				}
				else if (all(order==c("theta1","beta1","theta2","beta2"))){
# theta1<beta1<theta2<beta2
# A piece of the original arc is removed
					ang.middle<-(angles["beta1"]-angles["theta1"])/2
					v.new<-rotation(c(1,0),arcs[watch,"theta"]-ang.middle-ang.OX)
					arcs[watch,]<-c(arcs[watch,1],arcs[watch,2],arcs[watch,3],v.new[1],v.new[2],ang.middle)
				}
				else if (all(order==c("beta1","theta1","theta2","beta2"))){
# beta1<theta1<theta2<beta2
					nowatch<-c(nowatch,watch)
				}
				else  if (all(order==c("theta1","beta1","beta2","theta2"))){
# theta1<beta1<beta2<theta2
# The original arc turns into two arcs
					ang.middle<-(angles["beta1"]-angles["theta1"])/2
					v.new<-rotation(c(1,0),arcs[watch,"theta"]-ang.middle-ang.OX)							
					ang.middle2<-(angles["theta2"]-angles["beta2"])/2
					v.new2<-rotation(c(1,0),-arcs[watch,"theta"]+ang.middle2-ang.OX)

					arcs<-rbind(arcs,c(arcs[watch,1],arcs[watch,2],arcs[watch,3],v.new2[1],v.new2[2],ang.middle2))
					arcs[watch,]<-c(arcs[watch,1],arcs[watch,2],arcs[watch,3],v.new[1],v.new[2],ang.middle)
					n.arc<-n.arc+1
				}
				else if (all(order==c("beta1","theta1","beta2","theta2"))){
# beta1<theta1<beta2<theta2
# A piece of the original arc is removed
					ang.middle2<-(angles["theta2"]-angles["beta2"])/2
					v.new2<-rotation(c(1,0),-arcs[watch,"theta"]+ang.middle2-ang.OX)
					arcs[watch,]<-c(arcs[watch,1],arcs[watch,2],arcs[watch,3],v.new2[1],v.new2[2],ang.middle2)
				}
				}# close if v.int

			}# close if intersection
		} #close if j
		j<-j+1
	} #close while j
	watch<-watch+1
	j<-1
} # close while watch

ind<-1:dim(arcs)[1]
if (length(nowatch)>1){ind<-ind[-nowatch[-1]]}
ahull.arcs<-arcs[ind,]
length<-lengthahull(arcs[ind,])
ahull.arcs<-rbind(ahull.arcs,mat.noforget)
}
else{
ahull.arcs<-mat.noforget
length<-0}
ahull.obj<-list("arcs"=ahull.arcs,"length"=length,"complement"=compl,"alpha"=alpha,"ashape.obj"=ashape.obj,"x"=ashape.obj$x)
class(ahull.obj)<-"ahull"
invisible(ahull.obj)
}

