`lengthahull` <-
function(ahull.arcs){
# theta= angle of each arch of the alpha-hull
gamma<-2*ahull.arcs[,"theta"]
# Length of an arc
# Length of an arc of a circle with radius r y angle gamma
# l=(2*pi*r)/(2*pi/gamma)=gamma*r
length<-gamma*ahull.arcs[,"r"]
return(length=sum(length))
}

