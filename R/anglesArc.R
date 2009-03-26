`anglesArc` <-
function(v,theta){

# Angle of v in respect to OX (in [0,2*pi))
theta.OX<-ifelse(v[2]>=0,acos(v[1]),2*pi-acos(v[1]))
angs<-c(theta.OX-theta,theta.OX+theta)
return(angs)

}

