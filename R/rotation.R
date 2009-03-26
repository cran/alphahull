`rotation` <-
function(v,theta){
# Clockwise rotation
# theta.... angle of rotation
# v=(v1,v2)

v.rot<-numeric(2)
v.rot[1]<-cos(theta)*v[1]+sin(theta)*v[2]
v.rot[2]<--sin(theta)*v[1]+cos(theta)*v[2]

return(v.rot)
}

