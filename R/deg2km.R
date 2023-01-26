"deg2km"<-function (data){
  earth<-6371
  curv<-2 * pi * earth / 360
  data$x<-data$Long*curv*cos(pi*data$Lat/180)
  data$y<-data$Lat*curv
  data
}
