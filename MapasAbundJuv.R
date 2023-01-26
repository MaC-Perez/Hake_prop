# MAPAS DE DISTRIBUCION DE JUVENILES (Presencia-Ausencia)
#MAPAS DE DISTRIBUCION DE LA ABUNDANCIA DE JUVENILES
#funcion de geofun: deg2km (al final del script)
source("JuvMgGam/R/deg2km.R")
require(spatstat)
require(sm)
require(maptools)
source("JuvMgGam/R/findlimits.fun.R")
library(mgcv)
library(MASS)
library(PBSmapping)
data(worldLLhigh)
source("JuvMgGam/R/source_indicators.r")
# Arreglos
worldLLhigh2 <- worldLLhigh
worldLLhigh2$X <- worldLLhigh2$X-360
clr <- PBSclr()
names(clr)


#Lee el archivo con se?al acustica
mgacu <- read.csv("JuvMgGam/Acu93-2006.csv")
costa <- read.csv("JuvMgGam/costa.csv")
names(mgacu)

mgacu97 = mgacu[mgacu$Year==1997,,]
mgacu99 = mgacu[mgacu$Year==1999,,]
mgacu00 = mgacu[mgacu$Year==2000,,]
mgacu01 = mgacu[mgacu$Year==2001,,]
mgacu02 = mgacu[mgacu$Year==2002,,]
mgacu04 = mgacu[mgacu$Year==2004,,]
mgacu05 = mgacu[mgacu$Year==2005,,]
mgacu06 = mgacu[mgacu$Year==2006,,]

#Lee el archivo de lances de pesca con proporcion de juveniles
mgjuv <- read.csv("JuvMgGam/Pjuv19972006.csv")

mgj97 <- mgjuv[mgjuv$Year==1997,,]
mgj99 <- mgjuv[mgjuv$Year==1999,,]
mgj00 <- mgjuv[mgjuv$Year==2000,,]
mgj01 <- mgjuv[mgjuv$Year==2001,,]
mgj02 <- mgjuv[mgjuv$Year==2002,,]
mgj04 <- mgjuv[mgjuv$Year==2004,,]
mgj05 <- mgjuv[mgjuv$Year==2005,,]
mgj06 <- mgjuv[mgjuv$Year==2006,,]

yr <- "2006"
datj <- mgj06
dats <- mgacu06

#GAM Models for proporcionof juveniles
m1 <- gam(Pjuv~ s(Long,Lat)+s(Ppro),weights=frec,family=binomial,data=datj)
summary(m1)
plot(m1,select=2)
points(datj$Long,datj$Lat,col=3)
points(dats$Long,dats$Lat,col=5)

plot(m1,select=2)

gam.check(m1)
dats$Pjuv <- predict.gam(m1,dats,type="response")
dats$Padul <- 1-dats$Pjuv
dats$Sajuv <- dats$Sa*dats$Pjuv
dats$Saadul <- dats$Sa*dats$Padul
plot(dats$Lat,dats$Saadul)
points(dats$Lat,dats$Sajuv,col=2)
plot(dats$Lat,dats$Pjuv)

## ANALISIS transformacion SAjuv
ytransf <- dats$Sajuv^0.15
boxplot(ytransf)
hist(ytransf,breaks=20,freq=F)
lines(density(ytransf))
shapiro.test(ytransf)
qqnorm(ytransf)
qqline(ytransf)

m2 <- gam(Sajuv^0.15~s(Long,Lat)+s(Ppro),data=dats,family=gaussian(link="identity"))
gam.check(m2)
summary(m2)
plot(m2,select=1)
plot(m2,select=2)
m3 <- gam(Saadul^0.15~s(Long,Lat)+s(Ppro),data=dats,family=gaussian(link="identity"))
summary(m3)
plot(m3,select=2)

#LIMITES DEL SURVEY
#Convierte grados a km
dats <- deg2km(dats)
costa <- deg2km(costa)
names(dats)
plot(dats$x,dats$y,xlab="Easting (km)",ylab="Northing (km)",cex=0.4)
lines(costa$x,costa$y)

#Funcion de geofun: findlimits.fun (al final del script)
#Encuentra los limites del poligono de datos
spatstat.options(npixel=c(300,300))
survey.lim <- findlimits.fun(dats, dist=10, plot="limits")

#1997: survey.lim <- findlimits.fun(dats, dist=9, plot="limits")
#1999: survey.lim <- findlimits.fun(dats, dist=10, plot="limits")
#2000: survey.lim <- findlimits.fun(dats, dist=10, plot="limits")
#2001: survey.lim <- findlimits.fun(dats, dist=9, plot="limits")
#2002: survey.lim <- findlimits.fun(dats, dist=10, plot="limits")
#2004: survey.lim <- findlimits.fun(dats, dist=10, plot="limits")
#2005: survey.lim <- findlimits.fun(dats, dist=10, plot="limits")
#2006: survey.lim <- findlimits.fun(dats, dist=10, plot="limits")

## Lo primero es seleccionar con
## que limites nos vamos a quedar, en este caso, por ejemplo las ·reas
## que no esten contenidas dentro de otras areas

survey.lim <- survey.lim[c(1,3,4,5,8,9,10,11)]

#1997: survey.lim <- survey.lim[c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,17,18)]
#1999: survey.lim <- survey.lim[c(1,2,3,4,5,7,10,11,13)]
#2000: survey.lim <- survey.lim[c(1,2,4,5,6,11,13,14,15)]
#2001: survey.lim <- survey.lim[c(1,2,3,4,5,6,7,9,10,13)]
#2002: survey.lim <- survey.lim[c(1,2,3,4,5,6,8,9,10)]
#2004: survey.lim <- survey.lim[c(1,4,6,8)]
#2005: survey.lim <- survey.lim[c(1,2,7,8)]
#2006: survey.lim <- survey.lim[c(1,3,4,5,8,9,10,11)]

## quitando el vertice final (duplicado)
#cuantos poligonos? npol=?
npol=8
for(i in 1:npol)
{
survey.lim[[i]]$x <- survey.lim[[i]]$x[-length(survey.lim[[i]]$x)]
survey.lim[[i]]$y <- survey.lim[[i]]$y[-length(survey.lim[[i]]$y)]
}

survey.lim.win <- owin(poly=survey.lim)
plot(survey.lim.win)
area.owin(survey.lim.win)
dats$Sea.area <- as.numeric(dirichletWeights(ppp(x=dats$x,y=dats$y, window=survey.lim.win)))
sum(dats$Sea.area)

#AHORA SE MAPEA
#Grilla predictiva
predictive.grid <- expand.grid(Lat = seq(min(dats$Lat)-0.1, max(dats$Lat)+0.1, length=300),Long = seq(min(dats$Long)-0.1, max(dats$Long)+0.1, length=300))
predictive.grid <- deg2km(predictive.grid)
###Se recorta el area del survey
predictive.grid$Survey <- inside.owin(predictive.grid$x, predictive.grid$y, survey.lim.win)
#Graficamos para verificar
plot(survey.lim.win)
points(predictive.grid$x, predictive.grid$y, pch=".")
points(predictive.grid$x[predictive.grid$Survey], predictive.grid$y[predictive.grid$Survey], pch=".", col="red")
#Preddicciones de Pjuv en la grilla
#modelo para Ppro
modPpro <- gam(Ppro~s(Long,Lat),data=dats);summary(modPpro)
predictive.grid$Ppro <- predict.gam(modPpro,predictive.grid)

#Prediccion para el mapa de abundancia de Juveniles
predictive.grid$EstSajuv <- predict.gam(m2,predictive.grid,type="response")
predictive.grid$EstSajuv <- predictive.grid$EstSajuv^(1/0.15)
predictive.grid$EstSajuv[!predictive.grid$Survey] <- NA
predictive.grid$EstSaadul <- predict.gam(m3,predictive.grid,type="response")
predictive.grid$EstSaadul <- predictive.grid$EstSaadul^(1/0.15)
predictive.grid$EstSaadul[!predictive.grid$Survey] <- NA
predictive.grid$EstPjuv <- predict.gam(m1,predictive.grid,type="response")
predictive.grid$EstPjuv[!predictive.grid$Survey] <- NA
(paste("predgrilla",yr,".csv",sep=""))
#write.csv(predictive.grid,paste("predgrilla",yr,".csv",sep=""))

#Transformar los datos a una imagen
PrJuv.srf <- list(x=seq(min(predictive.grid$Long), max(predictive.grid$Long), length=300), y=seq(min(predictive.grid$Lat), max(predictive.grid$Lat), length=300), z=t(matrix(predictive.grid$EstSajuv, ncol = 300, byrow=TRUE)))
PrJuv.im <- im(mat=PrJuv.srf$z,xcol=PrJuv.srf$x, yrow=PrJuv.srf$y)

PrAdul.srf <- list(x=seq(min(predictive.grid$Long), max(predictive.grid$Long), length=300), y=seq(min(predictive.grid$Lat), max(predictive.grid$Lat), length=300), z=t(matrix(predictive.grid$EstSaadul, ncol = 300, byrow=TRUE)))
PrAdul.im <- im(mat=PrAdul.srf$z,xcol=PrAdul.srf$x, yrow=PrAdul.srf$y)

PrPjuv.srf <- list(x=seq(min(predictive.grid$Long), max(predictive.grid$Long), length=300), y=seq(min(predictive.grid$Lat), max(predictive.grid$Lat), length=300), z=t(matrix(predictive.grid$EstPjuv, ncol = 300, byrow=TRUE)))
PrPjuv.im <- im(mat=PrPjuv.srf$z,xcol=PrPjuv.srf$x, yrow=PrPjuv.srf$y)



image(PrJuv.im, col=terrain.colors(100), xlab="Longitude W", ylab="Latitude S")
image(PrJuv.im, col=topo.colors(100), xlab="Longitude W", ylab="Latitude S",main=yr)
image(Pr.im, col=gray(seq(1,0.1,l=100)), xlab="Longitude W", ylab="Latitude S")

image(PrAdul.im, col=terrain.colors(100), xlab="Longitude W", ylab="Latitude S")
image(PrAdul.im, col=topo.colors(100), xlab="Longitude W", ylab="Latitude S",main=yr)
image(PrAdul.im, col=gray(seq(1,0.1,l=100)), xlab="Longitude W", ylab="Latitude S")

image(PrPjuv.im, col=terrain.colors(100), xlab="Longitude W", ylab="Latitude S")
image(PrPjuv.im, col=topo.colors(100), xlab="Longitude W", ylab="Latitude S",main=yr)
image(PrPjuv.im, col=gray(seq(1,0.1,l=100)), xlab="Longitude W", ylab="Latitude S")


#contour(Pr.im,nlevels=4,xlim=c(-75,-73),ylim=c(-42,-38))
lines(costa$Long+.1,costa$Lat,lwd=2,col="black")
points(dats$Long+.1,dats$Lat,pch=".",col=2)

#Alternativa

par(mar=c(5,5,0.5,0.5)+0.5)
plotMap(worldLLhigh2,col="brown",bg="darkblue",ylim=c(-42,-29),xlim=c(-76,-69),projection="LL",axes=F,xlab="",ylab="",cex.axis=1.2)
axis(2,yaxp=c(-42,-29,5),las=1)
axis(1,xaxp=c(-76,-69,5),las=1)
axis(3,xaxp=c(-76,-69,5),las=1,tick = T,labels = F)
mtext(side=2,line=3.5,"Latitud Sur",cex=1.2)
mtext(side=1,line=1.8,"Longitud Oeste",cex=1)
image(PrJuv.im, col=topo.colors(100),add=TRUE,axes=FALSE)
lines(costa$Long,costa$Lat, pch=".",cex=1.5,col="black")
text(-71.34,-29.953,"Coquimbo",pos=4,cex=0.8,col="white")
text(-71.62,-33.039,"Valparaiso",pos=4,cex=0.8,col="white")
text(-72.45,-35.39,"Constitucion",pos=4,cex=0.7,col="white")
text(-73,-36.55,"Dichato",pos=4,cex=0.7,col="white")
text(-73.25,-37.2,"Golfo de Arauco",pos=4,cex=0.7,col="white")
text(-73.75,-37.473,"Pta. Lavapie",pos=4,cex=0.7,col="white")
text(-73.7,-37.78,"Lebu",pos=4,cex=0.7,col="white")
text(-75.3,-38.3,"Isla Mocha",pos=4,cex=0.7,col="white")
text(-73.4,-39.81798931,"Corral",pos=4,cex=0.8,col="white")
text(-76,-29.5,"Densidad",pos=4,cex=1,col="white")
text(-76,-30,"Juveniles",pos=4,cex=1,col="white")
text(-76,-30.5,yr,pos=4,cex=1,col="white")

#image(Pr.im, col=topo.colors(100),axes=F)


par(mar=c(5,5,0.5,0.5)+0.5)
plotMap(worldLLhigh2,col="brown",bg="darkblue",ylim=c(-42,-29),xlim=c(-76,-69),projection="LL",axes=F,xlab="",ylab="",cex.axis=1.2)
axis(2,yaxp=c(-42,-29,5),las=1)
axis(1,xaxp=c(-76,-69,5),las=1)
axis(3,xaxp=c(-76,-69,5),las=1,tick = T,labels = F)
mtext(side=2,line=3.5,"Latitud Sur",cex=1.2)
mtext(side=1,line=1.8,"Longitud Oeste",cex=1)
image(PrPjuv.im, col=topo.colors(100),add=TRUE,axes=FALSE)
lines(costa$Long,costa$Lat, pch=".",cex=1.5,col="black")
text(-71.34,-29.953,"Coquimbo",pos=4,cex=0.8,col="white")
text(-71.62,-33.039,"Valparaiso",pos=4,cex=0.8,col="white")
text(-72.45,-35.39,"Constitucion",pos=4,cex=0.7,col="white")
text(-73,-36.55,"Dichato",pos=4,cex=0.7,col="white")
text(-73.25,-37.2,"Golfo de Arauco",pos=4,cex=0.7,col="white")
text(-73.75,-37.473,"Pta. Lavapie",pos=4,cex=0.7,col="white")
text(-73.7,-37.78,"Lebu",pos=4,cex=0.7,col="white")
text(-75.3,-38.3,"Isla Mocha",pos=4,cex=0.7,col="white")
text(-73.4,-39.81798931,"Corral",pos=4,cex=0.7,col="white")
text(-76,-29.5,"Presencia",pos=4,cex=1,col="white")
text(-76,-30,"Juveniles",pos=4,cex=1,col="white")
text(-76,-30.5,yr,pos=4,cex=1,col="white")

par(mar=c(5,5,0.5,0.5)+0.5)
plotMap(worldLLhigh2,col="brown",bg="darkblue",ylim=c(-42,-29),xlim=c(-76,-69),projection="LL",axes=F,xlab="",ylab="",cex.axis=1.2)
axis(2,yaxp=c(-42,-29,5),las=1)
axis(1,xaxp=c(-76,-69,5),las=1)
axis(3,xaxp=c(-76,-69,5),las=1,tick = T,labels = F)
mtext(side=2,line=3.5,"Latitude (ºS)",cex=1.2)
mtext(side=1,line=1.8,"Longitude (ºW)",cex=1)
image(PrAdul.im, col=topo.colors(100),add=TRUE,axes=FALSE)
lines(costa$Long,costa$Lat, pch=".",cex=1.5,col="black")
text(-71.34,-29.953,"Coquimbo",pos=4,cex=0.8,col="white")
text(-71.62,-33.039,"Valparaíso",pos=4,cex=0.8,col="white")
text(-72.45,-35.39,"Constitución",pos=4,cex=0.7,col="white")
text(-72.8,-36,"Punta Nugurne",pos=4,cex=0.7,col="white")
text(-73,-36.55,"Dichato",pos=4,cex=0.7,col="white")
text(-73.25,-37.2,"Golfo de Arauco",pos=4,cex=0.7,col="white")
text(-73.75,-37.473,"Pta. Lavapie",pos=4,cex=0.7,col="white")
text(-73.7,-37.78,"Lebu",pos=4,cex=0.7,col="white")
text(-75.3,-38.3,"Isla Mocha",pos=4,cex=0.7,col="white")
text(-73.4,-39.81798931,"Corral",pos=4,cex=0.7,col="white")
text(-75,-29.5,"Presencia",pos=4,cex=1,col="white")
text(-75,-30,"Adultos",pos=4,cex=1,col="white")
text(-75,-30.5,yr,pos=4,cex=1,col="white")

