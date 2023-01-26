
#Lee el archivo con datos de acustica
mgacu <- read.csv("JuvMgGam/Acu93-2006.csv")
#lee linea de costa 
costa <- read.csv("JuvMgGam/costa.csv")
names(mgacu)

#seleccionamos acada ano 
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

#seleccionamos cada ano
mgj97 <- mgjuv[mgjuv$Year==1997,,]
mgj99 <- mgjuv[mgjuv$Year==1999,,]
mgj00 <- mgjuv[mgjuv$Year==2000,,]
mgj01 <- mgjuv[mgjuv$Year==2001,,]
mgj02 <- mgjuv[mgjuv$Year==2002,,]
mgj04 <- mgjuv[mgjuv$Year==2004,,]
mgj05 <- mgjuv[mgjuv$Year==2005,,]
mgj06 <- mgjuv[mgjuv$Year==2006,,]

#GAM years vs Latitud
#Presencia de Juveniles
library(mgcv)
library(MASS)
library(shape)
library(Cardinal)

# presencia/ausencia de juveniles con los datos de los lances de pesca 
m0pres <- gam(Pjuv~s(Year,Lat)+s(Ppro),family=binomial,weights=frec,data=mgjuv)
summary(m0pres)
plot(m0pres,select=1)

#Grilla de prediccion
yrs <- seq(1997,2006,length=200)
lat <- seq(min(mgjuv$Lat),max(mgjuv$Lat),length=200)
#prof <-seq(min(mgjuv$Ppro),max(mgjuv$Ppro),length=200) 
grilla <- expand.grid(Year=yrs,Lat=lat,Ppro=50)
grilla$Presmgayi <- predict.gam(m0pres,grilla,type="response")
#grilla$Dmgayi <- predict.gam(m1,grilla,type="response")
#grilla$Eggs <- predict.gam(m1lgnsc,grilla,type="response")*grilla$Peggs
PrMgayi.srf <- list(x=yrs,y=lat,z=t(matrix(log(grilla$Presmgayi+0.1),ncol=200,byrow=TRUE)))
#DenMgayi.srf<- list(x=yrs,y=lat,z=t(matrix(grilla$Peggs,ncol=200,byrow=TRUE)))


op <- par(no.readonly=TRUE)
ps.options(horizontal = FALSE, bg="white",onefile = FALSE, paper = "special")
postscript("Figura_1.eps",height = 7, width = 10) #tamaño es en pulgadas
#par(oma=c(bottom,left,top,right),mar=c(bottom,left,top,right))
#jpeg("Figura_1.jpeg",width=960,height==500,quality=100)
x11()
par(oma=c(4,2,2,1),mar=c(1,3,2,1))
layout(matrix(1:2,nrow=1,ncol=2,byrow=T),widths=c(2,0.3))
image(PrMgayi.srf,col=gradient.colors,axes=F)
axis(2,col="black",cex.axis=1.2)
axis(1,col="black",cex.axis=1.2);box(lwd=1)
mtext(2,text="Latitude (?S)",line=3,cex=1.4)
mtext(1,text="Years",line=3,cex=1.4)
emptyplot(main="")
colorlegend(posx=c(0.2,0.3),posy=c(0.05,0.9),zlim=c(0,1),left=F,col=gradient.colors,digit = 1, dz = 0.2,main="")
par(op)



datj <- mgj00
dats <- mgacu00
#GAM Models for proporcion of juveniles
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
plot(dats$Lat,dats$Sa)
points(dats$Lat,dats$Sajuv,col=2)
plot(dats$Lat,dats$Pjuv)

## ANALISIS transformacion SAjuv
hist(dats$Sajuv)
ytransf <- dats$Sajuv^0.15
hist(ytransf)

boxplot(ytransf)
hist(ytransf,breaks=20,freq=F)
lines(density(ytransf))
shapiro.test(ytransf)
qqnorm(ytransf)
qqline(ytransf)

m2 <- gam(Sajuv^0.15~s(Lat)+s(Ppro),data=dats,family=gaussian(link="identity"))
gam.check(m2)
summary(m2)
plot(m2,select=1)
plot(m2,select=2)


dats$EstSajuv <- predict.gam(m2,dats,type="response")
dats$EstSajuv <- dats$EstSajuv^(1/0.15)
plot(dats$Lat,dats$Sajuv,cex=0.4)
lines(dats$Lat,dats$EstSajuv,type="l")
plot(dats$Lat,dats$EstSajuv,type="p")



#BUEN INTENTO...pero NO graciaslibrary(mixdist)
bins <- cut(datxx$Lat,br=seq(-41.2,-29,0.2),labels=seq(-41,-29,0.2),include.lowest=T)
frJuv <- table(datxx$Juv,bins)
barplot(frJuv,col="black")
frJuv
plot(bins,datxx$Juv)
#SOLO UNIVARIADOlibrary(mixtools)
attach(datxx)
#Estimacion de multiples normales
## DATA_SECTION
datxx <- list()
datxx$Lat <- dats$Lat
datxx$Juv <- dats$EstSajuv
datxx <- data.frame(datxx)
plot(datxx$Lat,datxx$Juv)
sum(datxx$Juv,na.rm=T)
#Parameters
 latmed <- c(-37.8,-35.5,-33,-30.5)
 latsd <- c(0.5,0.5,0.5,0.2)
 latp <- c(0.1,0.3,0.2,0.4)
 latp <- latp/sum(latp)
y1 <- dnorm(datxx$Lat,latmed[1],latsd[1])*latp[1]
y2 <- dnorm(datxx$Lat,latmed[2],latsd[2])*latp[2]
y3 <- dnorm(datxx$Lat,latmed[3],latsd[3])*latp[3]
y4 <- dnorm(datxx$Lat,latmed[4],latsd[4])*latp[4]
yhat <-(y1+y2+y3+y4)
plot(datxx$Lat,yhat)
lines(datxx$Lat,y1,col=2)
lines(datxx$Lat,y2,col=3)

yhat <-yhat*mean(datxx$Juv,na.rm=T)
lines(datxx$Lat,yhat,col=2,lwd=2)

#ANALISIS SAadultos
hist(dats$Saadul^0.15)
m3 <- gam(Saadul^0.15~s(Lat),data=dats,family=gaussian(link="identity"))
summary(m3)
plot(m3)
dats$EstSaadul <- predict.gam(m3,dats,type="response")
dats$EstSaadul <- dats$EstSaadul^(1/0.15)
plot(dats$Lat,dats$Saadul,cex=0.4)
points(dats$Lat,dats$EstSaadul,type="p",col=3)
points(dats$EstSaadul,dats$EstSajuv)
plot(dats$EstSaadul,dats$EstSajuv,pch="x",cex=0.3,col=2)

###### NO Va
m4 <- gam(Sajuv^0.15~s(Saadul),data=dats,family=gaussian(link="identity"))
gam.check(m4)
plot(m4,select=1,res=T,cex=1.2)
summary(m4)
m2.gl <- glm(Sajuv^0.15~Lat+Ppro+Saadul,data=dats)
summary(m2.gl)
par(mfrow=c(1,3))
termplot(m2.gl)
################




##### PARA EL CALCULO DE AREAS DE DISTRIBUCION DE JUVE > 0.5
#Convierte grados a km
#funcion de geofun: deg2km (al final del script)
source("JuvMgGam/R/deg2km.R")
dats <- deg2km(dats)
costa <- deg2km(costa)
names(dats)
plot(dats$x,dats$y,xlab="Easting (km)",ylab="Northing (km)",cex=0.4)
lines(costa$x,costa$y)

#Funcion de geofun: findlimits.fun (al final del script)
#Encuentra los limites del poligono de datos
require(spatstat)
require(sm)
require(maptools)
source("JuvMgGam/R/findlimits.fun.R")

spatstat.options(npixel=c(300,300))

survey.lim <- findlimits.fun(dats, dist=10, plot="limits")
## Lo primero es seleccionar con
## que limites nos vamos a quedar, en este caso, por ejemplo las areas
## 1 y 2

survey.lim <- survey.lim[c(1,2,4,5,6,11,13,14,15)]

## quitando el vertice final (duplicado)
#cuantos poligonos? npol=9
npol=9
for(i in 1:npol)
{
survey.lim[[i]]$x <- survey.lim[[i]]$x[-length(survey.lim[[i]]$x)]
survey.lim[[i]]$y <- survey.lim[[i]]$y[-length(survey.lim[[i]]$y)]
}

survey.lim.win <- owin(poly=survey.lim)
plot(survey.lim.win)
area.owin(survey.lim.win)

points(dats$x, dats$y, pch=".")
points(dats$x[dats$Pjuv>0.5],dats$y[dats$Pjuv>0.5], pch=19, cex=0.3, col="blue")

dats$Sea.area <- as.numeric(dirichletWeights(ppp(x=dats$x,y=dats$y, window=survey.lim.win)))
sum(dats$Sea.area)

### AREA POSTIVA DE JUVENILES#############
spatstat.options(npixel=c(300,300))
positive.Juv.lim <- findlimits.fun(dats[dats$Pjuv>0.5,], dist=10, plot="limits")
positive.Juv.lim <- positive.Juv.lim[c(1,2,3,4,5,7,8,9,10,12,13,14,15)]

## quitando el vertice final (duplicado)
#Cuantos poligonos? npol=13
npol=13
for(i in 1:npol){
positive.Juv.lim[[i]]$x <- positive.Juv.lim[[i]]$x[-length(positive.Juv.lim[[i]]$x)]
positive.Juv.lim[[i]]$y <- positive.Juv.lim[[i]]$y[-length(positive.Juv.lim[[i]]$y)]
}


#descuentos de los agujeros 8 que estan en la posicion 3(no se realizo, fue un ejemplo)

#positive.Juv.lim[[3]]$x <- rev(positive.Juv.lim[[3]]$x)
#positive.Juv.lim[[3]]$y <- rev(positive.Juv.lim[[3]]$y)

##realizarlos continuamente!!!
### los 4 poligonos se guardan en positive.lim.win
positive.Juv.lim.win <- owin(poly=positive.Juv.lim)
x11()
plot(positive.Juv.lim.win)

#plot(survey.lim.win)
points(dats$x, dats$y, pch=".")
points(dats$x[dats$Pjuv>0.5],dats$y[dats$Pjuv>0.5], pch=19, cex=0.8, col="blue")

## Se puede intentar ahora erosionar un poco el ·rea externa para evitar
## incluir demasiadas estaciones negativas en la parte de fuera
#positive.Juv.lim.win <- erode.owin(positive.Juv.lim.win, 1)

plot(positive.Juv.lim.win)
##plot(survey.lim.win)

points(dats$x, dats$y, pch=".")
points(dats$x[dats$Pjuv>0],dats$y[dats$Pjuv>0], pch=19, cex=0.6, col="blue")

## Finalmente se aÒade una columna para ver si las estaciones estan
## dentro o fuera del area positiva
dats$PositiveJuv <- inside.owin(dats$x, dats$y, positive.Juv.lim.win)



### AREA POSITIVA DE ADULTOS
spatstat.options(npixel=c(300,300))
positive.Adul.lim <- findlimits.fun(dats[dats$Padul>0.5,], dist=10, plot="limits")

#poligonos de 1999
positive.Adul.lim <- positive.Adul.lim[c(1,2,3,5,7,9,10,12,17,19,20,21,22)]
## quitando el vertice final (duplicado)
#cuantospoligonos? npol=?
#1999 npol=11
npol=13
for(i in 1:npol){
positive.Adul.lim[[i]]$x <- positive.Adul.lim[[i]]$x[-length(positive.Adul.lim[[i]]$x)]
positive.Adul.lim[[i]]$y <- positive.Adul.lim[[i]]$y[-length(positive.Adul.lim[[i]]$y)]
}

positive.Adul.lim.win <- owin(poly=positive.Adul.lim)
x11()
plot(positive.Adul.lim.win)

## Finalmente se aÒade una columna para ver si las estaciones estan
## dentro o fuera del area positiva
dats$PositiveAdul <- inside.owin(dats$x, dats$y, positive.Adul.lim.win)

###################
## Sumario de area muestreada y area positiva Juveniles
###################
Area99<-sum(dats$Sea.area) 
Area99juv<-sum(dats$Sea.area[dats$PositiveJuv])
Area99adul<-sum(dats$Sea.area[dats$PositiveAdul])

### HASTA AQUI INDICES DE COBERTURA ESPACAL Proporcion Juv y Adul > 0.5



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
modPpro <- gam(Ppro~s(Long,Lat),data=dats)
predictive.grid$Ppro <- predict.gam(modPpro,predictive.grid)

#Prediccion para el mapa de proporcion de Juveniles
predictive.grid$EstPjuv <- predict.gam(m1,predictive.grid,type="response")
predictive.grid$EstPjuv[!predictive.grid$Survey] <- 0

#Transformar los datos a una imagen
Pr.srf <- list(x=seq(min(predictive.grid$Long), max(predictive.grid$Long), length=300), y=seq(min(predictive.grid$Lat), max(predictive.grid$Lat), length=300), z=t(matrix(predictive.grid$EstPjuv, ncol = 300, byrow=TRUE)))

Pr.im <- im(mat=Pr.srf$z,xcol=Pr.srf$x, yrow=Pr.srf$y)

image(Pr.im, col=terrain.colors(100), xlab="Longitude W", ylab="Latitude S")
image(Pr.im, col=topo.colors(100), xlab="Longitude W", ylab="Latitude S",main="1999")

image(Pr.im, col=gray(seq(1,0.1,l=10)), xlab="Longitude W", ylab="Latitude S")
contour(Pr.im,nlevels=4,xlim=c(-75,-73),ylim=c(-42,-38))
lines(costa$Long+.1,costa$Lat,lwd=2,col="black")
points(dats$Long+.1,dats$Lat,pch=".")
points(datj$Long+.1,datj$Lat,pch="+",col="blue")

##Prediccion de la abundancia
predictive.grid$EstSajuv <- predict.gam(m2,predictive.grid,type="response")
predictive.grid$EstSajuv <-predictive.grid$EstSajuv^(1/0.15) 
predictive.grid$EstSajuv[!predictive.grid$Survey] <- 0

#Transformar los datos a una imagen
Pr.srf <- list(x=seq(min(predictive.grid$Long), max(predictive.grid$Long), length=300), y=seq(min(predictive.grid$Lat), max(predictive.grid$Lat), length=300), z=t(matrix(predictive.grid$EstSajuv, ncol = 300, byrow=TRUE)))

Pr.im <- im(mat=Pr.srf$z,xcol=Pr.srf$x, yrow=Pr.srf$y)

image(Pr.im, col=terrain.colors(100), xlab="Longitude W", ylab="Latitude S")
image(Pr.im, col=topo.colors(100), xlab="Longitude W", ylab="Latitude S",main="1999")

image(Pr.im, col=gray(seq(1,0.1,l=10)), xlab="Longitude W", ylab="Latitude S")
lines(costa$Long+.1,costa$Lat,lwd=2,col="black")




mgj99 <- datj
mgacu99 <- dats
#write.table(mgacu99,"mgresults99.csv")
#write.table(predictive.grid,"PredGrilla99.csv")





#ANALISIS DE CAMBIOS EN LA DENSIDAD DE JUVENILES Y ADULTOS
#setwd("/Users/luiscubillos/Rwork/MGayi/JuvMgGam/Todos")
setwd("~/Documents/Tesis Pregrado/Todos")
library(lattice)
dir()
mg97 <- read.csv("Todos/mgac1997.csv")
mg99 <- read.csv("Todos/mgac1999.csv")
mg00 <- read.csv("Todos/mgac2000.csv")
mg01 <- read.csv("Todos/mgac2001.csv")
mg02 <- read.csv("Todos/mgac2002.csv")
mg04 <- read.csv("Todos/mgac2004.csv")
mg05 <- read.csv("Todos/mgac2005.csv")
mg06 <- read.csv("Todos/mgac2006.csv")

mgac <- rbind(mg97,mg99,mg00,mg01,mg02,mg04,mg05,mg06)

xyplot(Sajuv~Saadul|as.factor(Year),data=mgac,xlab="Densidad Adultos (t/mn2)",ylab="Densidad de Juveniles (t/mn2)")


m0 <- gam(Sajuv~s(Lat)+s(Ppro),data=mgac)
summary(m0)
m1 <- gam(Saadul~s(Lat)+s(Ppro),data=mgac)
summary(m1)
par(mfrow=c(2,1),mar=c(5,5,1,1))
plot(m0,select=1,xlab="Latitud Sur",ylab="Densidad Juveniles",shade=T,shade.col="gray80",)
plot(m1,select=1,xlab="Latitud Sur",ylab="Densidad Adultos",shade=T,shade.col="gray80",)
text(-42,500,"1997-2006",pos=4)

par(mfrow=c(2,1),mar=c(5,5,1,1))
plot(m0,select=2,xlab="Profundidad (m)",ylab="Densidad Juveniles",shade=T,shade.col="gray80",)
plot(m1,select=2,xlab="Profundidad (m)",ylab="Densidad Adultos",shade=T,shade.col="gray80",)
text(10,500,"1997-2006",pos=4)

## DENSIDAD MEDIA POR AÑO
#JUVENILES
m3 <- gam(Sajuv^0.15~as.factor(Year)+s(Long,Lat)+s(Ppro),data=mgac)
summary(m3)
plot.gam(m3,select=2)
#Extraccion de la densidad anual (coeficientes)
djuv <- coef(m3)[c(0,1,2,3,4,5,6,7,8)]
dj <- rep(NA,8)
for(i in 1:7){
  djuv[i+1]=djuv[i+1]+djuv[1]
}
dj <- djuv^(1/0.15)
dj

#Adultos
m4 <- gam(Saadul^0.15~as.factor(Year)+s(Long,Lat)+s(Ppro),data=mgac)
summary(m4)
plot.gam(m4,select=2)
#Extraccion de la densidad anual (coeficientes)
dadul <- coef(m4)[c(0,1,2,3,4,5,6,7,8)]
da <- rep(NA,8)
for(i in 1:7){
  dadul[i+1]=dadul[i+1]+dadul[1]
}
da <- dadul^(1/0.15)
da

yr <- c(1997,1999,2000,2001,2002,2004,2005,2006)
plot(yr,da,pch=19,xlab="",ylab="Densidad (t/mn2)")
lines(yr,da)
points(yr,dj,col=2,pch=19)
lines(yr,dj,lty=2,col=2)
legend(2004,150,legend=c("Adultos","Juveniles"),col=c(1,2),lty=c(1,2),pch=c(19,19))

#Relacion adultos-juveniles
plot(dj~da,col=2,cex=1.2,pch=19,xlab="Densidad adultos (t/mn2)",ylab="Densidad juveniles (t/mn2)",type="b",lty=2)
m5 <- lm(log(dj/da)~da)
summary(m5)
x <- seq(0,200,1)
y <- x*exp(coef(m5)[1]+coef(m5)[2]*x)
lines(x,y,lwd=2)
text(5,17,"2006",pos=3)
text(15,1,"2002",pos=4)
text(18,26,"2001",pos=4)
text(100,25,expression(J==A*exp(0.448-0.029*A)),pos=4)
text(100,22,expression(r^2 == 0.672),pos=4)


