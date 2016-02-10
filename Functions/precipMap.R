library(zonalDaymet)
library(maptools)
library(maps)
library(raster)
library(sp)
library(akima)
library(plotrix)
source("C:/Users/Evan/Desktop/Suckers/Figures/filled.contour2.R")#a more flexible filled.contour

#densityColors<-heat.colors(ceiling(max(density[,,1]))+1)

load("~/lee/dataStore/cleanData/niles.rdata")
rm(fish)
rm(skippedPasses)
rm(nonTrout)
siteData<-siteData[,.(lat=mean(siteLatitude),
                      long=mean(siteLongitude)),
                   by=site]
###Prepare data for big map
states<-c("Pennsylvania",
          "West Virginia",
          "Virginia",
          "Maryland",
          "Delaware",
          "New York",
          "New Jersey",
          "Massachusetts",
          "Connecticut",
          "Vermont",
          "New Hampshire",
          "Maine",
          "Rhode Island")

precip<-list()
for(s in states){
  precip[[s]]<-readRDS(paste0("~/lee/figures/stateData/",s,".rds"))
}
precip<-do.call(rbind,precip)


flowLines<-readShapeLines("~/shapefiles/NHDPlusV21_MA_02_NHDSnapshot_04/Hydrography/NHDFlowline.shp")
mainStem<-flowLines[grep("oyal",flowLines$GNIS_NAME),]
huc<-flowLines[grep("0205020600",flowLines$REACHCODE),]
huc<-crop(huc,extent( -77,-76.22213,41.23473,41.64562))
naId<-huc[is.na(huc$GNIS_ID),]
#these ids are in the watershed (indices on unique(huc$GNIS_ID))
ids<-c(3,4,5,6,7,10,11,12,14,16,17,18,19,20,21,23,24,25,29,35,37,38,
       39,40,43,46,70,71,72,73,74,75,76,77,78,79,
       80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97:148,
       183)
#these reach codes (as indices on unique(naId$REACHCODE)) for id==NA are in the watershed
reaches<-c(1,2,4,5,7,8,9,13,14,15,54,55,56,57,58,59,60,61,62,63,64,
           65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,
           84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,
           102,103,104,105,106,107,108,109,110,111,112,113,114,115,
           116,117,118,119,120,121,122,123,124,125,126,127,128,129,
           130,131,132,133,134,135,136,137,138,139,140,141,142,143,
           144,145,146,147,148,149,150,151,152,153,154,155,156,157,
           158,159,160,161,162,163,164,165,166,167,168,169,170,171,
           172,173,174,175,176,177,178,179,180,181,182,183,184,185,
           186,187)

loyalsock<-huc[huc$GNIS_ID %in% unique(huc$GNIS_ID)[ids]|
                 huc$REACHCODE %in% unique(naId$REACHCODE)[reaches],]

# Code for finding which are in the watershed
# for(i in reaches){
#   plot(boof)
#   plot(def,add=T,col='blue',lwd=2)
#   plot(huc[huc$REACHCODE==unique(naId$REACHCODE)[i],],col='red',lwd=6,add=T)
#   usr <- scan(what=character(),nmax=1,quiet=TRUE)
#   if(usr){loyal<-c(loyal,i)} else {print(i)}
# }

# flowLines<-crop(flowLines,with(limits,extent(xmin,xmax,ymin,ymax)))


####Prepare data for zoomed in map
pad=0.1

usgsGage<-data.table(long=-76.91278,
                     lat=41.325)

limits<-as.matrix(extent(loyalsock))

limits<-list(xmin=limits[1,1]-pad,
             xmax=limits[1,2]+pad,
             ymin=limits[2,1]-pad,
             ymax=limits[2,2]+pad)

sitePrecip<-precip[lat>limits$ymin&
                     lat<limits$ymax&
                     long>limits$xmin&
                     long<limits$xmax
                   ]
res<-0.001
sitePrecip<-interp(x=sitePrecip$long,y=sitePrecip$lat,z=sitePrecip$precip,
                   xo=seq(min(sitePrecip$long),max(sitePrecip$long),by=res),
                   yo=seq(min(sitePrecip$lat),max(sitePrecip$lat),by=res))


colors<-heat.colors(max(precip$precip+1))
colors<-colorRampPalette(c('black','white'))(max(precip$precip+1))

width<-6.8
r2<-(limits$xmax-limits$xmin)/(limits$ymax-limits$ymin)
r1<-diff(range(precip$long))/diff(range(precip$lat))
r<-r1/r2
x2<-width/(r+1)
x1<-width-x2
height<-x1/r1

tiff.par("~/lee/figures/precipMap.tif",width=width,height=height)
layout(matrix(c(1,2),nrow=c(1)),
       widths=c(x1,x2))
par(mar=c(0,0,0,0))
plot(lat~long,data=precip,col=colors[precip+1],pch=15,cex=0.03,
     axes=F,xlab="",ylab="")
map("state","Pennsylvania",add=T,col=gray(0.7),mar=c(0,0,0,0))
pad<-0.05
rect(limits$xmin+pad,limits$ymin+pad,
     limits$xmax-pad,limits$ymax-pad)
legendLabs<-0:max(precip$precip)
legendLabs[!legendLabs %in% c(0,100,200,300)]<-NA
color.legend(min(precip$long)+3,
             max(precip$lat)-3.5,
             min(precip$long)+3.5,
             max(precip$lat)-0.5,
             legend=legendLabs,cex=0.7,
             rect.col=colors,
             gradient="y")
text(min(precip$long)+3.25,max(precip$lat)-0.1,"Precipitation (mm)",cex=0.8)
marPad<-0.1

par(mar=c(x1*marPad,height*marPad,x1*marPad,height*marPad))
filled.contour2(sitePrecip$x,sitePrecip$y,sitePrecip$z,
                nlevels=327,
                #levels=0:max(precip$precip),
                col=colors,
                xlim=c(limits$xmin+pad,limits$xmax-pad),
                ylim=c(limits$ymin+pad,limits$ymax-pad),
                zlim=c(0,326),
                axes=F)
plot(loyalsock,col=gray(0.3),add=T)
plot(mainStem,col=gray(0.3),add=T,lwd=2)
points(siteData$long,siteData$lat,pch=19,cex=0.5)#,col=densityColors[round(density[,2,1]+1)])
points(-76.91278,41.325,pch="*",cex=2)
box(bty='o',lwd=1)
par(xpd=NA)
# arrows(limits$xmin,limits$ymax,limits$xmin-0.231,limits$ymax-0.167,length=0)
# arrows(limits$xmin,limits$ymin,limits$xmin-0.231,limits$ymin+0.136,length=0)
dev.off()



