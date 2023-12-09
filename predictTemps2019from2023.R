#R
#GJR 11/29/23
#Fit linear model predicting tree phloem temperatures from snotel stations (2023) and predict phloem temperature for 2019


#read ibutton data
data<-read.table('GuanellaNaylorRoadIbuttonFall2023.csv',header=T,row.names=NULL,sep=",",stringsAsFactors=F,quote="\"",comment.char="",skip=14)
names(data)<-c("date","time","unit","temperature")
naylor<-data

data<-read.table('GuanellaCampgroundIbuttonFall2023.csv',header=T,row.names=NULL,sep=",",stringsAsFactors=F,quote="\"",comment.char="",skip=14)
names(data)<-c("date","time","unit","temperature")
camp<-data

#cut out transport days for the ibuttons, then calculate daily min and max
getMinMax<-function(data) {
  
  #this is hacky, but only way I could get the posix objects to aggregate into a vector
  dateTime<-as.POSIXct(strptime( paste(data$date[1],data$time[1]) , "%m/%d/%Y %I:%M:%S %p"))
  for (i in 2:nrow(data)) {
    dateTime<-c(dateTime,  as.POSIXct(strptime( paste(data$date[i],data$time[i]) , "%m/%d/%Y %I:%M:%S %p")) )
  }
  data$dateTime<-dateTime
  data$doy<-as.POSIXlt(data$dateTime, format = "%d.%m")$yday+1 #Posix counts Jan1 as '0'
  
  #drop data prior to 8/6 and after 10/15 
  data<-data[data$doy >= 217 & data$doy <= 287,]

  
  #calculate daily minima and maxima
  minVec<-numeric()
  maxVec<-numeric()
  doyVec<-unique(data$doy)
  for(doy in doyVec) {
    minVec<-c(minVec,min(data$temperature[data$doy==doy]))
    maxVec<-c(maxVec,max(data$temperature[data$doy==doy]))
  }
  minMax<-data.frame(doy=doyVec,min=minVec,max=maxVec)
  out<-list(data,minMax)
  names(out)<-c("daily","minMax")
  return(out)
  
}

naylor<-getMinMax(naylor)
camp<-getMinMax(camp)



#plot time series
library(ggplot2)
ggplot(aes(x = dateTime, y = temperature), data = naylor$daily) + geom_line()
ggplot(aes(x = dateTime, y = temperature), data = camp$daily) + geom_line()


#fit linear models to 2023 data, make predictions for 2019 data
snotel<-read.csv('SNOTEL_EchoLakeAndJackwackerGulch.csv',skip=1,header=T)
data2023<-merge(snotel[snotel$Year==2023,],naylor$minMax)
names(data2023)[8:9]<-c("naylorMin","naylorMax")
data2023<-merge(data2023,camp$minMax)
names(data2023)[10:11]<-c("campMin","campMax")

fitPred<-function(response,predictors,data,min=TRUE) {
  f <- paste(paste(response), "~", paste(predictors, collapse=" + "))
  fm<-lm(f,data=data)
  predTemp<-predict(fm,newdata=snotel[snotel$Year==2019,])
  meanTemp<-mean(predTemp)
  CI<-c(meanTemp-2*(sd(predTemp)/sqrt(length(predTemp))),meanTemp+2*(sd(predTemp)/sqrt(length(predTemp))))
  out<-c(meanTemp,CI)
  if (min) {out<-c(out,sum(predTemp <= 5))} else {
    out<-c(out,sum(predTemp >= 10))
  }
  return(out)
}

#Predictions for 5 August to 14 October

#Naylor min, mean=4.668718, CI=3.631923  5.705513 31 days with min<=5
fitPred('naylorMin',c("Echo_min","Jack_min"),data2023,min=T)

#Naylor max, mean=12.24493, CI=11.29012 13.19975, 57 days with max<=10
fitPred('naylorMax',c("Echo_max","Jack_max"),data2023,min=F)

#Camp min, mean=4.373753, CI=3.317876  5.429629  33 days with min<=5
fitPred('campMin',c("Echo_min","Jack_min"),data2023,min=T)

#Camp max, mean=14.26872, CI=13.24017 15.29727, 63 days with max>=10
fitPred('campMax',c("Echo_max","Jack_max"),data2023,min=F)


#Predictions for 5 August to 27 September (first date of beetle collection in 2019)

fitPred<-function(response,predictors,data,min=TRUE) {
  f <- paste(paste(response), "~", paste(predictors, collapse=" + "))
  fm<-lm(f,data=data)
  predTemp<-predict(fm,newdata=snotel[snotel$Year==2019 & snotel$doy <=270,])
  meanTemp<-mean(predTemp)
  CI<-c(meanTemp-2*(sd(predTemp)/sqrt(length(predTemp))),meanTemp+2*(sd(predTemp)/sqrt(length(predTemp))))
  out<-c(meanTemp,CI)
  if (min) {out<-c(out,sum(predTemp <= 4))} else {
    out<-c(out,sum(predTemp >= 15))
  }
  return(list(stats=out,predTemp=predTemp))
}

#Naylor min, mean=6.434571, CI=5.733979  7.135163 8 days with min<=4
naylorPred<-fitPred('naylorMin',c("Echo_min","Jack_min"),data2023[data2023$doy<=270,],min=T)

#Camp min, mean=6.035329, CI=5.362951  6.707707  11 days with min<=4
campPred<-fitPred('campMin',c("Echo_min","Jack_min"),data2023,min=T)

meanPred<-rowMeans(data.frame(naylorPred$predTemp),campPred$predTemp)
mean(meanPred) #6.4
sum(meanPred <=4) #8
sum(meanPred <=2) #4

xbreaks<-0:20
#pdf('meanTree2019_5AugTo27SeptemberMinTemps.pdf')
hist(meanPred,las=1,main="Predicted Min., 8/5/23 - 9/27/23, mean Min. = 6.4C",breaks=xbreaks)
#dev.off()

#Naylor max, mean=13.83840, CI=13.22007 14.45673,  days with min<=5
naylorPred<-fitPred('naylorMax',c("Echo_max","Jack_max"),data2023[data2023$doy<=270,],min=F)

#Camp max, mean=14.26872, CI=13.24017 15.29727, 63 days with max>=10
campPred<-fitPred('campMax',c("Echo_max","Jack_max"),data2023,min=F)

meanPred<-rowMeans(data.frame(naylorPred$predTemp),campPred$predTemp)
mean(meanPred) #13.8
sum(meanPred >=15) #16

xbreaks<-0:20
#pdf('meanTree2019_5AugTo27SeptemberMaxTemps.pdf')
hist(meanPred,las=1,main="Predicted Max., 8/5/23 - 9/27/23, mean Max. = 13.8C",breaks=xbreaks)
#dev.off()
