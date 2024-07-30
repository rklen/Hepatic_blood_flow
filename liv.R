linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
}
meanRelError<-function(modelCurve,organCurve){
  r<-0
  for(i in 1:length(organCurve)){
    if(organCurve[i]>0){
      if(modelCurve[i]!=organCurve[i]){
        r<-r+abs(modelCurve[i]-organCurve[i])/organCurve[i]
      }
    }
  }
  r<-r/(length(organCurve)-1)
  return(r)
}
wilcoxSymbol<-function(x,y){
  s<-''
  p<-wilcox.test(x,y,alternative='less',paired=TRUE)$p.value
  if(p<0.05){s<-'*'}
  if(p<0.01){s<-'**'}
  if(p<0.001){s<-'***'}
  return(s)
}

filepath<-'C:/Users/Oona/Documents/Tpc/lok/Info for Oona.txt'
df<-read.table(filepath)
df<-df[c(2:9,11:28,30:34,36:61),]
studyNumbers<-df$V1
t<-c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,90,100,120,140,
     160,190,220,250,280)
times<-c(0:280)

#1st model
fit_model<-function(param,input,t1){
  f<-abs(param[1])
  k<-abs(param[2])
  Vb<-1/(1+exp(-param[3]))
  cT<-0
  for(i in 2:length(input)){
    if(i-t1-1<1){c0<-0}else{c0<-input[i-t1-1]}
    cT<-c(cT,f*c0+(1-k)*cT[i-1])
  }
  v<-(1-Vb)*cT+Vb*input
  return(v)
}
errorfunction<-function(param,input,organCurve,t1){
  return(sum((organCurve-fit_model(param,input,t1))^2))
}
optimizationAlg<-function(initialValues,input,organCurve,t1){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,input=input,
               organCurve=organCurve,t1=t1,stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=16)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/liv/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,6)
  input<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[1,])
    organCurve[l]<-linearInter(times[l],t,df[3,])
  }
  #plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(1,1,-log(9))
  errorsFort<-c()
  for(t1 in ((0:6)*10)){
    param<-optimizationAlg(initialValues,input,organCurve,t1)
    err<-errorfunction(param,input,organCurve,t1)
    errorsFort<-c(errorsFort,err)
  }
  t1_l<-c((0:6)*10)[which.min(errorsFort)]
  errorsFort<-c()
  for(t1 in max(0,t1_l-5):(t1_l+5)){
    param<-optimizationAlg(initialValues,input,organCurve,t1)
    err<-errorfunction(param,input,organCurve,t1)
    errorsFort<-c(errorsFort,err)
  }
  t1<-c(max(0,t1_l-5):(t1_l+5))[which.min(errorsFort)]
  param<-optimizationAlg(initialValues,input,organCurve,t1)
  points(times,fit_model(param,input,t1),type='l')
  df1[i,1:length(param)]<-param
  df1[i,6]<-t1
  df1[i,9]<-meanRelError(fit_model(param,input,t1),organCurve)
  err<-errorfunction(param,input/1000,organCurve/1000,t1)
  df1[i,10]<-err/(length(times)-1)
  df1[i,11]<-280*log(df1[i,10])+2*4
  df1[i,12]<-1-err/sum((organCurve/1000-mean(organCurve/1000))^2)
  mod1<-cbind(times[-1],fit_model(param,input,t1)[-1])[times=t,2]/1000
  oc1<-cbind(times[-1],organCurve[-1])[times=t,2]/1000
  df1[i,13]<-meanRelError(mod1,oc1)
  err<-sum((mod1-oc1)^2)
  df1[i,14]<-err/length(oc1)
  df1[i,15]<-length(oc1)*log(df1[i,12])+2*4
  df1[i,16]<-1-err/sum((oc1-mean(oc1))^2)
}
write.csv(df1,file=paste('liv_model1.csv',sep=''),row.names=F)

#2nd model
fit_model<-function(param,input,t1,t2,t3){
  fA<-abs(param[1])
  fp<-abs(param[2])
  k<-abs(param[3])
  E<-1/(1+exp(-param[4]))
  Vb<-1/(1+exp(-param[5]))
  cS<-0
  cT<-0
  for(i in 2:length(input)){
    if(i-t2-1<1){c0<-0}else{c0<-input[i-t2-1]}
    cS<-c(cS,fp*c0+(1-fp)*cS[i-1])
    if(i-t1-1<1){c0<-0}else{c0<-input[i-t1-1]}
    if(i-t3-1<1){cs1<-0}else{cs1<-cS[i-t3-1]}
    cT<-c(cT,fA*c0+E*fp*cs1+(1-k)*cT[i-1])
  }
  v<-(1-Vb)*cT+Vb*(fA*input+E*fp*cS)/(fA+E*fp)
  return(rbind(v,cS))
}
errorfunction<-function(param,input,organCurve,spleen,t1,t2,t3){
  l<-fit_model(param,input,t1,t2,t3)
  cPET<-l[1,]
  cS<-l[2,]
  return(sum((organCurve-cPET)^2)+sum((spleen-cS)^2))
}
optimizationAlg<-function(initialValues,input,organCurve,spleen,t1,t2,t3){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,input=input,
               organCurve=organCurve,spleen=spleen,t1=t1,t2=t2,t3=t3,
               stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=16)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/liv/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,6)
  input<-c()
  organCurve<-c()
  spleen<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[1,])
    organCurve[l]<-linearInter(times[l],t,df[3,])
    spleen[l]<-linearInter(times[l],t,df[4,])
  }
  #plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(0.5,0.6,1,log(9),-log(9))
  errorsFort<-c()
  for(l in 0:(7*7*7-1)){
    t1<-(l-l%%(7*7))/7/7*10
    t2<-(l%%(7*7)-l%%7)/7*10
    t3<-l%%7*10
    param<-optimizationAlg(initialValues,input,organCurve,spleen,t1,t2,t3)
    err<-errorfunction(param,input,organCurve,spleen,t1,t2,t3)
    errorsFort<-c(errorsFort,err)
  }
  l<-c(0:(7*7*7-1))[which.min(errorsFort)]
  t1_l<-(l-l%%(7*7))/7/7*10
  t2_l<-(l%%(7*7)-l%%7)/7*10
  t3_l<-l%%7*10
  errorsFort<-c()
  for(t1 in max(0,t1_l-5):(t1_l+5)){
    param<-optimizationAlg(initialValues,input,organCurve,spleen,t1,t2_l,t3_l)
    err<-errorfunction(param,input,organCurve,spleen,t1,t2_l,t3_l)
    errorsFort<-c(errorsFort,err)
  }
  t1<-c(max(0,t1_l-5):(t1_l+5))[which.min(errorsFort)]
  errorsFort<-c()
  for(t2 in max(0,t2_l-5):(t2_l+5)){
    param<-optimizationAlg(initialValues,input,organCurve,spleen,t1,t2,t3_l)
    err<-errorfunction(param,input,organCurve,spleen,t1,t2,t3_l)
    errorsFort<-c(errorsFort,err)
  }
  t2<-c(max(0,t2_l-5):(t2_l+5))[which.min(errorsFort)]
  errorsFort<-c()
  for(t3 in max(0,t3_l-5):(t3_l+5)){
    param<-optimizationAlg(initialValues,input,organCurve,spleen,t1,t2,t3)
    err<-errorfunction(param,input,organCurve,spleen,t1,t2,t3)
    errorsFort<-c(errorsFort,err)
  }
  t3<-c(max(0,t3_l-5):(t3_l+5))[which.min(errorsFort)]
  param<-optimizationAlg(initialValues,input,organCurve,spleen,t1,t2,t3)
  points(times,fit_model(param,input,t1,t2,t3)[1,],type='l')
  df1[i,1:length(param)]<-param
  df1[i,6]<-t1
  df1[i,7]<-t2
  df1[i,8]<-t3
  df1[i,9]<-meanRelError(fit_model(param,input,t1,t2,t3)[1,],organCurve)
  err<-sum((organCurve/1000-fit_model(param,input,t1,t2,t3)[1,]/1000)^2)
  df1[i,10]<-err/(length(times)-1)
  df1[i,11]<-280*log(df1[i,10])+2*8
  df1[i,12]<-1-err/sum((organCurve/1000-mean(organCurve/1000))^2)
  mod1<-cbind(times[-1],fit_model(param,input,t1,t2,t3)[1,][-1])[times=t,2]/1000
  oc1<-cbind(times[-1],organCurve[-1])[times=t,2]/1000
  df1[i,13]<-meanRelError(mod1,oc1)
  err<-sum((mod1-oc1)^2)
  df1[i,14]<-err/length(oc1)
  df1[i,15]<-length(oc1)*log(df1[i,12])+2*8
  df1[i,16]<-1-err/sum((oc1-mean(oc1))^2)
}
write.csv(df1,file=paste('liv_model2.csv',sep=''),row.names=F)

#3rd model
fit_model<-function(param,input,t1,t2,t3){
  fA<-abs(param[1])
  fp<-abs(param[2])
  k<-abs(param[3])
  E<-1/(1+exp(-param[4]))
  Vb<-1/(1+exp(-param[5]))
  cPV<-0
  cT<-0
  for(i in 2:length(input)){
    if(i-t2-1<1){c0<-0}else{c0<-input[i-t2-1]}
    cPV<-c(cPV,fp*c0+(1-fp)*cPV[i-1])
    if(i-t1-1<1){c0<-0}else{c0<-input[i-t1-1]}
    if(i-t3-1<1){cpv1<-0}else{cpv1<-cPV[i-t3-1]}
    cT<-c(cT,fA*c0+E*fp*cpv1+(1-k)*cT[i-1])
  }
  v<-(1-Vb)*cT+Vb*(fA*input+E*fp*cPV)/(fA+E*fp)
  return(v)
}
errorfunction<-function(param,input,organCurve,t1,t2,t3){
  return(sum((organCurve-fit_model(param,input,t1,t2,t3))^2))
}
optimizationAlg<-function(initialValues,input,organCurve,t1,t2,t3){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,input=input,
               organCurve=organCurve,t1=t1,t2=t2,t3=t3,
               stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=16)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/liv/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,6)
  input<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[1,])
    organCurve[l]<-linearInter(times[l],t,df[3,])
  }
  #plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(0.5,0.6,1,log(9),-log(9))
  errorsFort<-c()
  for(l in 0:(7*7*7-1)){
    t1<-(l-l%%(7*7))/7/7*10
    t2<-(l%%(7*7)-l%%7)/7*10
    t3<-l%%7*10
    param<-optimizationAlg(initialValues,input,organCurve,t1,t2,t3)
    err<-errorfunction(param,input,organCurve,t1,t2,t3)
    errorsFort<-c(errorsFort,err)
  }
  l<-c(0:(7*7*7-1))[which.min(errorsFort)]
  t1_l<-(l-l%%(7*7))/7/7*10
  t2_l<-(l%%(7*7)-l%%7)/7*10
  t3_l<-l%%7*10
  errorsFort<-c()
  for(t1 in max(0,t1_l-5):(t1_l+5)){
    param<-optimizationAlg(initialValues,input,organCurve,t1,t2_l,t3_l)
    err<-errorfunction(param,input,organCurve,t1,t2_l,t3_l)
    errorsFort<-c(errorsFort,err)
  }
  t1<-c(max(0,t1_l-5):(t1_l+5))[which.min(errorsFort)]
  errorsFort<-c()
  for(t2 in max(0,t2_l-5):(t2_l+5)){
    param<-optimizationAlg(initialValues,input,organCurve,t1,t2,t3_l)
    err<-errorfunction(param,input,organCurve,t1,t2,t3_l)
    errorsFort<-c(errorsFort,err)
  }
  t2<-c(max(0,t2_l-5):(t2_l+5))[which.min(errorsFort)]
  errorsFort<-c()
  for(t3 in max(0,t3_l-5):(t3_l+5)){
    param<-optimizationAlg(initialValues,input,organCurve,t1,t2,t3)
    err<-errorfunction(param,input,organCurve,t1,t2,t3)
    errorsFort<-c(errorsFort,err)
  }
  t3<-c(max(0,t3_l-5):(t3_l+5))[which.min(errorsFort)]
  param<-optimizationAlg(initialValues,input,organCurve,t1,t2,t3)
  points(times,fit_model(param,input,t1,t2,t3),type='l')
  df1[i,1:length(param)]<-param
  df1[i,6]<-t1
  df1[i,7]<-t2
  df1[i,8]<-t3
  df1[i,9]<-meanRelError(fit_model(param,input,t1,t2,t3),organCurve)
  err<-errorfunction(param,input/1000,organCurve/1000,t1,t2,t3)
  df1[i,10]<-err/(length(times)-1)
  df1[i,11]<-280*log(df1[i,10])+2*8
  df1[i,12]<-1-err/sum((organCurve/1000-mean(organCurve/1000))^2)
  mod1<-cbind(times[-1],fit_model(param,input,t1,t2,t3)[-1])[times=t,2]/1000
  oc1<-cbind(times[-1],organCurve[-1])[times=t,2]/1000
  df1[i,13]<-meanRelError(mod1,oc1)
  err<-sum((mod1-oc1)^2)
  df1[i,14]<-err/length(oc1)
  df1[i,15]<-length(oc1)*log(df1[i,12])+2*8
  df1[i,16]<-1-err/sum((oc1-mean(oc1))^2)
}
write.csv(df1,file=paste('liv_model3.csv',sep=''),row.names=F)

#4th model
fit_model<-function(param,input,pvInput,t1,t2){
  fA<-abs(param[1])
  fP<-abs(param[2])
  k<-abs(param[3])
  Vb<-1/(1+exp(-param[4]))
  cT<-0
  for(i in 2:length(input)){
    if(i-t1-1<1){c0<-0}else{c0<-input[i-t1-1]}
    if(i-t2-1<1){cpv1<-0}else{cpv1<-pvInput[i-t2-1]}
    cT<-c(cT,fA*c0+fP*cpv1+(1-k)*cT[i-1])
  }
  v<-(1-Vb)*cT+Vb*(fA*input+fP*pvInput)/(fA+fP)
  return(v)
}
errorfunction<-function(param,input,pvInput,organCurve,t1,t2){
  return(sum((organCurve-fit_model(param,input,pvInput,t1,t2))^2))
}
optimizationAlg<-function(initialValues,input,pvInput,organCurve,t1,t2){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,input=input,
               pvInput=pvInput,organCurve=organCurve,t1=t1,t2=t2,
               stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=16)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/liv/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,6)
  input<-c()
  pvInput<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[1,])
    pvInput[l]<-linearInter(times[l],t,df[6,])
    organCurve[l]<-linearInter(times[l],t,df[3,])
  }
  #plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(0.5,0.5,1,-log(9))
  errorsFort<-c()
  for(l in 0:(7*7-1)){
    t1<-(l%%(7*7)-l%%7)/7*10
    t2<-l%%7*10
    param<-optimizationAlg(initialValues,input,pvInput,organCurve,t1,t2)
    err<-errorfunction(param,input,pvInput,organCurve,t1,t2)
    errorsFort<-c(errorsFort,err)
  }
  l<-c(0:(7*7-1))[which.min(errorsFort)]
  t1_l<-(l%%(7*7)-l%%7)/7*10
  t2_l<-l%%7*10
  errorsFort<-c()
  for(t1 in max(0,t1_l-5):(t1_l+5)){
    param<-optimizationAlg(initialValues,input,pvInput,organCurve,t1,t2_l)
    err<-errorfunction(param,input,pvInput,organCurve,t1,t2_l)
    errorsFort<-c(errorsFort,err)
  }
  t1<-c(max(0,t1_l-5):(t1_l+5))[which.min(errorsFort)]
  errorsFort<-c()
  for(t2 in max(0,t2_l-5):(t2_l+5)){
    param<-optimizationAlg(initialValues,input,pvInput,organCurve,t1,t2)
    err<-errorfunction(param,input,pvInput,organCurve,t1,t2)
    errorsFort<-c(errorsFort,err)
  }
  t2<-c(max(0,t2_l-5):(t2_l+5))[which.min(errorsFort)]
  param<-optimizationAlg(initialValues,input,pvInput,organCurve,t1,t2)
  points(times,fit_model(param,input,pvInput,t1,t2),type='l')
  err<-errorfunction(param,input,pvInput,organCurve,t1,t2)
  df1[i,1:length(param)]<-param
  df1[i,6]<-t1
  df1[i,7]<-t2
  df1[i,9]<-meanRelError(fit_model(param,input,pvInput,t1,t2),organCurve)
  err<-errorfunction(param,input/1000,pvInput/1000,organCurve/1000,t1,t2)
  df1[i,10]<-err/(length(times)-1)
  df1[i,11]<-280*log(df1[i,10])+2*6
  df1[i,12]<-1-err/sum((organCurve/1000-mean(organCurve/1000))^2)
  mod1<-cbind(times[-1],fit_model(param,input,pvInput,t1,t2)[-1])[times=t,2]/1000
  oc1<-cbind(times[-1],organCurve[-1])[times=t,2]/1000
  df1[i,13]<-meanRelError(mod1,oc1)
  err<-sum((mod1-oc1)^2)
  df1[i,14]<-err/length(oc1)
  df1[i,15]<-length(oc1)*log(df1[i,12])+2*6
  df1[i,16]<-1-err/sum((oc1-mean(oc1))^2)
}
write.csv(df1,file=paste('liv_model4.csv',sep=''),row.names=F)

df1<-read.csv('C:/Users/Oona/Documents/Tpc/liv/liv_model4.csv')
y<-df1[,1]
print(c(round(mean(abs(y*60)),digits=3),round(sd(abs(y*60)),digits=3)))
print(c(round(mean(1/(1+exp(-y))),digits=3)*100,round(sd(1/(1+exp(-y))),digits=3)*100))
print(c(round(mean(y),digits=3),round(sd(y),digits=3)))

df1<-read.csv('C:/Users/Oona/Documents/Tpc/liv/liv_model4.csv')
#9 MRE 280, 10 MSE, 11 AIC, 12 R^2, 13 MRE 24, 14 MSE, 15 AIC, 16 R^2
print(median(df1[,13]))
print(median(df1[,14]))
print(median(df1[,15]))

i<-13
df1<-read.csv('C:/Users/Oona/Documents/Tpc/liv/liv_model4.csv')
y1<-df1[,i]
df1<-read.csv('C:/Users/Oona/Documents/Tpc/liv/liv_model1.csv')
y<-df1[,i]
wilcoxSymbol(y1,y)

filepath<-'C:/Users/Oona/Documents/Tpc/lok/Info for Oona.txt'
df<-read.table(filepath)
df<-df[2:61,]
df$V4<-c(90,94,79,85,73,82,81,108,63,92,64,117,99,59,89,100,103,
         101,78,86,87,77,114,81,90,78,65,94,121,61,95,64,60,92,
         101,69,62,64,132,98,130,67,123,80,NA,102,93,76,91,72,87,
         88,78,113,110,81,95,69,107,NA)
df$V5<-c(355,356,351,308,316,336,364,354,350,406,357,319,334,350,
         356,341,339,332,324,322,408,402,364,385,363,334,380,364,
         336,351,362,366,396,390,359,364,334,345,379,390,348,366,
         321,332,295,370,356,371,323,307,374,383,333,381,349,338,
         353,360,362,366)
df$V6<-c(168,180,172,169,180,172,170,173,166,174,158,185,164,166,
         159,177,172,170,179,176,174,178,184,160,164,167,163,164,
         176,152,183,161,158,192,180,164,158,161,167,171,182,167,
         190,171,NA,180,182,162,167,165,172,163,181,170,192,165,182,
         168,170,NA)
df$V7<-df[,4]/(df[,6]/100)^2
df<-df[c(1:8,10:27,29:33,35:60),]
df1<-read.csv('C:/Users/Oona/Documents/Tpc/liv/liv_model4.csv')
df$V8<-abs(df1[,1])*60
df$V9<-abs(df1[,2])*60
df$V10<-df$V8+df$V9

womendf<-df[df[,3]=="Female",]
mendf<-df[df[,3]=="Male",]
i<-10
print(median(df[,i]))
print(median(womendf[,i]))
print(median(mendf[,i]))
wilcox.test(womendf[,i],mendf[,i],paired=FALSE)
wilcox.test(womendf[,8]/womendf[,10],mendf[,8]/mendf[,10],paired=FALSE)

cor(df[,8],df[,9])
cor(womendf[,8],womendf[,9])
cor(mendf[,8],mendf[,9])

cor(df[,10],as.numeric(df[,2]))
cor(womendf[,10],as.numeric(womendf[,2]))
cor(mendf[,10],as.numeric(mendf[,2]))

cor(df[!is.na(df[,4]),10],df[!is.na(df[,4]),4])
cor(womendf[!is.na(womendf[,4]),10],womendf[!is.na(womendf[,4]),4])
cor(mendf[!is.na(mendf[,4]),10],mendf[!is.na(mendf[,4]),4])

cor(df[!is.na(df[,7]),10],df[!is.na(df[,7]),7])
cor(womendf[!is.na(womendf[,7]),10],womendf[!is.na(womendf[,7]),7])
cor(mendf[!is.na(mendf[,7]),10],mendf[!is.na(mendf[,7]),7])

plot(mendf[,8],mendf[,9],xlim=c(0,1),ylim=c(0,3),
     lwd=2,ylab='Portal HBF (mL/min/mL)',xlab='Arterial HBF (mL/min/mL)',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(womendf[,8],womendf[,9],xlim=c(0,1),ylim=c(0,3),
       pch=3,axes=FALSE,lwd=2)
fit<-lm(mendf[,9]~mendf[,8])
abline(fit,col='blue',lw=2)
fit<-lm(womendf[,9]~womendf[,8])
abline(fit,lw=2,lty=2)
legend('bottomright',legend=c('Women','Men'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)

plot(as.numeric(mendf[,2]),mendf[,10],xlim=c(40,85),ylim=c(0,3),
     lwd=2,ylab='Total HBF (mL/min/mL)',xlab='Age',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[,2]),womendf[,10],xlim=c(40,85),ylim=c(0,3),
       pch=3,axes=FALSE,lwd=2)
fit<-lm(mendf[,10]~as.numeric(mendf[,2]))
abline(fit,col='blue',lw=2)
fit<-lm(womendf[,10]~as.numeric(womendf[,2]))
abline(fit,lw=2,lty=2)
legend('topleft',legend=c('Women','Men'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)

plot(mendf[!is.na(mendf[,4]),4],mendf[!is.na(mendf[,4]),10],
     xlim=c(55,135),ylim=c(0,3),
     lwd=2,ylab='Total HBF (mL/min/mL)',xlab='Weight',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(womendf[!is.na(womendf[,4]),4],womendf[!is.na(womendf[,4]),10],
       xlim=c(40,85),ylim=c(0,3),
       pch=3,axes=FALSE,lwd=2)
fit<-lm(mendf[!is.na(mendf[,4]),10]~mendf[!is.na(mendf[,4]),4])
abline(fit,col='blue',lw=2)
fit<-lm(womendf[!is.na(womendf[,4]),10]~womendf[!is.na(womendf[,4]),4])
abline(fit,lw=2,lty=2)
legend('topright',legend=c('Women','Men'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)

filepath<-paste('C:/Users/Oona/Documents/Tpc/liv/array1_',studyNumbers[1],'.csv',sep='')
df<-read.csv(filepath)
df<-as.matrix(df)
colnames(df)<-NULL
rownames(df)<-NULL
df[,1]<-rep(0,6)
input<-c()
pvInput0<-c()
pvInput1<-c()
spleen<-c()
organCurve<-c()
for(l in 1:length(times)){
  input[l]<-linearInter(times[l],t,df[1,])/1000
  pvInput0[l]<-linearInter(times[l],t,df[2,])/1000
  pvInput1[l]<-linearInter(times[l],t,df[6,])/1000
  spleen[l]<-linearInter(times[l],t,df[4,])/1000
  organCurve[l]<-linearInter(times[l],t,df[3,])/1000
}
plot(times,input,type='l',
     xlim=c(0,280),ylim=c(0,100),
     lwd=2,ylab='Activity concentration (kBq/mL)',xlab='Time (s)',
     cex.axis=1.3,cex.lab=1.3)
points(times,pvInput0,type='l',
       xlim=c(0,280),ylim=c(0,100),lty=2,lwd=2,col='blue')
points(times,pvInput1,type='l',
       xlim=c(0,280),ylim=c(0,100),lwd=2,col='blue')
points(times,spleen,type='l',
       xlim=c(0,280),ylim=c(0,100),lwd=2,col='steelblue1')
points(times,organCurve,type='l',
       xlim=c(0,280),ylim=c(0,100),lwd=2,col='gray')
legend('topright',legend=c('Aorta',
                           'PV (orig. mean TAC)',
                           'PV (for the models)',
                           'Liver',
                           'Spleen'),
       lty=c(1,2,1,1,1),
       col=c('black','blue','blue','gray','steelblue1'),
       lwd=c(2,2,2,2,2),cex=1.4)
