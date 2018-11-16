## Exercise 10
## Iker Soto
setwd("Biocomp-Fall2018-181109-Exercise10/")
data<-read.csv("data.txt")
head(data)
library(ggplot2)
ggplot(data,aes(x=x,y=y))+geom_point()+theme_classic()

#Custom likelihood function
#Linear model
linelike<-function(p,x,y){
  A=p[1]
  B=p[2]
  sigma=exp(p[3])
  expected=A+B*x
  nl=-sum(dnorm(x=y, mean=expected,sd=sigma,log=TRUE))
  return(nl)
}

#Hump-shaped model
humplike<-function(p,x,y){
  A=p[1]
  B=p[2]
  C=p[3]
  sigma=exp(p[4])
  expected=A+(B*x)+C*(x^2)
  nh=-sum(dnorm(x=y, mean=expected,sd=sigma,log=TRUE))
  return(nh)
}

#Minimizing negative log likelihoods of linear model
Guess1=c(1,1,1)
fit=optim(par=Guess1,fn=linelike, x=data$x, y=data$y)
print(fit)

#Minimizing negative log likelihoods of hump-like model
Guess2=c(1,1,1,1)
fit1=optim(par=Guess2,fn=humplike, x=data$x,y=data$y)
print(fit1)

#Statistical analysis
test=2*(fit$value-fit1$value)
1-pchisq(q=test,df=1)
#This analysis gives a result of one, which means the linear model is the best fit for our data.

####

#Exercise 2 
#Custom function
ddCompetition<-function(t,y,p){
  N1=y[1]
  N2=y[2]
  r1=p[1]
  r2=p[2]
  a11=p[3]
  a12=p[4]
  a22=p[5]
  a21=p[6]
  dN1dt=r1*(1-(N1*a11)-(N2*a12))*N1
  dN2dt=r2*(1-(N2*a22)-(N1*a21))*N2
  return(list(c(dN1dt,dN2dt)))
}

#Case 1 Coexistence : a12<a11 and a21<a22
library(deSolve)
library(reshape2)
y0=c(0.1,0.1)
params=c(0.5,0.5,0.009,0.005,0.0075,0.006)
times=1:100
modelSimcomp=ode(y=y0,times=times,func=ddCompetition,parms=params)
modelOutputSimcomp=data.frame(time=modelSimcomp[,1],specie1=modelSimcomp[,2],specie2=modelSimcomp[,3])
modelOutputSimcomp=melt(modelOutputSimcomp,id.vars = "time")
ggplot(modelOutputSimcomp,aes(x=time,y=value))+geom_line(aes(color=variable))+theme_classic()

#Case 2 Specie 1 is outcompeted by Specie 2 a12>a11. Population os Specie 1 is driven to extinction
y0.1=c(0.1,0.1)
params1=c(0.5,0.5,0.005,0.009,0.0075,0.006)
times1=1:100
modelSimcomp1=ode(y=y0.1,times=times1,func=ddCompetition,parms=params1)
modelOutputSimcomp1=data.frame(time=modelSimcomp1[,1],specie1=modelSimcomp1[,2],specie2=modelSimcomp1[,3])
modelOutputSimcomp1=melt(modelOutputSimcomp1,id.vars = "time")
ggplot(modelOutputSimcomp1,aes(x=time,y=value))+geom_line(aes(color=variable))+theme_classic()

#Case 3 Specie 2 is outcompeted by Specie 1 a21>a22. Population of Specie 2 is driven to extinction.
y0.2=c(0.1,0.1)
params2=c(0.5,0.5,0.009,0.005,0.0075,0.0099)
times2=1:100
modelSimcomp2=ode(y=y0.2,times=times2,func=ddCompetition,parms=params2)
modelOutputSimcomp2=data.frame(time=modelSimcomp2[,1],specie1=modelSimcomp2[,2],specie2=modelSimcomp2[,3])
modelOutputSimcomp2=melt(modelOutputSimcomp2,id.vars = "time")
ggplot(modelOutputSimcomp2,aes(x=time,y=value))+geom_line(aes(color=variable))+theme_classic()

#Case 4 The coexistence criteria isn't statisfied for either population: a12>a11 & a21>a22. Population 
# of Specie 1 is driven to extinction
y0.3=c(0.1,0.1)
params3=c(0.5,0.5,0.005,0.009,0.006,0.0075)
times3=1:100
modelSimcomp3=ode(y=y0.3,times=times3,func=ddCompetition,parms=params3)
modelOutputSimcomp3=data.frame(time=modelSimcomp3[,1],specie1=modelSimcomp3[,2],specie2=modelSimcomp3[,3])
modelOutputSimcomp3=melt(modelOutputSimcomp3,id.vars = "time")
ggplot(modelOutputSimcomp3,aes(x=time,y=value))+geom_line(aes(color=variable))+theme_classic()
