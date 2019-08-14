
#-----------------------------------------------------------#
# Extract the conditional proabilities
#-----------------------------------------------------------#

rm(list=ls())
library(Matrix)
library(expm)


setwd("H:/Submitted/JRSSC/Code and Data/Workflow")


t_select<-c(1993:2014)

Y<-c()
for (t in t_select){
  
  
  
  #load(paste("2) Estimating the Models/results/coef_",t,".RData",sep=""))
  
  
  
  load(paste("2) Estimating the Models/normality/response_std_",t,".RData",sep=""))
  Y<-c(Y,response_std)


}

std<-function(x){(x-mean(x))/sd(x)}


pdf("2) Estimating the Models/normality/normality_response.pdf", width = 12,height=6)
par(mfrow=c(1,2))
plot(density(std(Y),adjust=2),main="Kernel Density Estimate - log(Y_obs)")
qqnorm(std(Y))
qqline(std(Y))
dev.off()




t_select<-c(1993:2014)

pdf("2) Estimating the Models/normality/normality_resid.pdf", width = 12,height=20)
par(mfrow=c(6,4))

for (t in t_select){
  
  
  load(paste("2) Estimating the Models/normality/residual_",t,".RData",sep=""))

  qqnorm(std(resi))
  qqline(std(resi))
}
dev.off()


