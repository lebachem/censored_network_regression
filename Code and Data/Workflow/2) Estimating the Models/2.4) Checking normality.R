#-----------------------------------------------------------#
# Extract the conditional proabilities
#-----------------------------------------------------------#

rm(list=ls())
library(Matrix)
library(expm)


setwd("H:/Submitted/JRSSC/Code and Data/Workflow")


t_select<-c(2014:1993)
Y1<-list()
Y1_std<-list()
res<-list()
res_std<-list()
for (t in t_select){
  
  
  
  load(paste("2) Estimating the Models/results/coef_",t,".RData",sep=""))
  
  # calcluate B
  B_est<-solve(A_est)
  
  # Calculate Sigma
  Sigma<-B_est%*%t(B_est)*theta_hat[length(theta_hat)]
  
  # Calculate the observed Sigma
  Sigma_obs<-Sigma[noncens_index,noncens_index]
  
  # Calculate the mean 
  mu<-B_est%*%X%*%theta_hat[4:(length(theta_hat)-1)]
  
  # Calculate the observed mean
  mu_obs<-mu[noncens_index]
  E <- eigen(Sigma) 
  V <- E$values 
  Q <- E$vectors 
  S_sqrtm <- Q%*%diag(1/sqrt(V))%*%t(Q)   
  
  resi<-S_sqrtm[noncens_index,noncens_index]%*%(Y[noncens_index]-mu_obs)
  resi_std<-(resi-mean(resi))/(sd(resi))
  
  qqnorm(resi_std)
  qqline(resi_std)
  
  save(resi,file=paste("2) Estimating the Models/normality/residual_",t,".RData",sep=""))
  
  res[[t-1992]]<-resi
  
  resi_std<-(resi-mean(resi))/(sd(resi))
  res_std[[t-1992]]<-resi_std
  save(resi_std,file=paste("2) Estimating the Models/normality/residual_std_",t,".RData",sep=""))
  
  response<-Y[Y>0]
  save(response,file=paste("2) Estimating the Models/normality/response_",t,".RData",sep=""))
  
  Y1[[t-1992]]<-response
  response_std<-(response-mean(response))/sd(response) 
  
  Y1_std[[t-1992]]<-response_std
  save(response_std,file=paste("2) Estimating the Models/normality/response_std_",t,".RData",sep=""))
  
  
}







