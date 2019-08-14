#-----------------------------------------------------------#
# Empirical Version of Louis formula
#-----------------------------------------------------------#

#remove the old stuff
rm(list = ls())
#load some packages
library(censReg)
library(tmvtnorm)
library(matrixcalc)
library(psych)
library(mvtnorm)
library(MASS)
library(TruncatedNormal)

setwd("~/Workflow")

# Select the year of interest
t_select<-c(1993:2014)
for (t in t_select){

  load(paste("2) Estimating the Models/results/coef_",t,".RData",sep=""))
  
  
  #-----------------------------------------------------------#
  # Pre-Loop Calculations 
  #-----------------------------------------------------------#
  
  # variance
  sigma_est<-theta_hat[length(theta_hat)]
  beta<-theta_hat[4:(length(theta_hat)-1)]
  
  # A_est 
  A_est<-diag(length(Y_cens))-theta_hat[1]*W1_fix-theta_hat[2]*W2_fix-theta_hat[3]*W3_fix
  # A_minus=A_est inverted
  A_minus<-solve(A_est)
  
  # Calculate the traces of the interaction terms of
  # A_minus and the weight matrices
  invAW1_der1<-tr(A_minus%*%W1_fix%*%A_minus%*%W1_fix)
  invAW1_der2<-tr(A_minus%*%W1_fix%*%A_minus%*%W2_fix)
  invAW1_der3<-tr(A_minus%*%W1_fix%*%A_minus%*%W3_fix)
  invAW2_der2<-tr(A_minus%*%W2_fix%*%A_minus%*%W2_fix)
  invAW2_der1<-invAW1_der2
  invAW2_der3<-tr(A_minus%*%W2_fix%*%A_minus%*%W3_fix)
  invAW3_der3<-tr(A_minus%*%W3_fix%*%A_minus%*%W3_fix)
  invAW3_der1<-invAW1_der3
  invAW3_der2<-invAW2_der3
  
  #-----------------------------------------------------------#
  # We need to rearrange everything according
  # to the missingnes patterns 
  #-----------------------------------------------------------#
  
  # A_rear
  A_rear<-Cov_reaarrange(A_est)$CovMat
  # X_rear
  X_rear<-rbind(X[noncens_index,],X[cens_index,])
  
  # W1_rear
  W1_rearr<-Cov_reaarrange(W1_fix)$CovMat
  # W2_rear
  W2_rearr<-Cov_reaarrange(W2_fix)$CovMat
  # W3_rear
  W3_rearr<-Cov_reaarrange(W3_fix)$CovMat
  
  ATA_rear<-Cov_reaarrange(t(A_est)%*%A_est)$CovMat
  
  Q1<- -W1_fix-t(W1_fix)+2*theta_hat[1]*t(W1_fix)%*%W1_fix+theta_hat[2]*(t(W1_fix)%*%W2_fix+t(W2_fix)%*%W1_fix)+theta_hat[3]*(t(W1_fix)%*%W3_fix+t(W3_fix)%*%W1_fix)
  Q2<- -W2_fix-t(W2_fix)+2*theta_hat[2]*t(W2_fix)%*%W2_fix+theta_hat[1]*(t(W1_fix)%*%W2_fix+t(W2_fix)%*%W1_fix)+theta_hat[3]*(t(W2_fix)%*%W3_fix+t(W3_fix)%*%W2_fix)
  Q3<- -W3_fix-t(W3_fix)+2*theta_hat[3]*t(W3_fix)%*%W3_fix+theta_hat[1]*(t(W1_fix)%*%W3_fix+t(W3_fix)%*%W1_fix)+theta_hat[2]*(t(W2_fix)%*%W3_fix+t(W3_fix)%*%W2_fix)
  
  Q1_rear<-Cov_reaarrange(Q1)$CovMat
  Q2_rear<-Cov_reaarrange(Q2)$CovMat
  Q3_rear<-Cov_reaarrange(Q3)$CovMat
  
  
  XB<-X_rear%*%beta
  
  
  tr_1<- -tr(A_minus%*%W1_fix)
  tr_2<- -tr(A_minus%*%W2_fix)
  tr_3<- -tr(A_minus%*%W3_fix)
  
  
  WW11<-Cov_reaarrange(t(W1_fix)%*%W1_fix)$CovMat
  WW22<-Cov_reaarrange(t(W2_fix)%*%W2_fix)$CovMat
  WW33<-Cov_reaarrange(t(W3_fix)%*%W3_fix)$CovMat
  WW12<-Cov_reaarrange(t(W1_fix)%*%W2_fix+t(W2_fix)%*%W1_fix)$CovMat
  WW13<-Cov_reaarrange(t(W1_fix)%*%W3_fix+t(W3_fix)%*%W1_fix)$CovMat
  WW32<-Cov_reaarrange(t(W3_fix)%*%W2_fix+t(W2_fix)%*%W3_fix)$CovMat
  
  
  dell_beta_beta<- -(1/(sigma_est))*t(X_rear)%*%X_rear
  
  
  
  
  

  E_step3<-function(A,Cov,Y,beta,samplesize=1000){
    out<-Cov_reaarrange(Cov)
    
    mu<-solve(A,X%*%beta)

    mu_o<-mu[out$noncens]
    mu_m<-mu[out$cens]
    
    
    mu_cond<-mu_m+out$Sigma_mo%*%solve(out$Sigma_oo)%*%(Y[out$noncens]-mu_o) 
    Sigma_cond<-out$Sigma_mm-out$Sigma_mo%*%solve(out$Sigma_oo)%*%t(out$Sigma_mo)
    
    mu_cond_save<-mu_cond
    sample<-mvrandn(l=rep(-Inf, length = length(mu_cond)),u=rep(minimum, length = length(mu_cond))-mu_cond,Sig=Sigma_cond,n=samplesize)
    sample<-sample+rep(mu_cond_save,samplesize)
    mu_cond<-t(sample)
    
    
    
    
    
    return(mu_cond)
  }
  
  
  
  
  result_full<-E_step3(A_est,A_minus%*%t(A_minus)*sigma_est,Y_cens,beta)
  
  
  ##
  
  I_louis<-list()
  re_save<-c()
  rep<-dim(result_full)[1]
  
  for (p in 1:rep){
    
    
    print(p)
    
    # insert the imputed values
    Y_imputed<-Y_cens
    Y_imputed[cens_index]<-result_full[p,]
    Y_imputed_rear<-c(Y_imputed[noncens_index],Y_imputed[cens_index])
    
    # calcluate the needed quantities
    S_star_A=t(Y_imputed_rear)%*%ATA_rear%*%Y_imputed_rear
    S_star_Q1=t(Y_imputed_rear)%*%Q1_rear%*%Y_imputed_rear
    S_star_Q2=t(Y_imputed_rear)%*%Q2_rear%*%Y_imputed_rear
    S_star_Q3=t(Y_imputed_rear)%*%Q3_rear%*%Y_imputed_rear
    
    
    dell_rho1<-tr_1-1/(2*sigma_est)*(S_star_Q1+2*t(XB)%*%W1_rearr%*%Y_imputed_rear)
    dell_rho2<-tr_2-1/(2*sigma_est)*(S_star_Q2+2*t(XB)%*%W2_rearr%*%Y_imputed_rear)
    dell_rho3<-tr_3-1/(2*sigma_est)*(S_star_Q3+2*t(XB)%*%W3_rearr%*%Y_imputed_rear)
    dell_beta<-1/(sigma_est)*(t(X_rear)%*%(A_rear%*%Y_imputed_rear-XB))
    dell_sigma<- -(n_edges/(2*sigma_est))+1/(2*sigma_est^2)*(S_star_A-2*t(XB)%*%A_rear%*%Y_imputed_rear+t(XB)%*%XB)
    re<-c(dell_rho1,dell_rho2,dell_rho3,dell_beta,dell_sigma)
    
    
    S_star_W1W1=t(Y_imputed_rear)%*%WW11%*%Y_imputed_rear
    S_star_W2W2=t(Y_imputed_rear)%*%WW22%*%Y_imputed_rear
    S_star_W3W3=t(Y_imputed_rear)%*%WW33%*%Y_imputed_rear
    S_star_W1W2=t(Y_imputed_rear)%*%WW12%*%Y_imputed_rear
    S_star_W2W1<-S_star_W1W2
    S_star_W1W3=t(Y_imputed_rear)%*%WW13%*%Y_imputed_rear
    S_star_W3W1<-S_star_W1W3
    S_star_W3W2=t(Y_imputed_rear)%*%WW32%*%Y_imputed_rear
    S_star_W2W3<-S_star_W3W2
    
    
    # first row
    dell_rho1_rho1<--invAW1_der1-(1/sigma_est)*( S_star_W1W1 )
    dell_rho1_rho2<--invAW1_der2-(1/(2*sigma_est))*(S_star_W1W2  )
    dell_rho1_rho3<--invAW1_der3-(1/(2*sigma_est))*( S_star_W1W3 )
    dell_rho1_beta<- -(1/sigma_est)*t(X_rear)%*%t(W1_rearr)%*%Y_imputed_rear
    dell_rho1_sigma<- (1/(2*sigma_est^2))*(S_star_Q1 + 2*t(XB)%*%W1_rearr%*%Y_imputed_rear)
    r1<-c(dell_rho1_rho1,dell_rho1_rho2,dell_rho1_rho3,dell_rho1_beta,dell_rho1_sigma)
    
    # second row
    dell_rho2_rho1<-dell_rho1_rho2
    dell_rho2_rho2<--invAW2_der2-(1/sigma_est)*( S_star_W2W2  )
    dell_rho2_rho3<--invAW2_der3-(1/(2*sigma_est))*( S_star_W2W3 )
    dell_rho2_beta<- -(1/sigma_est)*t(X_rear)%*%t(W2_rearr)%*%Y_imputed_rear
    dell_rho2_sigma<- (1/(2*sigma_est^2))*(S_star_Q2 + 2*t(XB)%*%W2_rearr%*%Y_imputed_rear)
    r2<-c(dell_rho2_rho1,dell_rho2_rho2,dell_rho2_rho3,dell_rho2_beta,dell_rho2_sigma)
    
    
    # third row
    dell_rho3_rho1<-dell_rho1_rho3
    dell_rho3_rho2<-dell_rho2_rho3
    dell_rho3_rho3<--invAW3_der3-(1/sigma_est)*( S_star_W3W3  )
    dell_rho3_beta<- -(1/sigma_est)*t(X_rear)%*%t(W3_rearr)%*%Y_imputed_rear
    dell_rho3_sigma<- (1/(2*sigma_est^2))*(S_star_Q3 + 2*t(XB)%*%W3_rearr%*%Y_imputed_rear)
    r3<-c(dell_rho3_rho1,dell_rho3_rho2,dell_rho3_rho3,dell_rho3_beta,dell_rho3_sigma)
    
    
    # fourth row to seventh row
    dell_beta_rho1<-dell_rho1_beta
    dell_beta_rho2<-dell_rho2_beta
    dell_beta_rho3<-dell_rho3_beta
    dell_beta_sigma<- -(1/sigma_est^2)*(t(X_rear)%*%(A_rear%*%Y_imputed_rear-XB))
    r4_7<-cbind(dell_beta_rho1,dell_beta_rho2,dell_beta_rho3,dell_beta_beta,dell_beta_sigma)
    
    # eigth row
    dell_sigma_rho1<-dell_rho1_sigma
    dell_sigma_rho2<-dell_rho2_sigma
    dell_sigma_rho3<-dell_rho3_sigma
    dell_sigma_beta<-dell_beta_sigma
    dell_sigma_sigma<-(n_edges/(2*sigma_est^2))-(1/(sigma_est^3))*(S_star_A-2*t(XB)%*%A_rear%*%Y_imputed_rear+t(XB)%*%XB )
    
    r8<-c(dell_sigma_rho1,dell_sigma_rho2,dell_sigma_rho3,dell_sigma_beta,dell_sigma_sigma)
    
    I_comp<- -rbind(r1,r2,r3,r4_7,r8)
    
    I_louis[[p]]<- I_comp
    re_save<-cbind(re_save,re)
    
  }
  
  

  incl<-1:1000
  S<-rowMeans(re_save[,incl])
  V<-var(t(re_save[,incl]))
  # Here we apply Louis Formula
  E_H<- (1/length(I_louis[incl]))*( Reduce("+",I_louis[incl])) 
  E_FI<-E_H-V
  AVCOV<-solve(E_FI)
  var<-diag(AVCOV)
  sd<-sqrt(var)
  print(sd)
  
  # Remove everything that is not needed anymore
  rm(list=ls()[-c(which(ls()=="theta_hat"),which(ls()=="sd"),which(ls()=="t"))])
  
  # Save the results
  save.image(paste("2) Estimating the Models/standard errors/sd_",t,".RData",sep=""))
  
}


