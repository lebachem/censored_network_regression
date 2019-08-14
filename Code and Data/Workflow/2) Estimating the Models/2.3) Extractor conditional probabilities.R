#-----------------------------------------------------------#
# Extract the conditional proabilities
#-----------------------------------------------------------#

rm(list=ls())
library(Matrix)



setwd("~/Workflow")

t_select<-c(1993:2014)
int_positive<-c()
for (t in t_select){

  load(paste("2) Estimating the Models/results/coef_",t,".RData",sep=""))
  

  A_inv<-chol2inv(chol(A_est))
  
  # Calculate mu hat
  mu<-A_inv%*%X%*%theta_hat[4:(length(theta_hat)-1)]
  # Calculate sigma hat
  Cov<-A_inv%*%t(A_inv)*theta_hat[length(theta_hat)]
  
  # Define Z_{-ij}
  Z_minus_ij<-Y
  # Fill the empty entries with the conditional expectation
  Z_minus_ij[cens_index]<-mu[cens_index]
  
  # Defne the function to create conditional expectations and variances
  Cov_reaarrange<-function(Cov,number){
    
    zero_ind<-number
    nonzero_ind<-1:n_edges
    nonzero_ind<-nonzero_ind[-number]
    n_zero=length(zero_ind)
    n_greater_zero=length(Y)-n_zero
    
    Cov<-Cov[c(nonzero_ind,zero_ind),c(nonzero_ind,zero_ind)] 
    
    out<-list(cens=zero_ind,noncens=nonzero_ind,CovMat=Cov,Sigma_oo=Cov[1:n_greater_zero,1:n_greater_zero],Sigma_om=Cov[1:n_greater_zero,(n_greater_zero+1):length(Y)],Sigma_mo=Cov[(n_greater_zero+1):length(Y),1:n_greater_zero],Sigma_mm=Cov[(n_greater_zero+1):length(Y),(n_greater_zero+1):length(Y)])
    return(out)
  }
  
  probs<-c()
  
  for (i in 1:length(Y)){
    out<-Cov_reaarrange(Cov,i)

    mu_o<-mu[out$noncens]
    mu_m<-mu[out$cens]
    Sigma_oo_inv<-    chol2inv(chol(out$Sigma_oo))
    
    mu_cond<-mu_m+out$Sigma_mo%*%Sigma_oo_inv%*%(Z_minus_ij[out$noncens]-mu_o) 
    Sigma_cond<-out$Sigma_mm-out$Sigma_mo%*%Sigma_oo_inv%*%out$Sigma_mo
    
    probs<-c(probs,1-pnorm(minimum,mean=mu_cond,sd=Sigma_cond))
    print(i)
  }
  
  # Remove everything that is not needed
  rm(list=ls()[-c(which(ls()=="probs"),which(ls()=="Y"),which(ls()=="t"),which(ls()=="t_select"))])
  
  # Save the results  
  save.image(paste("2) Estimating the Models/probs/probs_",t,".RData",sep=""))

} 
