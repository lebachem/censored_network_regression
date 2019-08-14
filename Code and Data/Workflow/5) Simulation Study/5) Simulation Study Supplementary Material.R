#-----------------------------------------------------------#
# Simulation Study showing unbiased results
#-----------------------------------------------------------#
#-----------------------------------------------------------#
# Simulation Study building on DGP1
#-----------------------------------------------------------#

rm(list=ls())

library(censReg)
library(tmvtnorm)
library(matrixcalc)
library(psych)
library(TruncatedNormal)
library(MASS)
library(mvtnorm)
library(pROC)

# set working directory
setwd("~/Workflow")

# Parameters for DGP1----
set.seed(123)
n_nodes=20
n_edges=n_nodes*(n_nodes-1)
rho_1=0.1
rho_2=0.2
rho_3=0.3
beta<-c(1,2,3,4,5)
mean_x<-1
sigma_x<-1
sigma<-1
censoring<-0.75
rep<-100
stop<-0.0001  
# Containers for the Results ----
share<-c()
share2<-c()
real<-c()
estimated<-c()
parameter_est<-c()
TP<-c()
FP<-c()
FN<-c()
TN<-c()

# Define some functions ----
# Constructur function for the endogenous effects
A_constructur<-function(Y,n_nodes,b_1,b_2,b_3){
  
  first<-rep(1:n_nodes,each=(n_nodes-1))
  second<-rep(1:n_nodes,times=(n_nodes))
  
  count<-length(second)
  for (i in 1:n_nodes){
    second<-second[-count]
    count<-count-(n_nodes+1)
  }
  
  W1<-diag(rep(0,length(Y)))
  for (i in 1:dim(W1)[1]){
    
    W1[i,which(first[i]==first)]<- 1
    
  }
  diag(W1)<-0
  W1<-W1%*%diag(1/rowSums(W1))*(b_2)
  
  
  
  W2<-diag(rep(0,length(Y)))
  for (i in 1:dim(W2)[1]){
    W2[i,which(second[i]==second)]<- 1
  }
  diag(W2)<-0
  W2<-W2%*%diag(1/rowSums(W2))*(b_3)
  
  W3<-diag(rep(0,length(Y)))
  for (i in 1:dim(W3)[1]){
    
    W3[i,which(second[i]==first&first[i]==second)]<- 1
  }
  diag(W3)<-0
  W3<-W3%*%diag(1/rowSums(W3))*(b_1)
  
  W<-W1+W2+W3
  
  A<-diag(rep(1,length(Y)))-W
  
  return(A)
}

# Define the Covariance-rearranging function
Cov_reaarrange<-function(Cov){
  
  zero_ind<-which(Y_cens==0)
  nonzero_ind<-which(Y_cens!=0)
  n_zero=length(zero_ind)
  n_greater_zero=length(Y)-n_zero
  
  Cov<-Cov[c(nonzero_ind,zero_ind),c(nonzero_ind,zero_ind)] 
  
  out<-list(cens=zero_ind,noncens=nonzero_ind,CovMat=Cov,Sigma_oo=Cov[1:n_greater_zero,1:n_greater_zero],Sigma_om=Cov[1:n_greater_zero,(n_greater_zero+1):length(Y)],Sigma_mo=Cov[(n_greater_zero+1):length(Y),1:n_greater_zero],Sigma_mm=Cov[(n_greater_zero+1):length(Y),(n_greater_zero+1):length(Y)])
  return(out)
}

E_step2<-function(A,Cov,Y,beta,trunc=F,exact=F,samplesize=1000,burn.in=100,thin=10){
  out<-Cov_reaarrange(Cov)
  
  # A must not be symmetric, hence we solve a linear equation system
  mu<-chol2inv(chol(A))%*%X%*%beta
  
  mu_o<-mu[out$noncens]
  mu_m<-mu[out$cens]
  
  Sigma_00_inv<-chol2inv(chol(out$Sigma_oo))
  mu_cond<-mu_m+out$Sigma_mo%*%Sigma_00_inv%*%(Y[out$noncens]-mu_o) 
  Sigma_cond<-out$Sigma_mm-out$Sigma_mo%*%Sigma_00_inv%*%t(out$Sigma_mo)
  
  
  if (trunc==T){
    
    mu_cond_save<-mu_cond
    sample<-mvrandn(l=rep(-Inf, length = length(mu_cond)),u=rep(minimum, length = length(mu_cond))-mu_cond,Sig=Sigma_cond,n=samplesize)
    sample<-sample+rep(mu_cond_save,samplesize)
    mu_cond<-rowMeans(sample)
    Sigma_cond<-var(t(sample))  
    
  }
  
  
  
  
  Y_imputed<-Y_cens
  Y_imputed[out$cens]<-mu_cond
  
  
  return(list(Y_imputed=Y_imputed,mu_cond=mu_cond,Sigma_cond=Sigma_cond,Y_obs=Y_imputed[out$noncens],mu_o=mu_o))
}

Q_grid<-  function(rho){
  
  #extract parameters
  b_1<-rho[1]
  b_2<-rho[2]
  b_3<-rho[3]
  
  #form A
  A<-diag(rep(1,length(Y)))-b_1*W1_fix-b_2*W2_fix-b_3*W3_fix
  
  #imputed Y star
  Y_star<-c(result$Y_imputed[noncens_index],result$Y_imputed[cens_index])
  
  #prepare S star
  AtA_oo<-Cov_reaarrange(t(A)%*%A)$Sigma_oo
  AtA_mo<-Cov_reaarrange(t(A)%*%A)$Sigma_mo
  AtA_om<-Cov_reaarrange(t(A)%*%A)$Sigma_om
  AtA_mm<-Cov_reaarrange(t(A)%*%A)$Sigma_mm
  
  E_star1=result$Y_obs%*%AtA_oo%*%(result$Y_obs)
  E_star2=result$Y_obs%*%AtA_om%*%result$mu_cond
  E_star3=t(result$mu_cond)%*%AtA_mo%*%result$Y_obs
  E_star4=tr(AtA_mm%*%result$Sigma_cond)+t(result$mu_cond)%*%AtA_mm%*%result$mu_cond
  S_star=E_star1+E_star2+E_star3+E_star4
  
  A_rear<-Cov_reaarrange(A)$CovMat
  
  
  concentrate<- sum(log(eigen(A)$values)) -(n_edges/2)*log( S_star - t(Y_star)%*%t(A_rear)%*%X_rear%*%solve(t(X_rear)%*%X_rear)%*%t(X_rear)%*%A_rear%*%Y_star   )
  
  
  return(-as.numeric(concentrate))
}

Q_grid_gradient<-function(rho){
  
  #extract parameters
  b_1<-rho[1]
  b_2<-rho[2]
  b_3<-rho[3]
  
  #form A
  A<-diag(rep(1,length(Y)))-b_1*W1_fix-b_2*W2_fix-b_3*W3_fix
  
  #imputed Y star
  Y_star<-c(result$Y_imputed[noncens_index],result$Y_imputed[cens_index])
  
  #prepare S star
  AtA_oo<-Cov_reaarrange(t(A)%*%A)$Sigma_oo
  AtA_mo<-Cov_reaarrange(t(A)%*%A)$Sigma_mo
  AtA_om<-Cov_reaarrange(t(A)%*%A)$Sigma_om
  AtA_mm<-Cov_reaarrange(t(A)%*%A)$Sigma_mm
  
  E_star1=result$Y_obs%*%AtA_oo%*%(result$Y_obs)
  E_star2=result$Y_obs%*%AtA_om%*%result$mu_cond
  E_star3=t(result$mu_cond)%*%AtA_mo%*%result$Y_obs
  E_star4=tr(AtA_mm%*%result$Sigma_cond)+t(result$mu_cond)%*%AtA_mm%*%result$mu_cond
  S_star=E_star1+E_star2+E_star3+E_star4
  
  
  Q1<- -W1_fix-t(W1_fix)+2*b_1*t(W1_fix)%*%W1_fix+b_2*(t(W1_fix)%*%W2_fix+t(W2_fix)%*%W1_fix)+b_3*(t(W1_fix)%*%W3_fix+t(W3_fix)%*%W1_fix)
  Q1_oo<-Cov_reaarrange(Q1)$Sigma_oo
  Q1_mo<-Cov_reaarrange(Q1)$Sigma_mo
  Q1_om<-Cov_reaarrange(Q1)$Sigma_om
  Q1_mm<-Cov_reaarrange(Q1)$Sigma_mm
  E_star1=result$Y_obs%*%Q1_oo%*%(result$Y_obs)
  E_star2=result$Y_obs%*%Q1_om%*%result$mu_cond
  E_star3=t(result$mu_cond)%*%Q1_mo%*%result$Y_obs
  E_star4=tr(Q1_mm%*%result$Sigma_cond)+t(result$mu_cond)%*%Q1_mm%*%result$mu_cond
  S_star_Q1=E_star1+E_star2+E_star3+E_star4
  
  
  
  Q2<- -W2_fix-t(W2_fix)+2*b_2*t(W2_fix)%*%W2_fix+b_1*(t(W1_fix)%*%W2_fix+t(W2_fix)%*%W1_fix)+b_3*(t(W2_fix)%*%W3_fix+t(W3_fix)%*%W2_fix)
  Q2_oo<-Cov_reaarrange(Q2)$Sigma_oo
  Q2_mo<-Cov_reaarrange(Q2)$Sigma_mo
  Q2_om<-Cov_reaarrange(Q2)$Sigma_om
  Q2_mm<-Cov_reaarrange(Q2)$Sigma_mm
  E_star1=result$Y_obs%*%Q2_oo%*%(result$Y_obs)
  E_star2=result$Y_obs%*%Q2_om%*%result$mu_cond
  E_star3=t(result$mu_cond)%*%Q2_mo%*%result$Y_obs
  E_star4=tr(Q2_mm%*%result$Sigma_cond)+t(result$mu_cond)%*%Q2_mm%*%result$mu_cond
  S_star_Q2=E_star1+E_star2+E_star3+E_star4
  
  
  Q3<--W3_fix-t(W3_fix)+2*b_3*t(W3_fix)%*%W3_fix+b_1*(t(W1_fix)%*%W3_fix+t(W3_fix)%*%W1_fix)+b_2*(t(W2_fix)%*%W3_fix+t(W3_fix)%*%W2_fix)
  Q3_oo<-Cov_reaarrange(Q3)$Sigma_oo
  Q3_mo<-Cov_reaarrange(Q3)$Sigma_mo
  Q3_om<-Cov_reaarrange(Q3)$Sigma_om
  Q3_mm<-Cov_reaarrange(Q3)$Sigma_mm
  E_star1=result$Y_obs%*%Q3_oo%*%(result$Y_obs)
  E_star2=result$Y_obs%*%Q3_om%*%result$mu_cond
  E_star3=t(result$mu_cond)%*%Q3_mo%*%result$Y_obs
  E_star4=tr(Q3_mm%*%result$Sigma_cond)+t(result$mu_cond)%*%Q3_mm%*%result$mu_cond
  S_star_Q3=E_star1+E_star2+E_star3+E_star4
  
  
  A_rear<-Cov_reaarrange(A)$CovMat
  A_minus<-solve(A)
  
  
  SS1<-t(Y_star)%*%Cov_reaarrange( -Z%*%W1_fix-t(W1_fix)%*%Z+2*b_1*t(W1_fix)%*%Z%*%W1_fix+b_2*(t(W1_fix)%*%Z%*%W2_fix+t(W2_fix)%*%Z%*%W1_fix)+b_3*(t(W1_fix)%*%Z%*%W3_fix+t(W3_fix)%*%Z%*%W1_fix))$Cov%*%Y_star
  SS2<-t(Y_star)%*%Cov_reaarrange( -Z%*%W2_fix-t(W2_fix)%*%Z+2*b_2*t(W2_fix)%*%Z%*%W2_fix+b_1*(t(W2_fix)%*%Z%*%W1_fix+t(W1_fix)%*%Z%*%W2_fix)+b_3*(t(W2_fix)%*%Z%*%W3_fix+t(W3_fix)%*%Z%*%W2_fix))$Cov%*%Y_star
  SS3<-t(Y_star)%*%Cov_reaarrange( -Z%*%W3_fix-t(W3_fix)%*%Z+2*b_3*t(W3_fix)%*%Z%*%W3_fix+b_1*(t(W3_fix)%*%Z%*%W1_fix+t(W1_fix)%*%Z%*%W3_fix)+b_2*(t(W2_fix)%*%Z%*%W3_fix+t(W3_fix)%*%Z%*%W2_fix))$Cov%*%Y_star
  
  SS_star<-S_star - t(Y_star)%*%t(A_rear)%*%X_rear%*%solve(t(X_rear)%*%X_rear)%*%t(X_rear)%*%A_rear%*%Y_star
  
  dell_rho1<--tr(A_minus%*%W1_fix)-(n_edges/2)*(1/( SS_star  ))*(S_star_Q1-SS1)
  dell_rho2<--tr(A_minus%*%W2_fix)-(n_edges/2)*(1/( SS_star  ))*(S_star_Q2-SS2)
  dell_rho3<--tr(A_minus%*%W3_fix)-(n_edges/2)*(1/( SS_star  ))*(S_star_Q3-SS3)
  
  return(c(-dell_rho1,-dell_rho2,-dell_rho3))
}

# Defne the function to create conditional expectations and variances
Cov_reaarrange2<-function(Cov,number){
  
  zero_ind<-number
  nonzero_ind<-1:n_edges
  nonzero_ind<-nonzero_ind[-number]
  n_zero=length(zero_ind)
  n_greater_zero=length(Y)-n_zero
  
  Cov<-Cov[c(nonzero_ind,zero_ind),c(nonzero_ind,zero_ind)] 
  
  out<-list(cens=zero_ind,noncens=nonzero_ind,CovMat=Cov,Sigma_oo=Cov[1:n_greater_zero,1:n_greater_zero],Sigma_om=Cov[1:n_greater_zero,(n_greater_zero+1):length(Y)],Sigma_mo=Cov[(n_greater_zero+1):length(Y),1:n_greater_zero],Sigma_mm=Cov[(n_greater_zero+1):length(Y),(n_greater_zero+1):length(Y)])
  return(out)
}

# start the loop ----
for (r in 1:rep){
  
  # Draw X
  X<-cbind(1,matrix(rnorm(n_edges*(length(beta)-1),mean_x,sigma_x),nrow=n_edges,ncol=length(beta)-1))
  
  # Define the weighting matrices
  W1_fix<-diag(rep(1,n_edges))-A_constructur(rep(0,n_edges),n_nodes,1,0,0)
  W2_fix<-diag(rep(1,n_edges))-A_constructur(rep(0,n_edges),n_nodes,0,1,0)
  W3_fix<-diag(rep(1,n_edges))-A_constructur(rep(0,n_edges),n_nodes,0,0,1)
  
  # Define all parameters
  A<-diag(rep(1,n_edges))-rho_1*W1_fix-rho_2*W2_fix-rho_3*W3_fix
  B<-solve(A)
  mu<-B%*%X%*%beta
  SIGMA<-sigma*B%*%t(B)
  Y_real<-c(rmvnorm(1,mean=mu,sigma = SIGMA))
  Y<-Y_real
  real<-cbind(real,Y_real)
  # prepare censoring
  Y[Y<quantile(Y,censoring)]<-0
  minimum<-min(Y[Y!=0])
  flag2<-which(Y==0)
  
  Y_save<-Y
  Y_cens<-Y
  
  
  cens_index<-which(Y_cens==0)
  noncens_index<-which(Y_cens!=0)
  
  X_rear=rbind(X[which(Y_cens!=0),],X[which(Y_cens==0),])
  Z<-X%*%solve(t(X)%*%X)%*%t(X)
  
  # Define the endogenous effects for the pseudoML starting values
  W1Y_fix<-W1_fix%*%Y
  W2Y_fix<-W2_fix%*%Y
  W3Y_fix<-W3_fix%*%Y
  
  # Calculate the starting values
  parameter_start<-coef(censReg(Y~-1+W1Y_fix+W2Y_fix+W3Y_fix+X))
  
  parameter=parameter_start
  parameter[1:3]<-0
  
  A_est<-diag(rep(1,length(Y)))
  
  
  
  #E-Step
  sigma_est<-parameter[length(parameter)]
  result<-E_step2(A_est,solve(A_est)%*%solve(t(A_est))*sigma_est,Y_cens,parameter[4:(length(parameter)-1)],trunc=T)
  
  plot(Y_real~result$Y_imputed)
  
  # M-Step
  M<-optim(c(0,0,0),Q_grid,gr=Q_grid_gradient,method = "BFGS")  
  A_est<-diag(rep(1,length(Y)))-M$par[1]*W1_fix-M$par[2]*W2_fix-M$par[3]*W3_fix
  parameter<-c(M$par,solve( t(X_rear)%*%X_rear,t(X_rear)%*%Cov_reaarrange(A_est)$CovMat%*%c(result$Y_imputed[noncens_index],result$Y_imputed[cens_index]) ),(1/n_edges)*exp(  (2/n_edges)*( M$value+sum(log(eigen(A_est)$values)))  ))
  
  
  
  diff<-100
  value_old<-parameter
  while(diff>stop){
    
    #E-Step
    sigma_est<-parameter[length(parameter)]
    result<-E_step2(A_est,solve(A_est)%*%solve(t(A_est))*sigma_est,Y_cens,parameter[4:(length(parameter)-1)],trunc=T)
    
    plot(Y_real~result$Y_imputed)
    abline(a=0,b=1)
    # M-Step
    M<-optim(M$par,Q_grid,gr=Q_grid_gradient,method = "BFGS")  
    A_est<-diag(rep(1,length(Y)))-M$par[1]*W1_fix-M$par[2]*W2_fix-M$par[3]*W3_fix
    parameter<-c(M$par,solve( t(X_rear)%*%X_rear,t(X_rear)%*%Cov_reaarrange(A_est)$CovMat%*%c(result$Y_imputed[noncens_index],result$Y_imputed[cens_index]) ),(1/n_edges)*exp(  (2/n_edges)*( M$value+sum(log(eigen(A_est)$values)))  ))
    
    diff<-t(value_old-parameter)%*%(value_old-parameter)
    
    value_old<-parameter
    print(diff)
    
  }
  
  
  
  #-----------------------------------------------------------#
  # Extract the conditional proabilities
  #-----------------------------------------------------------#
  
  
  theta_hat<-parameter
  
  int_positive<-c()
  
  A_inv<-chol2inv(chol(A_est))
  
  # Calculate mu hat
  mu<-A_inv%*%X%*%theta_hat[4:(length(theta_hat)-1)]
  # Calculate sigma hat
  Cov<-A_inv%*%t(A_inv)*theta_hat[length(theta_hat)]
  
  # Define Z_{-ij}
  Z_minus_ij<-Y
  # Fill the empty entries with the conditional expectation
  Z_minus_ij[cens_index]<-mu[cens_index]
  

  
  probs<-c()
  
  for (i in 1:length(Y)){
    out<-Cov_reaarrange2(Cov,i)
    
    mu_o<-mu[out$noncens]
    mu_m<-mu[out$cens]
    Sigma_oo_inv<-    chol2inv(chol(out$Sigma_oo))
    
    mu_cond<-mu_m+out$Sigma_mo%*%Sigma_oo_inv%*%(Z_minus_ij[out$noncens]-mu_o) 
    Sigma_cond<-out$Sigma_mm-out$Sigma_mo%*%Sigma_oo_inv%*%out$Sigma_mo
    
    probs<-c(probs,1-pnorm(minimum,mean=mu_cond,sd=Sigma_cond))
    print(i)
  }
  
  
  plot(probs~Y_real)
  abline(h=0.5)  
  
  
  ##############
  
  
  #-----------------------------------------------------------#
  # Latent Utility Analysis
  #-----------------------------------------------------------#
  
  
  
  probis<-probs
  trade<-Y
  
  n<-n_nodes
  
  Pi<-(probis)
  Gamma<-(Y>0)
  
  roc_t<-roc(response=Gamma, predictor = Pi)
  
  plot(roc_t)
  opt_cutoff<-coords(roc_t,"best")[1]
  
  
  plot(probs~Y_real,col=col)
  abline(h=opt_cutoff)  
  
  share<-c(share,  sum(probis>opt_cutoff))
  share2<-c(share2, sum(probis[flag2]>opt_cutoff))
  
  
  estimated<-cbind(estimated,result$Y_imputed)
  
  parameter_est<-cbind(parameter_est,parameter)
  
  
  A1<-sum(probis[flag2]>opt_cutoff)
  A2<-0.75*n_edges-A1
  TP<-c(TP, 0 )
  FP<-c(FP, A1 )
  FN<-c(FN,0)
  TN<-c(TN,A2)
  
  
}


# Create a plot that shows the recovery of the latent variables ----

library(mgcv)
library(RColorBrewer)
thres=45
mod<-gam(c(real)[estimated<thres]~s(c(estimated)[estimated<thres]))
p<-plot(mod)

pdf("5) Simulation Study/dgp1.pdf", width = 12,height=6)
par(mfrow=c(1,2))




plot(real[,1],estimated[,1],ylim=c(min(estimated),45),xlim=c(min(estimated),45),col="gray",xlab="Y|Y<c",ylab="E[Y|Y<c,X]")

for (i in 1:dim(real)[2]){
  points(real[,i],estimated[,i],col="gray")
}

abline(a=0,b=1,lty=2)

lines(p[[1]]$fit+coef(mod)[1],p[[1]]$x,lwd=2)

#X<-cbind(c(real[estimated<thres]),c(estimated[estimated<thres]))
#hpts <- chull(X)
#hpts <- c(hpts, hpts[1])
#lines(X[hpts, ],lty=2)

## some pretty colors

k <- 8
my.cols <- rev(brewer.pal(k, "Spectral"))

## compute 2D kernel density, see MASS book, pp. 130-131
z <- kde2d(c(real[estimated<thres]),c(estimated[estimated<thres]),n=100)

contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)


diff<-(parameter_est-c(rho_1,rho_2,rho_3,beta,1))[1:3,]
boxplot(t(diff),xaxt="n",ylim=c(-0.2,0.2),xlab="parameter",ylab="estimated-true")
abline(h=0)
axis(1,1:3,c("rho_1","rho_2","rho_3"))
dev.off()
