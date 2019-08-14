#-----------------------------------------------------------#
# Estimating the Models
#-----------------------------------------------------------#

#remove the old stuff
rm(list = ls())
#load some packages
library(censReg)
library(tmvtnorm)
library(matrixcalc)
library(psych)
library(TruncatedNormal)
library(MASS)


set.seed(123)
setwd("~/Workflow")


t_select<-c(1993:2014)
# This will take a very long time (hours per iteration), think therefore about running the loop in parallel!
for (t in t_select){
load("1) Structuring Data for the Model/Data_ready_for_regression.RData")

# For a demonstatration of the algorithm with a reduced data set one can select a subset with big exporters of
# the 59 countries under study, this goes much faster and provides similar results
#select<-select[c(2,3,5,7,17,18,24,25,35,44,52,57,58)]
#select
# Define Stopping Criteria for the Algorithm
stop<-0.1  
  
  # the corresponding function:
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
  
  # calculate the number of nodes and edges
  n_nodes<-length(select)
  n_edges<-n_nodes*(n_nodes-1)
  
  
  W1_fix<-diag(rep(1,n_edges))-A_constructur(rep(0,n_edges),n_nodes,1,0,0)
  W2_fix<-diag(rep(1,n_edges))-A_constructur(rep(0,n_edges),n_nodes,0,1,0)
  W3_fix<-diag(rep(1,n_edges))-A_constructur(rep(0,n_edges),n_nodes,0,0,1)
  
  # select the response network and apply the logarithm
  net<-log(SALW[[t-1991]])
  net[net==-Inf]<-0
  diag(net)<-0
  # select the relevant countries
  net<-net[select,select]
  
  # select the difference between the full net and the arms net
  # being the sport net
  if (t>1996){sportnet<-(log(1+SALW_plus_sport[[t-1991]]-SALW[[t-1991]])+log(1+SALW_plus_sport[[t-1991-1]]-SALW[[t-1991-1]])+log(1+SALW_plus_sport[[t-1991-2]]-SALW[[t-1991-2]])+log(1+SALW_plus_sport[[t-1991-3]]-SALW[[t-1991-3]])+log(1+SALW_plus_sport[[t-1991-4]]-SALW[[t-1991-4]]))/5  }
  if (t==1996){sportnet<- (log(1+SALW_plus_sport[[t-1991]]-SALW[[t-1991]])+log(1+SALW_plus_sport[[t-1991-1]]-SALW[[t-1991-1]])+log(1+SALW_plus_sport[[t-1991-2]]-SALW[[t-1991-2]])+log(1+SALW_plus_sport[[t-1991-3]]-SALW[[t-1991-3]]))/4 }
  if (t==1995){sportnet<- (log(1+SALW_plus_sport[[t-1991]]-SALW[[t-1991]])+log(1+SALW_plus_sport[[t-1991-1]]-SALW[[t-1991-1]])+log(1+SALW_plus_sport[[t-1991-2]]-SALW[[t-1991-2]]))/3 }
  if (t==1994){sportnet<- (log(1+SALW_plus_sport[[t-1991]]-SALW[[t-1991]])+log(1+SALW_plus_sport[[t-1991-1]]-SALW[[t-1991-1]]))/2 }
  if (t==1993){sportnet<- log(1+SALW_plus_sport[[t-1991]]-SALW[[t-1991]]) } 
  sportnet[sportnet==-Inf]<-0
  diag(sportnet)<-0
  # select the relevant countries
  sportnet<-sportnet[select,select]  
  
  # Selected the lagged network and truncated the respective moving average
  if (t>1996){lagnet<-(log(SALW[[t-1991-1]])+log(SALW[[t-1991-2]])+log(SALW[[t-1991-3]])+log(SALW[[t-1991-4]])+log(SALW[[t-1991-5]]))/5  }
  if (t==1996){lagnet<- (log(SALW[[t-1991-1]])+log(SALW[[t-1991-2]])+log(SALW[[t-1991-3]])+log(SALW[[t-1991-4]]))/4 }
  if (t==1995){lagnet<- (log(SALW[[t-1991-1]])+log(SALW[[t-1991-2]])+log(SALW[[t-1991-3]]))/3 }
  if (t==1994){lagnet<- (log(SALW[[t-1991-1]])+log(SALW[[t-1991-2]]))/2 }
  if (t==1993){lagnet<- log(SALW[[t-1991-1]]) }  
  lagnet[lagnet==-Inf]<-0
  diag(lagnet)<-0
  # select the relevant countries
  lagnet<-lagnet[select,select]
  
  
  # Include the arms trade network as the log of the moving average
  mcwnet<-log(1+amk[[t-1949]]+amk[[t-1950]]+amk[[t-1951]]+amk[[t-1952]]+amk[[t-1953]])
  mcwnet<-mcwnet[select,select]
  diag(mcwnet)<-0
  
  bin_mcwnet<-mcwnet
  bin_mcwnet[bin_mcwnet!=0]<-1
  
  # log(GDP)
  lgdp_i<-log(GDPimp[select,t-1989])
  lgdp_j<-log(GDPimp[select,t-1989])
  
  # Polity Scores
  polity_i<-autopolity[select,t-1989]
  polity_j<-autopolity[select,t-1989]
  poldiff_ij<-outer(polity_i,polity_j,FUN = function(x,y){abs(x-y)})
  
  # Formal alliances
  formal_alliance_ij<-daml[[t-1989]][select,select] 
  
  # logarithmic geographic distance
  DIST_ij<-log(cdist)
  diag( DIST_ij)<-0
  DIST_ij<-  DIST_ij[select,select]
  
  # intrastate conflicts
  intra_conf_i<-conf[select,t-1989]
  intra_conf_j<-conf[select,t-1989]
  
  # save the relevant net and bring the covariates in a newtwork form
  
  # response network
  A<-net
  # censored network
  A_cens<-net
  # GDPi
  GDP_i<-matrix(rep(lgdp_i,dim(A)[1]),ncol=dim(A)[1],nrow=dim(A)[1])
  # GDPj
  GDP_j<-t(matrix(rep(lgdp_i,dim(A)[1]),ncol=dim(A)[1],nrow=dim(A)[1]))
  # Intrastateconflict j
  CONF_j<-t(matrix(rep(intra_conf_j,dim(A)[1]),ncol=dim(A)[1],nrow=dim(A)[1]))
  
  
  
  # order the covariates for a regression framework
  lgdp_i<-c()
  lgdp_j<-c()
  conf_j<-c()
  lagnet_ij<-c()
  polity_ij<-c()
  geodist_ij<-c()
  alliance_ij<-c()
  sport_ij<-c()
  mcw_ij<-c()
  bin_mcw_ij<-c()
  
  for (i in 1: n_nodes){
    for (j in 1:n_nodes){
      if (i!=j){
        lgdp_i<- c(lgdp_i,GDP_i[i,j])
        lgdp_j<- c(lgdp_j,GDP_j[i,j])
        conf_j<- c(conf_j,CONF_j[i,j])
        lagnet_ij<- c(lagnet_ij,lagnet[i,j])
        polity_ij<- c(polity_ij,poldiff_ij[i,j])
        geodist_ij<- c(geodist_ij,DIST_ij[i,j])
        alliance_ij<- c(alliance_ij,formal_alliance_ij[i,j])
        sport_ij<- c(sport_ij,sportnet[i,j])
        mcw_ij<-c(mcw_ij,mcwnet[i,j])
        bin_mcw_ij<-c(bin_mcw_ij,bin_mcwnet[i,j])
        
      }
      
    }
  }
  
  diag<-seq(from=1,to=n_nodes^2,by=n_nodes+1)
  response<-c(t(A))
  Y<-response[-diag]
  minimum<-min(Y[Y>0])
  
  
  X<-cbind(1,lgdp_i,lgdp_j, conf_j,lagnet_ij, polity_ij,geodist_ij, alliance_ij,sport_ij,mcw_ij,bin_mcw_ij)
  
  
  # Define the Covariancve-rearranging function
  Cov_reaarrange<-function(Cov){
    
    zero_ind<-which(Y_cens==0)
    nonzero_ind<-which(Y_cens!=0)
    n_zero=length(zero_ind)
    n_greater_zero=length(Y)-n_zero
    
    Cov<-Cov[c(nonzero_ind,zero_ind),c(nonzero_ind,zero_ind)] 
    
    out<-list(cens=zero_ind,noncens=nonzero_ind,CovMat=Cov,Sigma_oo=Cov[1:n_greater_zero,1:n_greater_zero],Sigma_om=Cov[1:n_greater_zero,(n_greater_zero+1):length(Y)],Sigma_mo=Cov[(n_greater_zero+1):length(Y),1:n_greater_zero],Sigma_mm=Cov[(n_greater_zero+1):length(Y),(n_greater_zero+1):length(Y)])
    return(out)
  }
  
  
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
  
  
  #E-Step
  sigma_est<-parameter[length(parameter)]
  result<-E_step2(A_est,solve(A_est)%*%solve(t(A_est))*sigma_est,Y_cens,parameter[4:(length(parameter)-1)],trunc=T)
  
  plot(Y_save~result$Y_imputed)
  
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
    
    plot(Y_save~result$Y_imputed)
    
    # M-Step
    M<-optim(M$par,Q_grid,gr=Q_grid_gradient,method = "BFGS")  
    A_est<-diag(rep(1,length(Y)))-M$par[1]*W1_fix-M$par[2]*W2_fix-M$par[3]*W3_fix
    parameter<-c(M$par,solve( t(X_rear)%*%X_rear,t(X_rear)%*%Cov_reaarrange(A_est)$CovMat%*%c(result$Y_imputed[noncens_index],result$Y_imputed[cens_index]) ),(1/n_edges)*exp(  (2/n_edges)*( M$value+sum(log(eigen(A_est)$values)))  ))
    
    diff<-t(value_old-parameter)%*%(value_old-parameter)
    
    value_old<-parameter
    print(diff)
    
  }
  theta_hat<-parameter
  
  rm(list=ls()[-c(which(ls()=="theta_hat"),which(ls()=="n_edges"),which(ls()=="minimum"),which(ls()=="A_est"),which(ls()=="noncens_index"),which(ls()=="cens_index"),which(ls()=="X_rear"),which(ls()=="X"),which(ls()=="Cov_reaarrange"),which(ls()=="W1_fix"),which(ls()=="W2_fix"),which(ls()=="W3_fix"),which(ls()=="Y"),which(ls()=="Y_cens"),which(ls()=="t"),which(ls()=="t_select"))])
  
  
  save.image(paste("2) Estimating the Models/results/coef_",t,".RData",sep=""))
}




