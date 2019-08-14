#-----------------------------------------------------------#
# Latent Utility Analysis
#-----------------------------------------------------------#

rm(list=ls())

library(statnet)
library(igraph)
library(ggplot2)
library(pROC)

setwd("~/Workflow")
load("1) Structuring Data for the Model/Data_ready_for_regression.RData")


t_select<-c(1993:2014)
probis<-c()
trade<-c()
for (t in t_select){
  load(paste("2) Estimating the Models/probs/probs_",t,".RData",sep=""))
  probis<-cbind(probis,probs)
  trade<-cbind(trade,Y)
}

n<-length(select)

t_select<-c(1993:2014)

opt_cutoff<-c()
for (t in 1:length(t_select)){
  Pi<-(probis[,t])
  Gamma<-(trade[,t]>0)
  
  roc_t<-roc(response=Gamma, predictor = Pi)
  
  plot(roc_t)
  opt_cutoff<-c(opt_cutoff,coords(roc_t,"best")[1])
}

Pi<-c()

for (t in 1:length(t_select)){
  Pi<-cbind(Pi,(probis[,t]>opt_cutoff[t]))
  
}

Gamma<-(trade>0)
Omega<-Pi-Gamma

Omega_plus<-(Omega>0)
Omega_minus<-(Omega<0)

prob_bin_plus<-list()
prob_bin_minus<-list()

for (t in 1:length(t_select)){
  mat_plus<-matrix(0,ncol=n,nrow=n)
  mat_minus<-matrix(0,ncol=n,nrow=n)
  
  rownames(mat_plus)<-laenderliste[which(laenderliste[,1]%in%names(select)),5]
  colnames(mat_plus)<-laenderliste[which(laenderliste[,1]%in%names(select)),5]
  
  rownames(mat_minus)<-laenderliste[which(laenderliste[,1]%in%names(select)),5]
  colnames(mat_minus)<-laenderliste[which(laenderliste[,1]%in%names(select)),5]
  
  up<-1
  for (i in 1:n){
    for (j in 1:n){
      if (i != j){
        mat_plus[i,j]<-Omega_plus[up,t]
        mat_minus[i,j]<-Omega_minus[up,t]
        up<-up+1
      }
    }
  }
  prob_bin_plus[[t]]<-mat_plus
  prob_bin_minus[[t]]<-mat_minus
  
}



selfun<-function(prob_bin=prob_bin_plus,period, threshold,upper_threshold=30,MA2=F,MA3=F,MA4=F,outdegree=0,pos=5,iso=F,interactive_switch=F,retmat=F,scale=0.4,e.scale=1.1){
  
  from<-which(t_select==period[1])
  to<-which(t_select==period[length(period)])
  
  A_out<-matrix(0,ncol = n,nrow = n)
  if (MA2==F) {
    for (t in from:to){
      A_out<-A_out+prob_bin[[t]]
    }
    if (retmat==T){return(A_out)}
    cex_nodes1<-0.4+sqrt(rowSums(A_out))*scale
    cex_nodes2<-0.4+sqrt(colSums(A_out))*scale
    cex_nodes<-cex_nodes1
    cex_nodes[which(cex_nodes2>cex_nodes1)]<-cex_nodes2[which(cex_nodes2>cex_nodes1)]
    colors<-rep("green",dim(A)[1])
    colors[which(cex_nodes2>cex_nodes1)]<-"red"
    A_out[A_out<threshold]<-0
    A_out[A_out>upper_threshold]<-0
    
    plot(network(A_out),vertex.col=colors,edge.col="gray",vertex.cex=cex_nodes,label.col="black",displayisolates=iso,displaylabels=T,label.pos=pos,interactive=interactive_switch,edge.lwd=(A_out-threshold)*e.scale,usecurve=F)
  }
  
  if (MA2==T) {
    for (t in from:(to-1)){
      A_out<-A_out+prob_bin[[t]]*prob_bin[[t+1]]
    }
    
    
    plot(network(A_out),displayisolates=iso,displaylabels=T,label.pos=5,interactive=interactive_switch)
  }
  
  if (MA3==T) {
    for (t in from:(to-2)){
      A_out<-A_out+prob_bin[[t]]*prob_bin[[t+1]]*prob_bin[[t+2]]
    }
    
    
    plot(network(A_out),displayisolates=F,displaylabels=T,label.pos=5)
  }
  if (MA4==T) {
    for (t in from:(to-3)){
      A_out<-A_out+prob_bin[[t]]*prob_bin[[t+1]]*prob_bin[[t+2]]*prob_bin[[t+3]]
    }
    
    
    plot(network(A_out),displayisolates=F,displaylabels=T,label.pos=5)
  }
  
  
  
  if (outdegree>0) {
    for (t in from:to){
      A_out<-A_out+prob_bin[[t]]
    }
    kill<-which(rowSums(A_out)<outdegree )
    
    
    A_out[kill,]<-0
    plot(network(A_out),displayisolates=F,displaylabels=T,label.pos=5)
  }
  
}

Omega_plus<-selfun(prob_bin=prob_bin_plus,period=t_select,threshold=10,iso=T,e.scale=2.5,scale=0.25,retmat = T)
Omega_minus<-selfun(prob_bin=prob_bin_minus,period=t_select,threshold=10,iso=T,e.scale=2.5,scale=0.25,retmat = T)

name_comb<-c()
for (i in 1:n){
  for (j in 1:n){
    
    name_comb<-c(name_comb,paste(rownames(Omega_plus)[i],"-",colnames(Omega_plus)[j],sep=""))
    
  }
}

Omega_plus<-selfun(prob_bin=prob_bin_plus,period=t_select,threshold=10,iso=T,e.scale=2.5,scale=0.25,retmat = T)
Omega_minus<-selfun(prob_bin=prob_bin_minus,period=t_select,threshold=10,iso=T,e.scale=2.5,scale=0.25,retmat = T)

vecOmega_plus<-t(Omega_plus[Omega_plus>0])
vecOmega_minus<-t(Omega_minus[Omega_minus>0])

pdf("4) Latent Utility Analysis/uncover.pdf",height=8,width=11)

plot(table(vecOmega_minus),lwd=1,type="b",pch="-",xlim=c(1,23),xaxt="n",cex=1.5,ylim=c(-250,500),col="darkgray",ylab = "Frequency",xlab="Value of Omega+ and Omega-",yaxt="n")
lines(table(vecOmega_plus),lwd=1,type="b",pch="+",cex=1.5)
axis(1,1:22,1:22)
axis(2, seq(from=0,to=500,by=50), seq(from=0,to=500,by=50))
abline(h=0,col="lightgray",lty=2)

vecOmega_plus<-t(Omega_plus)
vecOmega_minus<-t(Omega_minus)

object<-which(vecOmega_plus==21)
leng<-length(object)
leng

text(x=22,y=-50,name_comb[object[1]]) 
segments(x0=22, y0=-50, x1=21, y1=3,
         col = par("fg"), lty = par("lty"), xpd = FALSE)

text(x=22,y=-70,name_comb[object[2]]) 
text(x=22,y=-90,name_comb[object[3]]) 
text(x=22,y=-110,name_comb[object[4]]) 


object<-which(vecOmega_plus==20)
leng<-length(object)
leng

text(x=21,y=50,name_comb[object[1]]) 
text(x=21,y=70,name_comb[object[2]])
text(x=21,y=90,name_comb[object[3]])
text(x=21,y=110,name_comb[object[4]])
text(x=21,y=130,name_comb[object[5]])

segments(x0=21, y0=50, x1=20, y1=2,
         col = par("fg"), lty = par("lty"), xpd = FALSE)


object<-which(vecOmega_plus==19)
leng<-length(object)
leng

text(x=19,y=-50,name_comb[object[1]]) 
text(x=19,y=-70,name_comb[object[2]]) 
text(x=19,y=-90,name_comb[object[3]]) 
text(x=19,y=-110,name_comb[object[4]]) 
text(x=19,y=-130,name_comb[object[5]]) 

segments(x0=19, y0=-50, x1=19, y1=10,
         col = par("fg"), lty = par("lty"), xpd = FALSE)

object<-which(vecOmega_plus==18)
leng<-length(object)
leng

text(x=18,y=50,name_comb[object[1]]) 
text(x=18,y=70,name_comb[object[2]]) 
text(x=18,y=90,name_comb[object[3]]) 
text(x=18,y=110,name_comb[object[4]]) 
text(x=18,y=130,name_comb[object[5]]) 
text(x=18,y=150,name_comb[object[6]]) 
text(x=18,y=170,name_comb[object[7]]) 
text(x=18,y=190,name_comb[object[8]]) 
text(x=18,y=210,name_comb[object[9]]) 
text(x=18,y=230,name_comb[object[10]]) 

segments(x0=18, y0=50, x1=18, y1=10,
         col = par("fg"), lty = par("lty"), xpd = FALSE)


object<-which(vecOmega_plus==17)
leng<-length(object)
leng

text(x=16.5,y=-50,name_comb[object[1]]) 
text(x=16.5,y=-70,name_comb[object[2]]) 
text(x=16.5,y=-90,name_comb[object[3]]) 
text(x=16.5,y=-110,name_comb[object[4]]) 
text(x=16.5,y=-130,name_comb[object[5]]) 
text(x=16.5,y=-150,name_comb[object[6]]) 
text(x=16.5,y=-170,name_comb[object[7]]) 
text(x=16.5,y=-190,name_comb[object[8]]) 
text(x=16.5,y=-210,name_comb[object[9]]) 
text(x=16.5,y=-230,name_comb[object[10]]) 
text(x=16.5,y=-250,name_comb[object[11]]) 
segments(x0=16.5, y0=-50, x1=17, y1=7,
         col = par("fg"), lty = par("lty"), xpd = FALSE)

object<-which(vecOmega_plus==16)
leng<-length(object)
leng

text(x=16,y=50,name_comb[object[1]]) 
text(x=16,y=70,name_comb[object[2]]) 
text(x=16,y=90,name_comb[object[3]]) 
text(x=16,y=110,name_comb[object[4]]) 
text(x=16,y=130,name_comb[object[5]]) 
text(x=16,y=150,name_comb[object[6]]) 
text(x=16,y=170,name_comb[object[7]]) 
text(x=16,y=190,name_comb[object[8]]) 
text(x=16,y=210,name_comb[object[9]]) 
text(x=16,y=230,name_comb[object[10]]) 
text(x=16,y=250,name_comb[object[11]]) 
text(x=16,y=270,name_comb[object[12]]) 
text(x=16,y=290,name_comb[object[13]]) 
text(x=16,y=310,name_comb[object[14]]) 
text(x=16,y=330,name_comb[object[15]]) 
text(x=16,y=350,name_comb[object[16]]) 
text(x=16,y=370,name_comb[object[17]]) 
text(x=16,y=390,name_comb[object[18]]) 
text(x=16,y=410,name_comb[object[19]]) 
text(x=16,y=430,name_comb[object[20]]) 
segments(x0=16, y0=50, x1=16, y1=11,
         col = par("fg"), lty = par("lty"), xpd = FALSE)


object<-which(vecOmega_minus==15)
leng<-length(object)
leng
text(x=15,y=-20,name_comb[object[1]],col="darkgray") 
segments(x0=15, y0=-20, x1=15, y1=3, lty = par("lty"), xpd = FALSE,col="darkgray")


object<-which(vecOmega_minus==14)
leng<-length(object)
leng
text(x=14,y=-90,name_comb[object[1]],col="darkgray") 
text(x=14,y=-110,name_comb[object[2]],col="darkgray") 
text(x=14,y=-130,name_comb[object[3]],col="darkgray") 
segments(x0=14, y0=-90, x1=14, y1=3, lty = par("lty"), xpd = FALSE,col="darkgray")


object<-which(vecOmega_minus==13)
leng<-length(object)
leng
text(x=13,y=70,name_comb[object[1]],col="darkgray") 
text(x=13,y=90,name_comb[object[2]],col="darkgray") 

segments(x0=13, y0=50, x1=13, y1=4, lty = par("lty"), xpd = FALSE,col="darkgray")



object<-which(vecOmega_minus==12)
leng<-length(object)
leng
text(x=12,y=-50,name_comb[object[1]],col="darkgray") 
text(x=12,y=-70,name_comb[object[2]],col="darkgray") 
text(x=12,y=-90,name_comb[object[3]],col="darkgray") 
text(x=12,y=-110,name_comb[object[4]],col="darkgray") 
text(x=12,y=-130,name_comb[object[5]],col="darkgray") 
text(x=12,y=-150,name_comb[object[6]],col="darkgray") 
text(x=12,y=-170,name_comb[object[7]],col="darkgray") 
text(x=12,y=-190,name_comb[object[8]],col="darkgray") 

segments(x0=12, y0=-50, x1=12, y1=3, lty = par("lty"), xpd = FALSE,col="darkgray")

object<-which(vecOmega_minus==11)
leng<-length(object)
leng
text(x=11,y=110,name_comb[object[1]],col="darkgray") 
text(x=11,y=130,name_comb[object[2]],col="darkgray") 
text(x=11,y=150,name_comb[object[3]],col="darkgray") 
text(x=11,y=170,name_comb[object[4]],col="darkgray") 
text(x=11,y=190,name_comb[object[5]],col="darkgray") 
text(x=11,y=210,name_comb[object[6]],col="darkgray") 
text(x=11,y=230,name_comb[object[7]],col="darkgray") 
text(x=11,y=250,name_comb[object[8]],col="darkgray") 
text(x=11,y=270,name_comb[object[9]],col="darkgray") 
text(x=11,y=290,name_comb[object[10]],col="darkgray") 
text(x=11,y=310,name_comb[object[11]],col="darkgray") 
text(x=11,y=330,name_comb[object[12]],col="darkgray") 
text(x=11,y=350,name_comb[object[13]],col="darkgray") 
text(x=11,y=370,name_comb[object[14]],col="darkgray") 
text(x=11,y=390,name_comb[object[15]],col="darkgray") 

segments(x0=11, y0=110, x1=11, y1=16, lty = par("lty"), xpd = FALSE,col="darkgray")

object<-which(vecOmega_minus==10)
leng<-length(object)
leng
text(x=8,y=90,name_comb[object[1]],col="darkgray") 
text(x=8,y=110,name_comb[object[2]],col="darkgray") 
text(x=8,y=130,name_comb[object[3]],col="darkgray") 
text(x=8,y=150,name_comb[object[4]],col="darkgray") 
text(x=8,y=170,name_comb[object[5]],col="darkgray") 
text(x=8,y=190,name_comb[object[6]],col="darkgray") 
text(x=8,y=210,name_comb[object[7]],col="darkgray") 
text(x=8,y=230,name_comb[object[8]],col="darkgray") 
text(x=8,y=250,name_comb[object[9]],col="darkgray") 
text(x=8,y=270,name_comb[object[10]],col="darkgray") 
text(x=8,y=290,name_comb[object[11]],col="darkgray") 
text(x=8,y=310,name_comb[object[12]],col="darkgray") 
text(x=8,y=330,name_comb[object[13]],col="darkgray") 
text(x=8,y=350,name_comb[object[14]],col="darkgray") 
text(x=8,y=370,name_comb[object[15]],col="darkgray") 
text(x=8,y=390,name_comb[object[16]],col="darkgray") 
text(x=8,y=410,name_comb[object[17]],col="darkgray") 
segments(x0=8, y0=90, x1=10, y1=16, lty = par("lty"), xpd = FALSE,col="darkgray")



dev.off()

n_edges<-length(select)*(length(select)-1)

n_plus<-c()
n_minus<-c()
for (t in 1:length(t_select)){
  n_plus<-c(n_plus,sum(prob_bin_plus[[t]])/n_edges)
  n_minus<-c(n_minus,sum(prob_bin_minus[[t]])/n_edges) 
}
density<-n_plus
dens_data<-data.frame(t=t_select,density,n_minus)


pdf("4) Latent Utility Analysis/unobserved_denstiy.pdf", width = 6,height=3)

ggplot(dens_data, aes(t)) + 
  geom_line(aes(y = density, colour = "Omega+")) + 
  geom_line(aes(y = n_minus, colour = "Omega-"))
dev.off()



plotter<-function(tframe=1:length(t_select),prob_bin=prob_bin_plus,sel="eigen",direction=F,main){
  e_c<-c()
  
  for (t in tframe){
    if (sel=="eigen"){
      e_c<-cbind(e_c,eigen_centrality(graph_from_adjacency_matrix(prob_bin[[t]]), directed =direction, scale = TRUE, weights = NULL,
                                      options = arpack_defaults)$vector)
    }
    if (sel=="eigen2"){
      e_c<-cbind(e_c,svd(prob_bin[[t]])$d)
    }
    
    
    
    if (sel=="transitivity"){
      e_c<-cbind(e_c,transitivity(graph_from_adjacency_matrix(prob_bin[[t]]),type="local",
                                  isolates = "zero",vids = T))
    }
    
    if (sel=="outdegree"){
      e_c<-cbind(e_c,rowSums(prob_bin[[t]])/max(rowSums(prob_bin[[t]])))
    }
    
    if (sel=="indegree"){
      e_c<-cbind(e_c,colSums(prob_bin[[t]])/max(colSums(prob_bin[[t]])))
    }
    
  }
  row.names(e_c)<-laenderliste[select,5]
  e_c<-e_c[order(apply(e_c,1,median)),]
  boxplot(t(e_c),las=2,main=main)
}

pdf("4) Latent Utility Analysis/unobserved_features1.pdf", width = 16,height=10)
par(mfrow=c(3,1))
par(cex.lab=2.3)
par(cex.axis=2.3)
par(cex.main=2.3)
plotter(prob_bin=prob_bin_plus,sel="eigen",direction = F,main="Eigencentrality Omega+")
plotter(prob_bin=prob_bin_plus,sel="outdegree",main="Outdegree Omega+")
plotter(prob_bin=prob_bin_plus,sel="indegree",main="Indegree Omega+")
dev.off()

pdf("4) Latent Utility Analysis/unobserved_features2.pdf", width = 16,height=10)
par(mfrow=c(3,1))
par(cex.lab=2.3)
par(cex.axis=2.3)
par(cex.main=2.3)
plotter(prob_bin=prob_bin_minus,sel="eigen",direction = F,main="Eigencentrality Omega-")
plotter(prob_bin=prob_bin_minus,sel="outdegree",main="Outdegree Omega-")
plotter(prob_bin=prob_bin_minus,sel="indegree",main="Indegree Omega-")
dev.off()

pdf("4) Latent Utility Analysis/unobserved_features3.pdf", width = 18,height=10)
par(mfrow=c(2,1))
par(cex.lab=2)
par(cex.axis=2)
par(cex.main=2)
plotter(prob_bin=prob_bin_plus,sel="eigen",direction = F,main="Eigencentrality Omega+")
plotter(prob_bin=prob_bin_minus,sel="eigen",direction = F,main="Eigencentrality Omega-")
dev.off()




