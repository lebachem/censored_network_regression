#-----------------------------------------------------------#
# Simulation Study
#-----------------------------------------------------------#

rm(list=ls())

library(censReg)
library(tmvtnorm)
library(matrixcalc)
library(psych)
library(TruncatedNormal)
library(MASS)
library(Matrix)
library(mvtnorm)
library(statnet)
library(igraph)
library(ggplot2)
library(pROC)

# set working directory
setwd("~/Workflow")


# Exract information from DGP1 ----
load("5) Simulation Study/data_dgp1.RData")

# FPR ----
FPR_dgp1<-FP/(0.75*n_edges)


# Extract information from DGP2 ----
load("5) Simulation Study/data_dgp2.RData")
FPR_dgp2<-FP/(0.75*n_edges)
TPR_dgp2<-TP/(0.1*n_edges)
FDR_dgp2<-FP/(FP+TP)



pdf("5) Simulation Study/sim_res.pdf", width = 9,height=4.5)
par(mar=c(2,4,2,2), mfrow=c(1,2),
    oma = c(1, 1, 1, 1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
boxplot(cbind(FPR_dgp1,FPR_dgp2),ylim=c(0,0.1),ylab="False Positive Rate (FPR)",xaxt="n")
axis(1,1:2,c("DGP1","DGP2"))
#dev.off()

boxplot(TPR_dgp2,ylim=c(0,1),ylab="True Positive Rate (TPR)",xaxt="n")
axis(1,1,c("DGP2"))
boxplot(FDR_dgp2,ylim=c(0,1),ylab="False Discovery Rate (FDR)",xaxt="n")
axis(1,1,c("DGP2"))
dev.off()



