#remove the old stuff
rm(list = ls())
#load some packages
library(statnet)
library(xtable)

# set working directory
setwd("~/Workflow")

# load the data
load("Data.RData")



# find the biggest exporters in the sense that they have most outgoing ties
rS<-c()
for (i in 1:23){
  bin<-SALW[[i]]
  bin[bin>0]<-1 # discretize the net
  rS<-cbind(rS,rowSums(bin)) # calculate rowsums
}


rS_save<-rS # save the year-wise rowSums
rS<-rowSums(rS) # calculate the yearsums of the year-wise rowSums

# find the biggest importers in the sense that they have most ingoing ties
cS<-c()
for (i in 1:23){
  bin<-SALW[[i]]
  bin[bin>0]<-1
  cS<-cbind(cS,colSums(bin))
}

cS_save<-cS # save the year-wise colSums
cS<-rowSums(cS) # calculate the yearsums of the year-wise rowSums

sender_product<-apply(cS_save, 1, prod) # is only >0 of there was at least one export
receiver_product<-apply(rS_save, 1, prod) # is only >0 of there was at least one import

# pre-select those countries we have defined to be existant
select_ex<-which(rowSums(nEX[,(1992-1949):(2014-1949)])==23)

#select those countries that are among the top-seller and top-receiver 
select_rS<-which(rS>quantile(rS[rS>0],0.6)&cS>quantile(cS[cS>0],0.6))


#select_rS<-which(sender_product>0|receiver_product>0)
select<-select_ex[which(select_ex%in%select_rS)]

# Search for isolates
for (i in 1:23){
  net<-SALW[[i]][select,select]
  print(which(rowSums(net)==0&colSums(net)==0))
}

# Luxembourg is an isolate in some years
select<-select[-which(names(select)=="Luxembourg")]

par(mfrow=c(1,1))
for (t in c(1992-1991,1999-1991,2006-1991,2014-1991)){
  pdf(paste("1) Structuring Data for the Model/networks_",t+1991,".pdf",sep=""),height=25,width=10)
  
net<-SALW[[t]][select,select]
cex<- 2.2
net_sel<-net*0.2
cex_nodes<-3.2
colnames(net)<-laenderliste[select,5]
rownames(net)<-laenderliste[select,5]
net_50<-network(net)


plot(net_50,displayisolates=F,displaylabels=T,label.cex=1.4,edge.lwd=cex,vertex.col="gray",edge.col="black",vertex.cex=cex_nodes,label.pos=5,interactive=F)

dev.off()
}

# Now 59 countries remain
length(select)


density<-c()
share_bin_trade<-c()
share_vol_trade<-c()
export_vol<-c()
import_vol<-c()
for (i in 1:23){
  bin<-SALW[[i]]
  share_vol_trade<-c(share_vol_trade,sum(bin[select,select])/sum(bin))
  export_vol<-cbind(export_vol,rowSums(bin[select,select]))
  import_vol<-cbind(import_vol,colSums(bin[select,select]))
  bin[bin>0]<-1
  density<-c(density,sum(bin[select,select])/(length(select)*(length(select-1))))
  share_bin_trade<-c(share_bin_trade,sum(bin[select,select])/sum(bin))
 

}

rowSums(export_vol)[order(rowSums(export_vol),decreasing = T)]
USA<-export_vol[which(names(select)=="United States"),]
DEU<-export_vol[which(names(select)=="Germany"),]
ITA<-export_vol[which(names(select)=="Italy"),]
rest<-colSums(export_vol)-USA-DEU-ITA


d <- data.frame(year=rep(1992:2014,each=4),countries=rep(c("other","USA","DEU","ITA"),23),vol=c(rbind(rest,USA,DEU,ITA)))

p<-ggplot(d, aes(x=year,y=vol,group=countries,fill=countries)) + geom_area( position = "stack")
p_exp<-p+theme(legend.position="bottom")

rowSums(import_vol)[order(rowSums(import_vol),decreasing = T)]
USA<-import_vol[which(names(select)=="United States"),]
SAU<-import_vol[which(names(select)=="Saudi Arabia"),]
DEU<-import_vol[which(names(select)=="Germany"),]
rest<-colSums(import_vol)-USA-DEU-SAU


d <- data.frame(year=rep(1992:2014,each=4),countries=rep(c("other","USA","DEU","SAU"),23),vol=c(rbind(rest,USA,DEU,SAU)))

p<-ggplot(d, aes(x=year,y=vol,group=countries,fill=countries)) + geom_area( position = "stack")
p_imp<-p+theme(legend.position="bottom")

dens_data<-data.frame(t=1992:2014,density)
p_dens<-ggplot(dens_data,aes(x=t,y=density))+geom_line()
p_dens

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} 


pdf("1) Structuring Data for the Model/features.pdf", width = 8,height=4)
multiplot(p_exp,p_dens,cols=2)
dev.off()

par(mfrow=c(1,3))
plot(density,type="l")
plot(share_bin_trade,type="l")
plot(share_vol_trade,type="l")

rm(list=ls()[-which(ls()=="select")])


# Based on the Variable select we now put all covariates of interest in the form needed for regression
load("Data.RData")
save.image("1) Structuring Data for the Model/Data_ready_for_regression.RData")

