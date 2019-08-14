#-----------------------------------------------------------#
# Plot the Coefficients against time
#-----------------------------------------------------------#

library(ggplot2)
## Now we combine parameters and bootstrap
# Load in the data
rm(list=ls())

setwd("~/Workflow")


time_frame<-c(1993:2014)

res<-c()
se<-c()
for (i in time_frame){
  load(paste("2) Estimating the Models/standard errors/sd_",i,".RData",sep=""))

  res<-cbind(res,theta_hat)
  se<-cbind(se,sd)
  
  
}


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



library(ggplot2)

ts<-res[1,]
pv<- 2*pnorm(-1*abs(res[1,]/se[1,]),0,1)
Kiu<-res[1,]-2*se[1,]
Kio<-res[1,]+2*se[1,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Rho_1","KIu","KIo", 'pi')
p1<-qplot(Year, Rho_1, data=data) + geom_line(color='black')+ggtitle("Reciprocity")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_1<-p1+geom_abline(intercept =0,slope=0)+ylim(-0.05, 0.15)
###


ts<-res[2,]
pv<- 2*pnorm(-1*abs(res[2,]/se[2,]),0,1)
Kiu<-res[2,]-2*se[2,]
Kio<-res[2,]+2*se[2,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Rho_2","KIu","KIo", 'pi')
p1<-qplot(Year, Rho_2, data=data) + geom_line(color='black')+ggtitle("Exporter Effect")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))


p_2<-p1+geom_abline(intercept =0,slope=0)+ylim(0, 0.8)

ts<-res[3,]
pv<- 2*pnorm(-1*abs(res[3,]/se[3,]),0,1)
Kiu<-res[3,]-2*se[3,]
Kio<-res[3,]+2*se[3,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Rho_3","KIu","KIo", 'pi')
p1<-qplot(Year, Rho_3, data=data) + geom_line(color='black')+ggtitle("Importer Effect")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_3<-p1+geom_abline(intercept =0,slope=0)+ylim(0, 0.8)


ts<-res[4,]
pv<- 2*pnorm(-1*abs(res[4,]/se[4,]),0,1)
Kiu<-res[4,]-2*se[4,]
Kio<-res[4,]+2*se[4,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("Constant")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_4<-p1+geom_abline(intercept =0,slope=0)

ts<-res[5,]
pv<- 2*pnorm(-1*abs(res[5,]/se[5,]),0,1)
Kiu<-res[5,]-2*se[5,]
Kio<-res[5,]+2*se[5,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("log(GDP), Exporter")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_5<-p1+geom_abline(intercept =0,slope=0)+ylim(-0.3, 0.8)


ts<-res[6,]
pv<- 2*pnorm(-1*abs(res[6,]/se[6,]),0,1)
Kiu<-res[6,]-2*se[6,]
Kio<-res[6,]+2*se[6,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("log(GDP), Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_6<-p1+geom_abline(intercept =0,slope=0)+ylim(-0.3, 0.8)



ts<-res[7,]
pv<- 2*pnorm(-1*abs(res[7,]/se[7,]),0,1)
Kiu<-res[7,]-2*se[7,]
Kio<-res[7,]+2*se[7,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("Conflict, Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_7<-p1+geom_abline(intercept =0,slope=0)


ts<-res[8,]
pv<- 2*pnorm(-1*abs(res[8,]/se[8,]),0,1)
Kiu<-res[8,]-2*se[8,]
Kio<-res[8,]+2*se[8,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("Path Dependence, Exporter-Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_8<-p1+geom_abline(intercept =0,slope=0)

ts<-res[9,]
pv<- 2*pnorm(-1*abs(res[9,]/se[9,]),0,1)
Kiu<-res[9,]-2*se[9,]
Kio<-res[9,]+2*se[9,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("Absolute Difference Polity Score, Exporter Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_9<-p1+geom_abline(intercept =0,slope=0)


ts<-res[10,]
pv<- 2*pnorm(-1*abs(res[10,]/se[10,]),0,1)
Kiu<-res[10,]-2*se[10,]
Kio<-res[10,]+2*se[10,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("log(Distance), Exporter-Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_10<-p1+geom_abline(intercept =0,slope=0)

ts<-res[11,]
pv<- 2*pnorm(-1*abs(res[11,]/se[11,]),0,1)
Kiu<-res[11,]-2*se[11,]
Kio<-res[11,]+2*se[11,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("Formal Alliance, Exporter-Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_11<-p1+geom_abline(intercept =0,slope=0)


ts<-res[12,]
pv<- 2*pnorm(-1*abs(res[12,]/se[12,]),0,1)
Kiu<-res[12,]-2*se[12,]
Kio<-res[12,]+2*se[12,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("log(Civilian Weapons), Exporter->Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_12<-p1+geom_abline(intercept =0,slope=0)

ts<-res[13,]
pv<- 2*pnorm(-1*abs(res[13,]/se[13,]),0,1)
Kiu<-res[13,]-2*se[13,]
Kio<-res[13,]+2*se[13,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("Export of MCW, Exporter->Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_13<-p1+geom_abline(intercept =0,slope=0)


ts<-res[14,]
pv<- 2*pnorm(-1*abs(res[14,]/se[14,]),0,1)
Kiu<-res[14,]-2*se[14,]
Kio<-res[14,]+2*se[14,]


year<-time_frame

pv[pv>0.1]<-3
pv[pv>0.05 & pv<0.1]<-2
pv[pv<0.05]<-1

data<-data.frame(year,ts, Kiu,Kio, pv)
colnames(data)<-c("Year","Beta","KIu","KIo", 'pi')
p1<-qplot(Year, Beta, data=data) + geom_line(color='black')+ggtitle("I(Trade of MCW), Exporter->Importer")
p1<-p1+geom_ribbon(aes(ymin=KIu,ymax=KIo), alpha=0.3)
p1<-p1+scale_colour_manual(values=c('1'='springgreen4', '2'='orange1', '3'='red'))
p1<-p1+theme(legend.position='none', axis.text=element_text(size=10),axis.title=element_text(size=11),
             plot.title=element_text(size=14, face='bold'))+ geom_point(aes(colour=factor(pi)))

p_14<-p1+geom_abline(intercept =0,slope=0)



pdf("3) Plotting the Results/coef.pdf", width = 12,height=16)
multiplot(p_5,p_7,p_9,p_8,p_14,p_4,p_2,p_6,p_10,p_11,p_12,p_13,p_1,p_3,cols=2)
dev.off()






