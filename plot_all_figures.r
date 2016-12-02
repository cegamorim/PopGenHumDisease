#Eduardo Amorim | guerraamorim@gmail.com | cg2827@columbia.edu
#Nov 29, 2016
#Amorim et al. (2017) - The Population Genetics of Human Disease: The Case of Recessive Lethal Mutations

##################################################################################################
##################################################################################################
#############        Load external packages needed to use color for violin plots     #############
##################################################################################################
##################################################################################################

library(sm)
library(vioplot)
source("colorVioplot.r")
library(plotrix)

##################################################################################################
##################################################################################################
#################        Load data (all available as supplementary material)     #################
##################################################################################################
##################################################################################################

## User may save the data in a directory such as ../FinalData/

## Empirical data:
counsyl<-read.table("../FinalData/counsyl_exac_matched.tab", header=TRUE)
ExAC_noGCI<-read.table("../FinalData/ExAC_final373_filtered_GCIfilter.tab", header=TRUE)
ExAC<-read.table("../FinalData/ExAC_final385_filtered.tab.txt", header=TRUE)

## Deleterious allele frequencies from simulations of single site:
sims_EffectHet<-read.table("../FinalData/sims_currAF_Eur_varU_averageRate_effectHets.tab.txt", header=FALSE)
sims_AverU<-read.table("../FinalData/sims_currAF_Eur_varU_averageRate.tab.txt", header=FALSE)
sims_fourTypes<-read.table("../FinalData/sims_currAF_Eur_varU.tab", header=TRUE)
constSize<-read.table("../FinalData/constPopSize.txt",header=TRUE)
sim1M<-read.table("../FinalData/1M_eurDistrib.txt",header=FALSE)
sim2M<-read.table("../FinalData/2M_eurDistrib.txt",header=FALSE)
sim5M<-read.table("../FinalData/5M_eurDistrib.txt",header=FALSE)
## Simulations by gene (equivalent to the sum of deleterious allele frequency per gene)
eurbygene<-read.table("../FinalData/sims_euroByGene.tab",header=TRUE)

#################################################################################################
#################################################################################################
################################      Editing   datasets     ####################################
#################################################################################################
#################################################################################################

# Simulations:
# Replacing 0 per 1e-7 
# Average U - Tennessen model
sims_AverULOG<-replace(sims_AverU,sims_AverU==0,1e-7)
# Average U - Tennessen model with effect in hets
sims_EffectHetLOG<-replace(sims_EffectHet,sims_EffectHet==0,1e-7)             
# Average U - Constant population size
constSizeLOG<-replace(constSize,constSize==0,1e-7)
# Average U - More intense growth
sim1MLOG<-replace(sim1M,sim1M==0,1e-7)
sim2MLOG<-replace(sim2M,sim2M==0,1e-7)
sim5MLOG<-replace(sim5M,sim5M==0,1e-7)
#Subsample
sims_AverUPOI<-as.data.frame(apply(sims_AverU, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sims_AverUPOI)<-c("V1")
sim1MPOI<-as.data.frame(apply(sim1M, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sim1MPOI)<-c("V1")
sim2MPOI<-as.data.frame(apply(sim2M, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sim2MPOI)<-c("V1")
sim5MPOI<-as.data.frame(apply(sim5M, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sim5MPOI)<-c("V1")

# 4 types of mutations
# First subsample to match ExAC's average sample size (32,881)
simsPOI<-as.data.frame(apply(sims_fourTypes[,1:4], c(1,2), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
simsPOILOG<-replace(simsPOI,simsPOI==0,1e-7)                 
# Gene level analysis
eurbygenePOI<-as.data.frame(apply(eurbygene[,1:32], c(1,2), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
eurbygenePOILOG<-replace(eurbygenePOI,eurbygenePOI==0,1e-7)

# Empirical
# Original: counsyl & ExAC_noGCI & ExAC
# Replacing q = 0 to 1e-7:
counsylLOG<-replace(counsyl,counsyl==0,1e-7)
ExAC_noGCILOG<-replace(ExAC_noGCI,ExAC_noGCI==0,1e-7)
ExACLOG<-replace(ExAC,ExAC==0,1e-7)

#Subset empirical data by mutation type
XnonCpGtv <- subset(ExAC_noGCI,type=="nonCpGtv")
XnonCpGti <- subset(ExAC_noGCI,type=="nonCpGti")
XCpGtv <- subset(ExAC_noGCI,type=="CpGtv")
XCpGti <- subset(ExAC_noGCI,type=="CpGti")
XnonCpGtvLOG <- subset(ExAC_noGCILOG,type=="nonCpGtv")
XnonCpGtiLOG <- subset(ExAC_noGCILOG,type=="nonCpGti")
XCpGtvLOG <- subset(ExAC_noGCILOG,type=="CpGtv")
XCpGtiLOG <- subset(ExAC_noGCILOG,type=="CpGti")

#Set predictions based on mut-selection balance models for lethal recessive alleles
aver_u=1.2e-8
FIN=function(N){aver_u*sqrt(2*pi*N)} #Here, N is the number of diploid individuals, so as the input used for the c++ code used for simulations of constant pop size
INF=sqrt(aver_u)
FIN20K=FIN(20000)
FINTENN=FIN(512000)

############################################################################################
################################        Results 1:          ################################
################################  Compare Counsyl and ExAC  ################################
############################################################################################

# Wilcoxon Signed-Rank Test (nonparametric equivalent of paired t-test)
wilcox.test(counsyl$freqEur_Counsyl,counsyl$EuroExAC,paired = TRUE) # p-value = 0.3411 (can't reject the hypothesis of being from the same population)

############
## Fig S1 ##
############

ExAC_Obs<-ExAC$eur_freq[which(ExAC$eur_freq != 0)]
ExAC_0<-(((length(ExAC$eur_freq)-length(ExAC_Obs))/length(ExAC$eur_freq))*1)-6.5
CX_Obs<-counsyl$EuroExAC[which(counsyl$EuroExAC != 0)]
CX_0<-(((length(counsyl$EuroExAC)-length(CX_Obs))/length(counsyl$EuroExAC))*1)-6.5
CC_Obs<-counsyl$freqEur_Counsyl[which(counsyl$freqEur_Counsyl != 0)]
CC_0<-(((length(counsyl$freqEur_Counsyl)-length(CC_Obs))/length(counsyl$freqEur_Counsyl))*1)-6.5

#pdf("FigS1_compare_ExAC_Counsyl.pdf",5,6.7)
setEPS()
postscript("FigS1_compare_ExAC_Counsyl.eps")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-6.5,-2),axes=FALSE,ann=FALSE)
vioplot(log10(CX_Obs), log10(CC_Obs),add=TRUE,col=c("white"),drawRect = FALSE)
axis(side=1,at=1:2,labels=c("ExAC","Counsyl"))
axis(side=2,at=c(-6.25,-5:0),labels=c(0,-5:0))
axis.break(axis = 2,breakpos = -5.5)
title(ylab="Log10(Frequency)")
rect(0.6,-6.5,1.4,CX_0)
rect(1.6,-6.5,2.4,CC_0)
x<-format(round((length(counsyl$EuroExAC)-length(CX_Obs))/length(counsyl$EuroExAC), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-6.39,paste(x,"%",sep=""),cex=1)
x<-format(round((length(counsyl$freqEur_Counsyl)-length(CC_Obs))/length(counsyl$freqEur_Counsyl), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-6.4,paste(x,"%",sep=""),cex=1)
points(x=1,y=log10(mean(counsyl$EuroExAC)), pch="-", cex=3,col="red")
points(x=2,y=log10(mean(counsyl$freqEur_Counsyl)), pch="-", cex=3,col="red")
dev.off()


#pdf("FigS1_compare_ExAC_Counsyl_HIST.pdf",5.7,8.7)
#par(mfrow=c(3,1))
#hist(ExAC$eur_freq,breaks = 14,xlim=c(0,0.002),freq=FALSE,main = "ExAC",xlab = "Empirical allele frequency",ylim = c(0,9000))
#hist(counsyl$EuroExAC,breaks = 14,xlim=c(0,0.002),freq=FALSE,main="ExAC (subset)",xlab = "Empirical allele frequency",ylim = c(0,9000))
#hist(counsyl$freqEur_Counsyl,breaks = 14,xlim=c(0,0.002),freq=FALSE,main="Counsyl",xlab = "Empirical allele frequency",ylim = c(0,9000))
#dev.off()


############################################################################################
################################        Results 2:          ################################
################################    Comparing expectations  ################################
############################################################################################

#############
### Fig 1 ###
#############

# Fig 2 
#pdf("Fig1_compare_expectations.pdf",7.5,7)
setEPS()
postscript("Fig1_compare_expectations.eps")
histSIM<-hist(log10(sims_AverULOG$V1),plot = FALSE)
plot(histSIM$counts,log="y", type='h', lwd=20, lend=1,col="grey",axes=FALSE,xlab = "Log10(Allele frequency)",ylab = "Density")
axis(side=1,at=c(1,6,11,16,21,26,31),labels=c("0","-6","-5","-4","-3","-2","-1"))
axis.break(axis = 1,breakpos = 3) #depends on "library(plotrix)"
abline(v=10.53991814,col="red",lwd = 2) #SIM
abline(v=9.143930767,col="green",lwd = 2) #FIN
abline(v=16.19795289,col="blue",lwd = 2) #INF
axis(side=2,at=c(10,100,1000,10000,100000,1000000),labels=c("1e+01","1e+02","1e+03","1e+04","1e+05","1e+06"))     
dev.off()

############
## Fig S2 ##
############

#pdf("FigS2_compare_FIN_SIM.pdf",6,8.5)
setEPS()
postscript("FigS2_compare_FIN_SIM.eps",width=6,height=8.5)
par(mfrow=c(2,1))
# Panel A - Expectation under the FIN model as a function of Ne.
SmallN<-7310 #This is the smallest N used for SIM
LargeN<-512000 #This is the largest N used for SIM (corresponds to current European pop size)
plot(0:1,0:1,type="n",xlim=c(0,LargeN),ylim=c(0,FIN(LargeN)),axes=FALSE,ann=FALSE)
plot(FIN,SmallN,LargeN,add=TRUE)
axis(side=1)
axis(side=2)
title(ylab="Deleterious allele frequency",xlab="Population size")
x<-mean(sims_AverU$V1)
size=function(q){((q/aver_u)^2)/(2*pi)}
abline(v=size(x),col="red")
title(outer=FALSE,adj=0,main="A",cex=2,col="black",font=2,line=1)
# Panel B - Simulation with constant size                                                               

SIMCONST_Obs<-constSize$X0[which(constSize$X0 != 0)]
SIMCONST_0<-(((length(constSize$X0)-length(SIMCONST_Obs))/length(constSize$X0))*1)-8

SIM_Obs<-sims_AverU$V1[which(sims_AverU$V1 != 0)]
SIM_0<-(((length(sims_AverU$V1)-length(SIM_Obs))/length(sims_AverU$V1))*1)-8

plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(SIM_Obs),log10(SIMCONST_Obs),col=c("grey"),add=TRUE,drawRect = FALSE)
axis(side=1,at=1:2,labels=c("SIM","FIN"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.7)
title(ylab="Log10(Frequency)")
points(x=1,y=log10(mean(sims_AverU$V1)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(constSize$X0)), pch="-", cex=3,col="red")
rect(0.6,-8,1.4,SIM_0,col="grey")
x<-format(round((length(sims_AverU$V1)-length(SIM_Obs))/length(sims_AverU$V1), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.6,paste(x,"%",sep=""),cex=1)
rect(1.6,-8,2.4,SIMCONST_0,col="grey")
x<-format(round((length(constSize$X0)-length(SIMCONST_Obs))/length(constSize$X0), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.56,paste(x,"%",sep=""),cex=1)
title(outer=FALSE,adj=0,main="B",cex=2,col="black",font=2,line=1)

dev.off()

#######################
## Other comparisons ##
#######################

#Comparing q across different models:
mean(sims_AverU$V1)/FIN20K #SIM mean is 1.90 larger than FIN20K
INF/mean(sims_AverU$V1) #SIM mean is 13.5 larger than INF
INF/FIN20K #FIN mean is 25.75 larger than INF

#Compare Tennessen demography to constant pop size, using predicted Ne from coupling SIM q to FIN equation
x<-mean(constSize$X0)
y<-mean(sims_AverU$V1)
z<-size(x)
k<-size(y)
x/y
z/k #constant size yields slightl larger size and frequency

#What is the population size in FIN that corresponds to the frequencies observed in SIM? 
size(x) # N = 72347.74 # This is going to be the number of DIPLOID individuals

#Mass at frequency = 0 is larger for Constant size (~98%) than for SIM (~87%)
length(which(sims_AverU$V1 == 0))/length(sims_AverU$V1)
length(which(constSize$X0 == 0))/length(constSize$X0)

# Distributions differ:
ks.test(sims_AverU$V1,constSize$X0) # p-value < 2.2e-16 # Warning: p-value will be approximate in the presence of ties

##########################################################################################
######################             Results 3           ###################################
######################      Comparing ExAC to SIM      ###################################
##########################################################################################

#############
### Fig 2 ###
#############

obsSim_CpGti <- simsPOI$CpGti[which(simsPOI$CpGti != 0)]
obsSim_CpGtv <- simsPOI$CpGtv[which(simsPOI$CpGtv != 0)]
obsSim_nonCpGti <- simsPOI$nonCpGti[which(simsPOI$nonCpGti != 0)]
obsSim_nonCpGtv <- simsPOI$nonCpGtv[which(simsPOI$nonCpGtv != 0)]
ZeroCpGti<-(((2000000-length(obsSim_CpGti))/2000000)*1)-8
ZeroCpGtv<-(((2000000-length(obsSim_CpGtv))/2000000)*1)-8
ZerononCpGti<-(((2000000-length(obsSim_nonCpGti))/2000000)*1)-8
ZerononCpGtv<-(((2000000-length(obsSim_nonCpGtv))/2000000)*1)-8
obsExAC_CpGti <- XCpGti$eur_freq[which(XCpGti$eur_freq != 0)]
obsExAC_CpGtv <- XCpGtv$eur_freq[which(XCpGtv$eur_freq != 0)]
obsExAC_nonCpGti <- XnonCpGti$eur_freq[which(XnonCpGti$eur_freq != 0)]
obsExAC_nonCpGtv <- XnonCpGtv$eur_freq[which(XnonCpGtv$eur_freq != 0)]
ZeroExACCpGti<-(((length(XCpGti$eur_freq)-length(obsExAC_CpGti))/length(XCpGti$eur_freq))*1)-8
ZeroExACCpGtv<-(((length(XCpGtv$eur_freq)-length(obsExAC_CpGtv))/length(XCpGtv$eur_freq))*1)-8
ZeroExACnonCpGti<-(((length(XnonCpGti$eur_freq)-length(obsExAC_nonCpGti))/length(XnonCpGti$eur_freq))*1)-8
ZeroExACnonCpGtv<-(((length(XnonCpGtv$eur_freq)-length(obsExAC_nonCpGtv))/length(XnonCpGtv$eur_freq))*1)-8

#pdf("Fig2_ExAC_SIM.pdf",8.7,8)
setEPS()
postscript("Fig2_ExAC_SIM.eps",width=7.5,height=7)
par(mfrow=c(2,2))
#CpGti
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(obsSim_CpGti),log10(obsExAC_CpGti),col=c("white","grey","white"),add=TRUE,drawRect = FALSE)
title("CpGti (n=101)\np-value=0.18")
axis(side=1,at=1:2,labels=c("SIM","ExAC"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#mtext("Proportion of non-segregating sites", side=4, line=3, cex=0.6,las=3)
#Means:
points(x=1,y=log10(mean(simsPOI$CpGti)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(XCpGti$eur_freq)), pch="-", cex=3,col="red")
#Mass at zero
rect(0.6,-8,1.4,ZeroCpGti,col="grey")
rect(1.6,-8,2.4,ZeroExACCpGti)
x<-format(round(((2000000-(length(obsSim_CpGti)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.6,paste(x,"%",sep=""),cex=1)
x<-format(round((length(XCpGti$eur_freq)-length(obsExAC_CpGti))/length(XCpGti$eur_freq), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.34,paste(x,"%",sep=""),cex=1)

#CpGtv
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(obsSim_CpGtv),log10(obsExAC_CpGtv),col=c("white","grey","white"),add=TRUE,drawRect = FALSE)
title("CpGtv (n=13)\np-value=0.04")
axis(side=1,at=1:2,labels=c("SIM","ExAC"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#Mass at zero
rect(0.6,-8,1.4,ZeroCpGtv,col="grey")
rect(1.6,-8,2.4,ZeroExACCpGtv)
x<-format(round(((2000000-(length(obsSim_CpGtv)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.5,paste(x,"%",sep=""),cex=1)
x<-format(round((length(XCpGtv$eur_freq)-length(obsExAC_CpGtv))/length(XCpGtv$eur_freq), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.75,paste(x,"%",sep=""),cex=1)
#Means:
points(x=1,y=log10(mean(simsPOI$CpGtv)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(XCpGtv$eur_freq)), pch="-", cex=3,col="red")

#nonCpGti
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(obsSim_nonCpGti),log10(obsExAC_nonCpGti),col=c("white","grey","white"),add=TRUE,drawRect = FALSE)
title("nonCpGti (n=155)\np-value=1e-4")
axis(side=1,at=1:2,labels=c("SIM","ExAC"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#Mass at zero
rect(0.6,-8,1.4,ZerononCpGti,col="grey")
rect(1.6,-8,2.4,ZeroExACnonCpGti)
x<-format(round(((2000000-(length(obsSim_nonCpGti)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.56,paste(x,"%",sep=""),cex=1)
x<-format(round((length(XnonCpGti$eur_freq)-length(obsExAC_nonCpGti))/length(XnonCpGti$eur_freq), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.7,paste(x,"%",sep=""),cex=1)
#Means:
points(x=1,y=log10(mean(simsPOI$nonCpGti)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(XnonCpGti$eur_freq)), pch="-", cex=3,col="red")

#nonCpGtv
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(obsSim_nonCpGtv),log10(obsExAC_nonCpGtv),col=c("white","grey","white"),add=TRUE,drawRect = FALSE)
title("nonCpGtv (n=104)\np-value=6e-5")
axis(side=1,at=1:2,labels=c("SIM","ExAC"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#Mass at zero
rect(0.6,-8,1.4,ZerononCpGtv,col="grey")
rect(1.6,-8,2.4,ZeroExACnonCpGtv)
x<-format(round(((2000000-(length(obsSim_nonCpGtv)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.5,paste(x,"%",sep=""),cex=1)
x<-format(round((length(XnonCpGtv$eur_freq)-length(obsExAC_nonCpGtv))/length(XnonCpGtv$eur_freq), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.7,paste(x,"%",sep=""),cex=1)
#Means:
points(x=1,y=log10(mean(simsPOI$nonCpGtv)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(XnonCpGtv$eur_freq)), pch="-", cex=3,col="red")

dev.off()

##########################################################################################
#Significance of deviation

#How many obs in the empirical data for each mut
obsA<-nrow(XCpGti)
obsB<-nrow(XCpGtv)
obsC<-nrow(XnonCpGti)
obsD<-nrow(XnonCpGtv)

#how many repetitions
#x=100000
x=10

#Sampling from simulations
sig_CpGti <- c()
for (i in 1:x){
  simSampled<-simsPOI[sample(1:nrow(simsPOI), obsA, replace=FALSE),]
  sig_CpGti <- append(sig_CpGti,mean(simSampled$CpGti))
  }
sig_CpGtv <- c()
for (i in 1:x){
  simSampled<-simsPOI[sample(1:nrow(simsPOI), obsB, replace=FALSE),]
  sig_CpGtv <- append(sig_CpGtv,mean(simSampled$CpGtv))
}
sig_nonCpGti <- c()
for (i in 1:x){
  simSampled<-simsPOI[sample(1:nrow(simsPOI), obsC, replace=FALSE),]
  sig_nonCpGti <- append(sig_nonCpGti,mean(simSampled$nonCpGti))
}
sig_nonCpGtv <- c()
for (i in 1:x){
  simSampled<-simsPOI[sample(1:nrow(simsPOI), obsD, replace=FALSE),]
  sig_nonCpGtv <- append(sig_nonCpGtv,mean(simSampled$nonCpGtv))
}

#Significance (note that it may vary from the main text because it is subject to stochascity of the resampling scheme)
perc.rank <- function(x, xo)  (1-(length(x[x <= xo])/length(x)))*2
perc.rank(x = sig_CpGti, xo = mean(XCpGti$eur_freq) )   
perc.rank(x = sig_CpGtv, xo = mean(XCpGtv$eur_freq) )
perc.rank(x = sig_nonCpGti, xo = mean(XnonCpGti$eur_freq) )
perc.rank(x = sig_nonCpGtv, xo = mean(XnonCpGtv$eur_freq) )

##########################################################################################
#Compare empirical to expected
mean(XCpGti$eur_freq)/mean(simsPOI$CpGti)
mean(XCpGtv$eur_freq)/mean(simsPOI$CpGtv)
mean(XnonCpGti$eur_freq)/mean(simsPOI$nonCpGti)
mean(XnonCpGtv$eur_freq)/mean(simsPOI$nonCpGtv)


##########################################################################################
######################              Results 4            #################################
######################        Gene level analysis        #################################
##########################################################################################

#Editing data
geneExpEur<-colMeans(eurbygenePOI) # get mean allele frequency for each gene
aggGeneList<-aggregate(. ~ gene, data=ExAC_noGCI,FUN=sum)
geneTable<-rbind(geneExpEur,aggGeneList$eur_freq) #make matrix out of the 4 variables
geneTable<-t(geneTable) #transpose
colnames(geneTable)<-c("ExpEur","ObsEur")
geneTable<-as.data.frame(geneTable) #Convert matrix to data frame

#Calculating the significance (27 genes): # notice that I'm not computing it for 5 genes that have only 1 mut described (See Methods)
perc.rank <- function(x, xo)  (1-(length(x[x <= xo])/length(x)))*2
perc.rank(x = eurbygenePOILOG$ASPA, xo = geneTable[1,2] )
perc.rank(x = eurbygenePOILOG$ASS1, xo = geneTable[2,2] )
perc.rank(x = eurbygenePOILOG$CFTR, xo = geneTable[4,2] )
perc.rank(x = eurbygenePOILOG$CLN5, xo = geneTable[5,2] )
perc.rank(x = eurbygenePOILOG$DHCR7, xo = geneTable[6,2] )
perc.rank(x = eurbygenePOILOG$ERCC8, xo = geneTable[7,2] )
perc.rank(x = eurbygenePOILOG$FAH, xo = geneTable[8,2] )
perc.rank(x = eurbygenePOILOG$GAA, xo = geneTable[10,2] )
perc.rank(x = eurbygenePOILOG$GALC, xo = geneTable[11,2] )
perc.rank(x = eurbygenePOILOG$GAN, xo = geneTable[12,2] )
perc.rank(x = eurbygenePOILOG$GBE1, xo = geneTable[13,2] )
perc.rank(x = eurbygenePOILOG$HEXA, xo = geneTable[14,2] )
perc.rank(x = eurbygenePOILOG$HSD17B4, xo = geneTable[15,2] )
perc.rank(x = eurbygenePOILOG$IDUA, xo = geneTable[16,2] )
perc.rank(x = eurbygenePOILOG$LAMB3, xo = geneTable[18,2] )
perc.rank(x = eurbygenePOILOG$NPC1, xo = geneTable[19,2] )
perc.rank(x = eurbygenePOILOG$PEX7, xo = geneTable[20,2] )
perc.rank(x = eurbygenePOILOG$POLG, xo = geneTable[22,2] )
perc.rank(x = eurbygenePOILOG$POMGNT1, xo = geneTable[23,2] )
perc.rank(x = eurbygenePOILOG$PPT1, xo = geneTable[24,2] )
perc.rank(x = eurbygenePOILOG$PRF1, xo = geneTable[25,2] )
perc.rank(x = eurbygenePOILOG$SLC22A5, xo = geneTable[26,2] )
perc.rank(x = eurbygenePOILOG$SMARCAL1, xo = geneTable[27,2] )
perc.rank(x = eurbygenePOILOG$SMPD1, xo = geneTable[28,2] )
perc.rank(x = eurbygenePOILOG$STAR, xo = geneTable[29,2] )
perc.rank(x = eurbygenePOILOG$TK2, xo = geneTable[31,2] )
perc.rank(x = eurbygenePOILOG$TPP1, xo = geneTable[32,2] )

howManyObs <- function(x)  x[which(x != 0)]
ASPAObs<-howManyObs(eurbygenePOI$ASPA)
ASS1Obs<-howManyObs(eurbygenePOI$ASS1)
CFTRObs<-howManyObs(eurbygenePOI$CFTR)
CLN5Obs<-howManyObs(eurbygenePOI$CLN5)
DHCR7Obs<-howManyObs(eurbygenePOI$DHCR7)
ERCC8Obs<-howManyObs(eurbygenePOI$ERCC8)
FAHObs<-howManyObs(eurbygenePOI$FAH)
GAAObs<-howManyObs(eurbygenePOI$GAA)
GALCObs<-howManyObs(eurbygenePOI$GALC)
GANObs<-howManyObs(eurbygenePOI$GAN)
GBE1Obs<-howManyObs(eurbygenePOI$GBE1)
HEXAObs<-howManyObs(eurbygenePOI$HEXA)
HSD17B4Obs<-howManyObs(eurbygenePOI$HSD17B4)
IDUAObs<-howManyObs(eurbygenePOI$IDUA)
LAMB3Obs<-howManyObs(eurbygenePOI$LAMB3)
NPC1Obs<-howManyObs(eurbygenePOI$NPC1)
PEX7Obs<-howManyObs(eurbygenePOI$PEX7)
POLGObs<-howManyObs(eurbygenePOI$POLG)
POMGNT1Obs<-howManyObs(eurbygenePOI$POMGNT1)
PPT1Obs<-howManyObs(eurbygenePOI$PPT1)
PRF1Obs<-howManyObs(eurbygenePOI$PRF1)
SLC22A5Obs<-howManyObs(eurbygenePOI$SLC22A5)
SMARCAL1Obs<-howManyObs(eurbygenePOI$SMARCAL1)
SMPD1Obs<-howManyObs(eurbygenePOI$SMPD1)
STARObs<-howManyObs(eurbygenePOI$STAR)
TK2Obs<-howManyObs(eurbygenePOI$TK2)
TPP1Obs<-howManyObs(eurbygenePOI$TPP1)

howManyZeros <- function(x, xo) (((length(x)-length(xo))/length(x))*1)-6
ASPAzeros<-howManyZeros(eurbygenePOI$ASPA,ASPAObs)
ASS1zeros<-howManyZeros(eurbygenePOI$ASS1,ASS1Obs)
CFTRzeros<-howManyZeros(eurbygenePOI$CFTR,CFTRObs)
CLN5zeros<-howManyZeros(eurbygenePOI$CLN5,CLN5Obs)
DHCR7zeros<-howManyZeros(eurbygenePOI$DHCR7,DHCR7Obs)
ERCC8zeros<-howManyZeros(eurbygenePOI$ERCC8,ERCC8Obs)
FAHzeros<-howManyZeros(eurbygenePOI$FAH,FAHObs)
GAAzeros<-howManyZeros(eurbygenePOI$GAA,GAAObs)
GALCzeros<-howManyZeros(eurbygenePOI$GALC,GALCObs)
GANzeros<-howManyZeros(eurbygenePOI$GAN,GANObs)
GBE1zeros<-howManyZeros(eurbygenePOI$GBE1,GBE1Obs)
HEXAzeros<-howManyZeros(eurbygenePOI$HEXA,HEXAObs)
HSD17B4zeros<-howManyZeros(eurbygenePOI$HSD17B4,HSD17B4Obs)
IDUAzeros<-howManyZeros(eurbygenePOI$IDUA,IDUAObs)
LAMB3zeros<-howManyZeros(eurbygenePOI$LAMB3,LAMB3Obs)
NPC1zeros<-howManyZeros(eurbygenePOI$NPC1,NPC1Obs)
PEX7zeros<-howManyZeros(eurbygenePOI$PEX7,PEX7Obs)
POLGzeros<-howManyZeros(eurbygenePOI$POLG,POLGObs)
POMGNT1zeros<-howManyZeros(eurbygenePOI$POMGNT1,POMGNT1Obs)
PPT1zeros<-howManyZeros(eurbygenePOI$PPT1,PPT1Obs)
PRF1zeros<-howManyZeros(eurbygenePOI$PRF1,PRF1Obs)
SLC22A5zeros<-howManyZeros(eurbygenePOI$SLC22A5,SLC22A5Obs)
SMARCAL1zeros<-howManyZeros(eurbygenePOI$SMARCAL1,SMARCAL1Obs)
SMPD1zeros<-howManyZeros(eurbygenePOI$SMPD1,SMPD1Obs)
STARzeros<-howManyZeros(eurbygenePOI$STAR,STARObs)
TK2zeros<-howManyZeros(eurbygenePOI$TK2,TK2Obs)
TPP1zeros<-howManyZeros(eurbygenePOI$TPP1,TPP1Obs)

#############
### Fig 3 ###
#############

#Genes are ordered by significance of the deviation (ie p-value above - See Table S4)

#pdf("Fig3_geneLevel.pdf",9,5)
setEPS()
postscript("Fig3_geneLevel.eps",width=7.5,height=4.17)
par(mar=c(5.1,4.1,4.1,5.1))
plot(0:1,0:1,type="n",xlim=c(1,27),ylim=c(-6,-2),axes=FALSE,ann=FALSE)
vioplot(log10(CFTRObs),log10(IDUAObs),log10(SMPD1Obs),log10(POLGObs),log10(FAHObs),log10(LAMB3Obs),log10(ERCC8Obs),log10(ASPAObs),log10(HSD17B4Obs),log10(PPT1Obs),log10(DHCR7Obs),log10(PEX7Obs),log10(PRF1Obs),log10(ASS1Obs),log10(SMARCAL1Obs),log10(NPC1Obs),log10(GBE1Obs),log10(GALCObs),log10(TK2Obs),log10(STARObs),log10(CLN5Obs),log10(TPP1Obs),log10(GAAObs),log10(SLC22A5Obs),log10(POMGNT1Obs),log10(GANObs),log10(HEXAObs),col=c("grey"),drawRect=FALSE,add=TRUE)
#vioplot(log10(eurbygenePOILOG$CFTR),log10(eurbygenePOILOG$IDUA),log10(eurbygenePOILOG$SMPD1),log10(eurbygenePOILOG$POLG),log10(eurbygenePOILOG$FAH),log10(eurbygenePOILOG$LAMB3),log10(eurbygenePOILOG$ERCC8),log10(eurbygenePOILOG$ASPA),log10(eurbygenePOILOG$HSD17B4),log10(eurbygenePOILOG$PPT1),log10(eurbygenePOILOG$DHCR7),log10(eurbygenePOILOG$PEX7),log10(eurbygenePOILOG$PRF1),log10(eurbygenePOILOG$ASS1),log10(eurbygenePOILOG$SMARCAL1),log10(eurbygenePOILOG$NPC1),log10(eurbygenePOILOG$GBE1),log10(eurbygenePOILOG$GALC),log10(eurbygenePOILOG$TK2),log10(eurbygenePOILOG$STAR),log10(eurbygenePOILOG$CLN5),log10(eurbygenePOILOG$TPP1),log10(eurbygenePOILOG$GAA),log10(eurbygenePOILOG$SLC22A5),log10(eurbygenePOILOG$POMGNT1),log10(eurbygenePOILOG$GAN),log10(eurbygenePOILOG$HEXA),col=c("grey"),add=TRUE)
axis(1, at=seq(1, 27, by=1), labels = FALSE)
text(seq(1, 27, by=1), par("usr")[3] - 0.2, labels = c("CFTR","IDUA","SMPD1","POLG","FAH","LAMB3","ERCC8","ASPA","HSD17B4","PPT1","DHCR7","PEX7","PRF1","ASS1","SMARCAL1","NPC1","GBE1","GALC","TK2","STAR","CLN5","TPP1","GAA","SLC22A5","POMGNT1","GAN","HEXA"), srt = 45, pos = 1, xpd = TRUE,cex = 0.7)
axis(side=2,at=c(-5.75,-5:-2),labels=c(0,-5:-2))
axis.break(axis = 2,breakpos = -5.25)
title(ylab="Log10(Frequency)")

#Means:
points(x=1,y=log10(geneTable[4,2]),pch="-",cex=2,col="darkmagenta")
points(x=2,y=log10(geneTable[16,2]),pch="-",cex=2,col="darkmagenta")
points(x=3,y=log10(geneTable[28,2]),pch="-",cex=2,col="darkmagenta")
points(x=4,y=log10(geneTable[22,2]),pch="-",cex=2,col="darkmagenta")
points(x=5,y=log10(geneTable[8,2]),pch="-",cex=2,col="darkmagenta")
points(x=6,y=log10(geneTable[18,2]),pch="-",cex=2,col="darkmagenta")
points(x=7,y=log10(geneTable[7,2]),pch="-",cex=2,col="darkmagenta")
points(x=8,y=log10(geneTable[1,2]),pch="-",cex=2,col="darkmagenta")
points(x=9,y=log10(geneTable[15,2]),pch="-",cex=2,col="darkmagenta")
points(x=10,y=log10(geneTable[24,2]),pch="-",cex=2,col="darkmagenta")
points(x=11,y=log10(geneTable[6,2]),pch="-",cex=2,col="darkmagenta")
points(x=12,y=log10(geneTable[20,2]),pch="-",cex=2,col="darkmagenta")
points(x=13,y=log10(geneTable[25,2]),pch="-",cex=2,col="darkmagenta")
points(x=14,y=log10(geneTable[2,2]),pch="-",cex=2,col="darkmagenta")
points(x=15,y=log10(geneTable[27,2]),pch="-",cex=2,col="darkmagenta")
points(x=16,y=log10(geneTable[19,2]),pch="-",cex=2,col="darkmagenta")
points(x=17,y=log10(geneTable[13,2]),pch="-",cex=2,col="darkmagenta")
points(x=18,y=log10(geneTable[11,2]),pch="-",cex=2,col="darkmagenta")
points(x=19,y=log10(geneTable[31,2]),pch="-",cex=2,col="darkmagenta")
points(x=20,y=log10(geneTable[29,2]),pch="-",cex=2,col="darkmagenta")
points(x=21,y=log10(geneTable[5,2]),pch="-",cex=2,col="darkmagenta")
points(x=22,y=log10(geneTable[32,2]),pch="-",cex=2,col="darkmagenta")
points(x=23,y=log10(geneTable[10,2]),pch="-",cex=2,col="darkmagenta")
points(x=24,y=log10(geneTable[26,2]),pch="-",cex=2,col="darkmagenta")
points(x=25,y=log10(geneTable[23,2]),pch="-",cex=2,col="darkmagenta")
points(x=26,y=log10(geneTable[12,2]),pch="-",cex=2,col="darkmagenta")
points(x=27,y=log10(geneTable[14,2]),pch="-",cex=2,col="darkmagenta")


rect(0.6,-6,1.4,CFTRzeros,col="grey")
rect(1.6,-6,2.4,IDUAzeros,col="grey")
rect(2.6,-6,3.4,SMPD1zeros,col="grey")
rect(3.6,-6,4.4,POLGzeros,col="grey")
rect(4.6,-6,5.4,FAHzeros,col="grey")
rect(5.6,-6,6.4,LAMB3zeros,col="grey")
rect(6.6,-6,7.4,ERCC8zeros,col="grey")
rect(7.6,-6,8.4,ASPAzeros,col="grey")
rect(8.6,-6,9.4,HSD17B4zeros,col="grey")
rect(9.6,-6,10.4,PPT1zeros,col="grey")
rect(10.6,-6,11.4,DHCR7zeros,col="grey")
rect(11.6,-6,12.4,PEX7zeros,col="grey")
rect(12.6,-6,13.4,PRF1zeros,col="grey")
rect(13.6,-6,14.4,ASS1zeros,col="grey")
rect(14.6,-6,15.4,SMARCAL1zeros,col="grey")
rect(15.6,-6,16.4,NPC1zeros,col="grey")
rect(16.6,-6,17.4,GBE1zeros,col="grey")
rect(17.6,-6,18.4,GALCzeros,col="grey")
rect(18.6,-6,19.4,TK2zeros,col="grey")
rect(19.6,-6,20.4,STARzeros,col="grey")
rect(20.6,-6,21.4,CLN5zeros,col="grey")
rect(21.6,-6,22.4,TPP1zeros,col="grey")
rect(22.6,-6,23.4,GAAzeros,col="grey")
rect(23.6,-6,24.4,SLC22A5zeros,col="grey")
rect(24.6,-6,25.4,POMGNT1zeros,col="grey")
rect(25.6,-6,26.4,GANzeros,col="grey")
rect(26.6,-6,27.4,HEXAzeros,col="grey")
axis(side=4,at=c(-6,-5.75,-5.5,-5.25,-5),labels=c("0","25%","50%","75%","100%"),las=1)
#mtext("Prop. not segr.", side=4, line=3, cex.lab=0.5,las=2)
dev.off()

#Is the combined p-value significant? Answer: Yes! P-value: 4.424856e-15 
pvals<-c(0.002,0.002,0.01,0.018,0.02,0.02,0.03,0.032,0.044,0.048,0.052,0.054,0.132,0.144,0.144,0.194,0.218,0.3,0.318,0.454,0.456,0.464,0.476,0.496,0.62,0.75,0.85)
pchisq(-2*sum(log(pvals)),length(pvals),lower.tail=FALSE)

##########################################################################################
######################             Results 5            ##################################
######################      Selective effect in hets    ##################################
##########################################################################################

############
## Fig S3 ##
############

howManyObs <- function(x)  x[which(x != 0)]
recObs<-howManyObs(sims_AverU$V1)
effObs<-howManyObs(sims_EffectHet$V1)

howManyZeros <- function(x, xo) (((length(x)-length(xo))/length(x))*1)-8
recZero<-howManyZeros(sims_AverU$V1,recObs)
effZero<-howManyZeros(sims_EffectHet$V1,effObs)

#pdf("FigS3_effectInHets.pdf",5.3,7)
setEPS()
postscript("FigS3_effectInHets.eps",width=5.3,height=7)
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(recObs),log10(effObs),col=c("grey"),drawRect=FALSE,add=TRUE)
axis(side=1,at=1:2,labels=c("h = 0","h = 1%"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#Means:
points(x=1,y=log10(mean(sims_AverU$V1)), pch="-",cex=2,col="red")
points(x=2,y=log10(mean(sims_EffectHet$V1)), pch="-", cex=2,col="red")
#Mass at zero
rect(0.6,-8,1.4,recZero,col="grey")
rect(1.6,-8,2.4,effZero,col="grey")
x<-format(round(((2000000-(length(recObs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.6,paste(x,"%",sep=""),cex=1)
x<-format(round(((2000000-(length(effObs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.55,paste(x,"%",sep=""),cex=1)

dev.off()

#How does the effect in het impact the mean allele freq
mean(sims_EffectHet$V1)/mean(sims_AverU$V1) #H=0.01 decreases mean freq of deleterious alleles to 42%

#Distributions are different according to WMW test
wilcox.test(sims_EffectHet$V1,sims_AverU$V1) # p-value < 2.2e-16
ks.test(sims_EffectHet$V1,sims_AverU$V1) # p-value < 2.2e-16


##############################################################################################
######################             Results 6                ##################################
######################      Consider more growth for SIM    ##################################
##############################################################################################

############
## Fig S4 ##
############

howManyObs <- function(x)  x[which(x != 0)]
Tobs<-howManyObs(sims_AverUPOI$V1)
Aobs<-howManyObs(sim1MPOI$V1)
Bobs<-howManyObs(sim2MPOI$V1)
Cobs<-howManyObs(sim5MPOI$V1)
eobs<-howManyObs(ExAC$eur_freq)

howManyZeros <- function(x, xo) (((length(x)-length(xo))/length(x))*1)-7
Tzero<-howManyZeros(sims_AverUPOI$V1,Tobs)
Azero<-howManyZeros(sim1MPOI$V1,Aobs)
Bzero<-howManyZeros(sim2MPOI$V1,Bobs)
Czero<-howManyZeros(sim5MPOI$V1,Cobs)
ezero<-howManyZeros(ExAC$eur_freq,eobs)

#pdf("FigS4_SIM_more_growth.pdf",6,7)
setEPS()
postscript("FigS4_SIM_more_growth.eps",width=6,height=7)
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=c(-7,-2),axes=FALSE,ann=FALSE)
vioplot(log10(Tobs),log10(Aobs),log10(Bobs),log10(Cobs),log10(eobs),drawRect=FALSE,col=c("grey","grey","grey","grey","grey","white"),add=TRUE)
axis(side=1,at=1:5,labels=c("A","B","C","D","ExAC"))
axis(side=2,at=c(-6.65,-5:-2),labels=c(0,-5:-2))
axis.break(axis = 2,breakpos = -5.7)
#title(ylab="Log10(Frequency)",xlab="Current population size")
title(ylab="Log10(Frequency)",)
#Mass at zero
rect(0.6,-7,1.4,Tzero,col="grey")
rect(1.6,-7,2.4,Azero,col="grey")
rect(2.6,-7,3.4,Bzero,col="grey")
rect(3.6,-7,4.4,Czero,col="grey")
rect(4.6,-7,5.4,ezero,col="white")
#Means:
points(x=1,y=log10(mean(sims_AverUPOI$V1)), pch="-",cex=2,col="red")
points(x=2,y=log10(mean(sim1MPOI$V1)), pch="-", cex=2,col="red")
points(x=3,y=log10(mean(sim2MPOI$V1)), pch="-", cex=2,col="red")
points(x=4,y=log10(mean(sim5MPOI$V1)), pch="-", cex=2,col="red")
points(x=5,y=log10(mean(ExAC$eur_freq)), pch="-", cex=2,col="red")

x<-format(round(((2000000-(length(Tobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-6.65,paste(x,"%",sep=""),cex=1)
x<-format(round(((2000000-(length(Aobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-6.65,paste(x,"%",sep=""),cex=1)
x<-format(round(((2000000-(length(Bobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(3,-6.65,paste(x,"%",sep=""),cex=1)
x<-format(round(((2000000-(length(Cobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(4,-6.65,paste(x,"%",sep=""),cex=1)
x<-format(round(((nrow(ExAC)-(length(eobs)))/nrow(ExAC)), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(5,-6.73,paste(x,"%",sep=""),cex=1)

dev.off()

wilcox.test(sims_AverU$V1,sim1M$V1)
wilcox.test(sims_AverU$V1,sim2M$V1)
wilcox.test(sims_AverU$V1,sim5M$V1)

mean(sims_AverU$V1)/mean(sim1M$V1) # 0.9800684
mean(sims_AverU$V1)/mean(sim2M$V1) # 0.9861226
mean(sims_AverU$V1)/mean(sim5M$V1) # 0.9709141

length(which(sims_AverU$V1 == 0))/length(sims_AverU$V1) # 0.868952
length(which(sim1M$V1 == 0))/length(sim1M$V1) # 0.789466 
length(which(sim2M$V1 == 0))/length(sim2M$V1) # 0.677886

wilcox.test(sim1M$V1,sims_AverU$V1)
ks.test(sim1M$V1,sims_AverU$V1)

wilcox.test(sim2M$V1,sims_AverU$V1)
ks.test(sim2M$V1,sims_AverU$V1)


##############################################################################################
######################             Results 7                   ###############################
######################      Reploting Fig 2 without 2 genes    ###############################
##############################################################################################

############
## Fig S5 ##
############

#Remove 32 genes from ExAC dataset
exac_wo2genes<-subset(ExAC_noGCI,gene!="CFTR" & gene!="DHCR7")
x<-subset(ExAC_noGCI,gene=="CFTR")
y<-subset(ExAC_noGCI,gene=="DHCR7")
nrow(x)
nrow(y)

XnonCpGtv <- subset(exac_wo2genes,type=="nonCpGtv")
XnonCpGti <- subset(exac_wo2genes,type=="nonCpGti")
XCpGtv <- subset(exac_wo2genes,type=="CpGtv")
XCpGti <- subset(exac_wo2genes,type=="CpGti")
obsSim_CpGti <- simsPOI$CpGti[which(simsPOI$CpGti != 0)]
obsSim_CpGtv <- simsPOI$CpGtv[which(simsPOI$CpGtv != 0)]
obsSim_nonCpGti <- simsPOI$nonCpGti[which(simsPOI$nonCpGti != 0)]
obsSim_nonCpGtv <- simsPOI$nonCpGtv[which(simsPOI$nonCpGtv != 0)]
ZeroCpGti<-(((2000000-length(obsSim_CpGti))/2000000)*1)-8
ZeroCpGtv<-(((2000000-length(obsSim_CpGtv))/2000000)*1)-8
ZerononCpGti<-(((2000000-length(obsSim_nonCpGti))/2000000)*1)-8
ZerononCpGtv<-(((2000000-length(obsSim_nonCpGtv))/2000000)*1)-8
obsExAC_CpGti <- XCpGti$eur_freq[which(XCpGti$eur_freq != 0)]
obsExAC_CpGtv <- XCpGtv$eur_freq[which(XCpGtv$eur_freq != 0)]
obsExAC_nonCpGti <- XnonCpGti$eur_freq[which(XnonCpGti$eur_freq != 0)]
obsExAC_nonCpGtv <- XnonCpGtv$eur_freq[which(XnonCpGtv$eur_freq != 0)]
ZeroExACCpGti<-(((length(XCpGti$eur_freq)-length(obsExAC_CpGti))/length(XCpGti$eur_freq))*1)-8
ZeroExACCpGtv<-(((length(XCpGtv$eur_freq)-length(obsExAC_CpGtv))/length(XCpGtv$eur_freq))*1)-8
ZeroExACnonCpGti<-(((length(XnonCpGti$eur_freq)-length(obsExAC_nonCpGti))/length(XnonCpGti$eur_freq))*1)-8
ZeroExACnonCpGtv<-(((length(XnonCpGtv$eur_freq)-length(obsExAC_nonCpGtv))/length(XnonCpGtv$eur_freq))*1)-8

#pdf("FigS5_wo_2genes.pdf",8.7,8)
setEPS()
postscript("FigS5_wo_2genes.eps",width=7.5,height=7)
par(mfrow=c(2,2))
#CpGti
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(obsSim_CpGti),log10(obsExAC_CpGti),col=c("white","grey","white"),add=TRUE,drawRect = FALSE)
title("CpGti (n=80)\np-value=0.30")
axis(side=1,at=1:2,labels=c("SIM","ExAC"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#mtext("Proportion of non-segregating sites", side=4, line=3, cex=0.6,las=3)
#Means:
points(x=1,y=log10(mean(simsPOI$CpGti)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(XCpGti$eur_freq)), pch="-", cex=3,col="red")
#Mass at zero
rect(0.6,-8,1.4,ZeroCpGti,col="grey")
rect(1.6,-8,2.4,ZeroExACCpGti)
x<-format(round(((2000000-(length(obsSim_CpGti)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.6,paste(x,"%",sep=""),cex=0.85)
x<-format(round((length(XCpGti$eur_freq)-length(obsExAC_CpGti))/length(XCpGti$eur_freq), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.3,paste(x,"%",sep=""),cex=0.85)

#CpGtv
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(obsSim_CpGtv),log10(obsExAC_CpGtv),col=c("white","grey","white"),add=TRUE,drawRect = FALSE)
title("CpGtv (n=9)\np-value=0.02")
axis(side=1,at=1:2,labels=c("SIM","ExAC"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#Mass at zero
rect(0.6,-8,1.4,ZeroCpGtv,col="grey")
rect(1.6,-8,2.4,ZeroExACCpGtv)
x<-format(round(((2000000-(length(obsSim_CpGtv)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.5,paste(x,"%",sep=""),cex=0.85)
x<-format(round((length(XCpGtv$eur_freq)-length(obsExAC_CpGtv))/length(XCpGtv$eur_freq), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.6,paste(x,"%",sep=""),cex=0.85)
#Means:
points(x=1,y=log10(mean(simsPOI$CpGtv)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(XCpGtv$eur_freq)), pch="-", cex=3,col="red")

#nonCpGti
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(obsSim_nonCpGti),log10(obsExAC_nonCpGti),col=c("white","grey","white"),add=TRUE,drawRect = FALSE)
title("nonCpGti (n=83)\np-value=0.01")
axis(side=1,at=1:2,labels=c("SIM","ExAC"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#Mass at zero
rect(0.6,-8,1.4,ZerononCpGti,col="grey")
rect(1.6,-8,2.4,ZeroExACnonCpGti)
x<-format(round(((2000000-(length(obsSim_nonCpGti)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.5,paste(x,"%",sep=""),cex=0.85)
x<-format(round((length(XnonCpGti$eur_freq)-length(obsExAC_nonCpGti))/length(XnonCpGti$eur_freq), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.6,paste(x,"%",sep=""),cex=0.85)
#Means:
points(x=1,y=log10(mean(simsPOI$nonCpGti)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(XnonCpGti$eur_freq)), pch="-", cex=3,col="red")

#nonCpGtv
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(-8,-1),axes=FALSE,ann=FALSE)
vioplot(log10(obsSim_nonCpGtv),log10(obsExAC_nonCpGtv),col=c("white","grey","white"),add=TRUE,drawRect = FALSE)
title("nonCpGtv (n=59)\np-value<0.01")
axis(side=1,at=1:2,labels=c("SIM","ExAC"))
axis(side=2,at=c(-7.5,-6:0),labels=c(0,-6:0))
axis.break(axis = 2,breakpos = -6.6)
title(ylab="Log10(Frequency)")
#Mass at zero
rect(0.6,-8,1.4,ZerononCpGtv,col="grey")
rect(1.6,-8,2.4,ZeroExACnonCpGtv)
x<-format(round(((2000000-(length(obsSim_nonCpGtv)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-7.5,paste(x,"%",sep=""),cex=0.85)
x<-format(round((length(XnonCpGtv$eur_freq)-length(obsExAC_nonCpGtv))/length(XnonCpGtv$eur_freq), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-7.7,paste(x,"%",sep=""),cex=0.85)
#Means:
points(x=1,y=log10(mean(simsPOI$nonCpGtv)), pch="-",cex=3,col="red")
points(x=2,y=log10(mean(XnonCpGtv$eur_freq)), pch="-", cex=3,col="red")

dev.off()

##########################################################################################
#Significance of deviation

#How many obs in the empirical data for each mut
obsA<-nrow(XCpGti)
obsB<-nrow(XCpGtv)
obsC<-nrow(XnonCpGti)
obsD<-nrow(XnonCpGtv)

#how many repetitions
#x=1000
x=10

#Sampling from simulations
sig_CpGti <- c()
for (i in 1:x){
  simSampled<-simsPOI[sample(1:nrow(simsPOI), obsA, replace=FALSE),]
  sig_CpGti <- append(sig_CpGti,mean(simSampled$CpGti))
}
sig_CpGtv <- c()
for (i in 1:x){
  simSampled<-simsPOI[sample(1:nrow(simsPOI), obsB, replace=FALSE),]
  sig_CpGtv <- append(sig_CpGtv,mean(simSampled$CpGtv))
}
sig_nonCpGti <- c()
for (i in 1:x){
  simSampled<-simsPOI[sample(1:nrow(simsPOI), obsC, replace=FALSE),]
  sig_nonCpGti <- append(sig_nonCpGti,mean(simSampled$nonCpGti))
}
sig_nonCpGtv <- c()
for (i in 1:x){
  simSampled<-simsPOI[sample(1:nrow(simsPOI), obsD, replace=FALSE),]
  sig_nonCpGtv <- append(sig_nonCpGtv,mean(simSampled$nonCpGtv))
}
#Significance
perc.rank <- function(x, xo)  (1-(length(x[x <= xo])/length(x)))*2
perc.rank(x = sig_CpGti, xo = mean(XCpGti$eur_freq) )   #0.28
perc.rank(x = sig_CpGtv, xo = mean(XCpGtv$eur_freq) )   #0.044
perc.rank(x = sig_nonCpGti, xo = mean(XnonCpGti$eur_freq) ) #0.008
perc.rank(x = sig_nonCpGtv, xo = mean(XnonCpGtv$eur_freq) ) # 0

##############################################################################################
##############################################################################################
######################      Comparing coverage between subsets    ############################
##############################################################################################
##############################################################################################

#pdf("FigS6_coverage.pdf",6,7)
setEPS()
postscript("FigS6_coverage.eps",width=6,height=7)
zero<-read.table("../FinalData/zeros_coverage.txt",header=FALSE)
nonzero<-read.table("../FinalData/nonZeros_coverage.txt",header=FALSE)
data<-data.frame(zero)
boxplot(zero,at=2,xlim=c(0,3))
boxplot(nonzero,at=1,add=TRUE)
axis(side = 1,at = c(1,2),labels = c("A","B"))
title(ylab="Coverage")
mean(zero$V1)
mean(nonzero$V1)
dev.off()
