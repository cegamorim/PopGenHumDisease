#############################################################################################################
### Amorim et al. (2017) The population genetics of human disease: the case of recessive, lethal mutations###
#############################################################################################################

#Author: 	Eduardo Amorim (with contributions from Ziyue Gao)
#Date: 		May, 2017
#Contact: 	guerraamorim@gmail.com

##################################################################################################
#load external packages needed to use color for violin plots
##################################################################################################
library(sm)
library(vioplot)
source("colorVioplot.r")
library(plotrix)

##################################################################################################
#################        Load data (all available as supplementary material)     #################
##################################################################################################

######################
### Empirical data ###
######################

counsyl<-read.table("counsyl_exac_matched.tab", header=TRUE)  					#Matched sites for Counsyl and ExAC                
ExAC<-read.table("ExAC_final385_filtered.tab.txt", header=TRUE) 				#ExAC - all sites             
ExAC_noGCI<-read.table("ExAC_final373_filtered_GCIfilter.tab", header=TRUE) 	#ExAC - sites not in CpG island

#Subset ExAC data by mutation type
XnonCpGtv <- subset(ExAC_noGCI,type=="nonCpGtv")
XnonCpGti <- subset(ExAC_noGCI,type=="nonCpGti")
XCpGtv <- subset(ExAC_noGCI,type=="CpGtv")
XCpGti <- subset(ExAC_noGCI,type=="CpGti")

######################
### Simulated data ###
######################

# Simulated allele frequencies for site-level analysis (as opposed to gene-level analysis) using average mutation rate of 1.5e-8 estimated for exons (Neale et al. 2012)
sims_AverUExons<-read.table("uFixed_exon_rescaled.txt", header=FALSE)

# Simulated allele frequencies with 1.5x larger mutation rate 
sims_LargerAverU<-read.table("sims_currAF_Eur_varU_1.5x_Larger_average_u.txt", header=FALSE)

# Simulated allele frequencies considering a small effect in hets 
sims_EffectHet<-read.table("eurDistrib_effectHets.txt", header=FALSE)

# Simulated allele frequencies using average mutation rate u from Kong et al. 2012 for each mutation type (CpGti, CpGtv, nonCpGti, and nonCpGtv)
sims_fourTypes<-read.table("4types_rescalled.tab", header=TRUE)

# Simulated allele frequencies using a constant population size for European populations
constSize<-read.table("eurDistrib_constantsize.txt",header=TRUE)

# Simulated allele frequencies with different demography, considering larger Ne in the present and more intense growth
sim1M<-read.table("1M_eurDistrib.txt",header=FALSE) # this is the same as sims_AverUExons, but here I used uVar
sim2M<-read.table("2M_eurDistrib.txt",header=FALSE)
sim5M<-read.table("5M_eurDistrib.txt",header=FALSE)
sim10M<-read.table("10M_eurDistrib.txt",header=FALSE)
simIG<-read.table("nophase1_eurDistrib.txt",header=FALSE)

# Simulations by gene (equivalent to the sum of deleterious allele frequency per gene)
eurbygene<-read.table("euroByGene.tab",header=TRUE) ################################### running as of Aug 8  ###################################

#######################################################################################################
####################        Subsample simulated datasets to match ExAC's n          ###################
#######################################################################################################

sims_AverUPOI<-as.data.frame(apply(sims_AverUExons, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sims_AverUPOI)<-c("V1")
sim1MPOI<-as.data.frame(apply(sim1M, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sim1MPOI)<-c("V1")
sim2MPOI<-as.data.frame(apply(sim2M, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sim2MPOI)<-c("V1")
sim5MPOI<-as.data.frame(apply(sim5M, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sim5MPOI)<-c("V1")
sim10MPOI<-as.data.frame(apply(sim10M, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sim10MPOI)<-c("V1")
simIGPOI<-as.data.frame(apply(simIG, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(simIGPOI)<-c("V1")
eurbygenePOI<-as.data.frame(apply(eurbygene[,1:32], c(1,2), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
sims_LargerAverUPOI<-as.data.frame(apply(sims_LargerAverU, c(1), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))
colnames(sims_LargerAverUPOI)<-c("V1")
simsPOI<-as.data.frame(apply(sims_fourTypes[,1:4], c(1,2), function(x) (rpois(1,(round(x*512000*2))*32881/512000))/65762))

########################################################################################
################################      MSB models        ################################
########################################################################################

aver_u=1.5e-8 							          #for exons, following Nealle et al. 2012
FIN=function(N){aver_u*sqrt(2*pi*N)} 	#Here, N is the number of diploid individuals, so as the input used for the c++ code used for simulations
INF=sqrt(aver_u)
FINTENN=FIN(512000)
FIN20K=FIN(20000) 
size=function(q){((q/aver_u)^2)/(2*pi)}

############################################################################################
################################            Fig 1           ################################
################################    Comparing expectations  ################################
############################################################################################

#setEPS()
#postscript("Fig1_TennessenU_rescalled_exons.eps",width=7.5,height=7)
tiff(filename = "Fig1.tif", width = 10, height = 10, units = "cm",res = 450)
SIMmean<-mean(sims_AverUExons$V1[1:100000])
sims_AverULOG=replace(sims_AverUExons$V1,sims_AverUExons$V1==0,1e-7)
histSIM<-hist(log10(sims_AverULOG[1:100000]),plot = FALSE)
plot(histSIM$counts,log="y", type='h', lwd=7, lend=1,col="grey",axes=FALSE,xlab = "Log10(Allele frequency)",ylab = "Density",cex.lab=0.8)
axis(side=1,at=c(1,6,11,16,21,26,31),labels=c("0","-6","-5","-4","-3","-2","-1"),cex.axis=0.8)
axis.break(axis = 1,breakpos = 3)                               #depends on "library(plotrix)"
abline(v=((1-((log10(SIMmean)+5)*-1))*5)+6,col="red",lwd = 1)       #SIM
abline(v=((1-((log10(FIN20K)+5)*-1))*5)+6,col="green",lwd = 1)      #FIN
abline(v=((1-((log10(INF)+3)*-1))*5)+16,col="blue",lwd = 1)         #INF
axis(side=2,at=c(10,100,1000,10000,100000,1000000),labels=c("1e+01","1e+02","1e+03","1e+04","1e+05","1e+06"),cex.axis=0.8)     
dev.off()

##########################################################################################
######################               Fig 2             ###################################
######################      Comparing ExAC to SIM      ###################################
##########################################################################################

# K for each mutation type
obsA<-nrow(XCpGti)
obsB<-nrow(XCpGtv)
obsC<-nrow(XnonCpGti)
obsD<-nrow(XnonCpGtv)

# How many repetitions
x=100000

#Sampling K allele frequencies from simulated dataset for each mutation type
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

#Significance of the deviation (notice that exact value may vary from the main text because it is subject to stochascity of the sampling scheme)
perc.rank <- function(x, xo)  (1-( (length(x[x <= xo])+1) / (length(x)+1) ))*2
perc.rank(x = sig_CpGti, xo = mean(XCpGti$eur_freq) )
perc.rank(x = sig_CpGtv, xo = mean(XCpGtv$eur_freq) )
perc.rank(x = sig_nonCpGti, xo = mean(XnonCpGti$eur_freq) )
perc.rank(x = sig_nonCpGtv, xo = mean(XnonCpGtv$eur_freq) )

#Set data points for panel B
x<-mean(XCpGti$eur_freq)/mean(sig_CpGti)
y<-mean(XCpGtv$eur_freq)/mean(sig_CpGtv)
z<-mean(XnonCpGti$eur_freq)/mean(sig_nonCpGti)
k<-mean(XnonCpGtv$eur_freq)/mean(sig_nonCpGtv)
foldIncrease<-c(k,z,y,x)
x=1.12e-7
y=9.59e-9
z=6.18e-9
k=3.76e-9
mutationRates<-c(k,z,y,x)

#setEPS()
#postscript("Fig2_rescaled.eps",width=6.8,height=6.35)

tiff(filename = "Fig2.tif", width = 14, height = 18, units = "cm",res = 450)

layout(matrix(c(1,2,3,4,5,5,5,5), 4, 2, byrow = TRUE))

SIMsig<-hist(sig_CpGti,plot=F,breaks=24)
SIMsig$counts=SIMsig$counts+1
SIMsig$counts=log10(SIMsig$counts)
SIMsig$counts=SIMsig$counts/max(SIMsig$counts)
plot(SIMsig,ylab='Log10(Density)',xlab='Mean allele frequencies',main="CpGti (K=101)\np-val=0.59",cex.main=0.9,col="grey",yaxt = "n",cex.main=1.1)#,xlim = c(0,1e-4))
abline(v=mean(XCpGti$eur_freq),col="blue",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1),cex.axis=0.7)

SIMsig<-hist(sig_CpGtv,plot=F,breaks=24)
SIMsig$counts=SIMsig$counts+1
SIMsig$counts=log10(SIMsig$counts)
SIMsig$counts=SIMsig$counts/max(SIMsig$counts)
plot(SIMsig,ylab='Log10(Density)',xlab='Mean allele frequencies',main="CpGtv (K=13)\np-val=0.08",cex.main=0.9,col="grey",yaxt = "n",cex.main=1.1)#,xlim = c(0,1e-4))
abline(v=mean(XCpGtv$eur_freq),col="blue",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1),cex.axis=0.7)

SIMsig<-hist(sig_nonCpGti,plot=F,breaks=24)
SIMsig$counts=SIMsig$counts+1
SIMsig$counts=log10(SIMsig$counts)
SIMsig$counts=SIMsig$counts/max(SIMsig$counts)
xuplim=mean(XnonCpGti$eur_freq)*1.1
plot(SIMsig,ylab='Log10(Density)',xlab='Mean allele frequencies',main="nonCpGti (K=155)\np-val<1e-4",col="grey",cex.main=0.9,yaxt = "n", xlim = c(0,xuplim),cex.main=1.1)
abline(v=mean(XnonCpGti$eur_freq),col="blue",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1),cex.axis=0.7)

SIMsig<-hist(sig_nonCpGtv,plot=F,breaks=24)
SIMsig$counts=SIMsig$counts+1
SIMsig$counts=log10(SIMsig$counts)
SIMsig$counts=SIMsig$counts/max(SIMsig$counts)
xuplim=mean(XnonCpGtv$eur_freq)*1.1
plot(SIMsig,ylab='Log10(Density)',xlab='Mean allele frequencies',main="nonCpGtv (K=104)\np-val<1e-4",col="grey",cex.main=0.9,yaxt = "n",cex.main=1.1)#, xlim = c(0,xuplim))
abline(v=mean(XnonCpGtv$eur_freq),col="blue",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1),cex.axis=0.7)

plot(log10(mutationRates),foldIncrease,pch=19,xlab="Log10(u)",ylab="Fold increase",cex.axis=1.2,cex.lab=1.3)#,main="B") 
labels=c("\nnonCpGtv","nonCpGti","CpGtv","")
text(log10(mutationRates),foldIncrease, labels=labels, cex= 1.2, pos=4,col=c("black","black","black","white"))
labels=c("","","","CpGti\n")
text(log10(mutationRates),foldIncrease, labels=labels, cex= 1.2, pos=2,col=c("white","white","white","black"))

mtext("A", side = 1, adj=0, line=-48, cex=1.3, font=2)
mtext("B", side = 1, adj=0, line=-20, cex=1.3, font=2)

dev.off()

#############################################################################################
###########################             Figure 3            	 ############################
###########################        Gene level analysis        	 ############################
#############################################################################################

geneExpEur<-colMeans(eurbygenePOI) 							
aggGeneList<-aggregate(. ~ gene, data=ExAC_noGCI,FUN=sum)
geneTable<-rbind(geneExpEur,aggGeneList$eur_freq) 			
geneTable<-t(geneTable) 
colnames(geneTable)<-c("ExpEur","ObsEur")
geneTable<-as.data.frame(geneTable) 						

#Significance of the deviation 

perc.rank <- function(x, xo)  (1-(length(x[x <= xo])/length(x)))*2
perc.rank(x = eurbygenePOI$ASPA, xo = geneTable[1,2] )
perc.rank(x = eurbygenePOI$ASS1, xo = geneTable[2,2] )
perc.rank(x = eurbygenePOI$CFTR, xo = geneTable[4,2] )
perc.rank(x = eurbygenePOI$CLN5, xo = geneTable[5,2] )
perc.rank(x = eurbygenePOI$DHCR7, xo = geneTable[6,2] )
perc.rank(x = eurbygenePOI$ERCC8, xo = geneTable[7,2] )
perc.rank(x = eurbygenePOI$FAH, xo = geneTable[8,2] )
perc.rank(x = eurbygenePOI$GAA, xo = geneTable[10,2] )
perc.rank(x = eurbygenePOI$GALC, xo = geneTable[11,2] )
perc.rank(x = eurbygenePOI$GAN, xo = geneTable[12,2] )
perc.rank(x = eurbygenePOI$GBE1, xo = geneTable[13,2] )
perc.rank(x = eurbygenePOI$HEXA, xo = geneTable[14,2] )
perc.rank(x = eurbygenePOI$HSD17B4, xo = geneTable[15,2] )
perc.rank(x = eurbygenePOI$IDUA, xo = geneTable[16,2] )
perc.rank(x = eurbygenePOI$LAMB3, xo = geneTable[18,2] )
perc.rank(x = eurbygenePOI$NPC1, xo = geneTable[19,2] )
perc.rank(x = eurbygenePOI$PEX7, xo = geneTable[20,2] )
perc.rank(x = eurbygenePOI$POLG, xo = geneTable[22,2] )
perc.rank(x = eurbygenePOI$POMGNT1, xo = geneTable[23,2] )
perc.rank(x = eurbygenePOI$PPT1, xo = geneTable[24,2] )
perc.rank(x = eurbygenePOI$PRF1, xo = geneTable[25,2] )
perc.rank(x = eurbygenePOI$SLC22A5, xo = geneTable[26,2] )
perc.rank(x = eurbygenePOI$SMARCAL1, xo = geneTable[27,2] )
perc.rank(x = eurbygenePOI$SMPD1, xo = geneTable[28,2] )
perc.rank(x = eurbygenePOI$STAR, xo = geneTable[29,2] )
perc.rank(x = eurbygenePOI$TK2, xo = geneTable[31,2] )
perc.rank(x = eurbygenePOI$TPP1, xo = geneTable[32,2] )

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

#setEPS()
#postscript("Fig3.eps",width=6,height=3.5)
tiff(filename = "Fig3.tif", width = 16, height = 10, units = "cm",res = 450)
par(mar=c(5.1,4.1,4.1,5.1))

plot(0:1,0:1,type="n",xlim=c(1,27),ylim=c(-6,-2),axes=FALSE,ann=FALSE)
vioplot(log10(CFTRObs),log10(IDUAObs),log10(POLGObs),log10(SMPD1Obs),log10(LAMB3Obs),log10(FAHObs),
        log10(PPT1Obs),log10(ASPAObs),log10(ERCC8Obs),log10(DHCR7Obs),log10(HSD17B4Obs),
        log10(PEX7Obs),log10(PRF1Obs),log10(ASS1Obs),log10(SMARCAL1Obs),log10(GBE1Obs),log10(NPC1Obs),log10(GALCObs),
        log10(TK2Obs),log10(TPP1Obs),log10(CLN5Obs),log10(GAAObs),log10(STARObs),
        log10(SLC22A5Obs),log10(GANObs),log10(HEXAObs),log10(POMGNT1Obs),col=c("grey"),drawRect=FALSE,add=TRUE)

axis(1, at=seq(1, 27, by=1), labels = FALSE)
text(seq(1, 27, by=1), par("usr")[3] - 0.2,
     labels = c("","","","","","","","ASPA","ERCC8","DHCR7","HSD17B4","PEX7","PRF1","ASS1","SMARCAL1",
                "GBE1","NPC1","GALC","TK2","TPP1","CLN5","GAA","STAR","SLC22A5","GAN","HEXA","POMGNT1"), srt = 45, pos = 1, xpd = TRUE,cex = 0.6)
text(seq(1, 27, by=1), par("usr")[3] - 0.2, labels = c("CFTR","IDUA","POLG","SMPD1","LAMB3","FAH","PPT1","","","","","","","","","","","","","","","","","","","",""), srt = 45, pos = 1, xpd = TRUE,cex = 0.6,font=2)
axis(side=2,at=c(-5.75,-5:-2),labels=c(0,-5:-2))
axis.break(axis = 2,breakpos = -5.25)
title(ylab="Log10(Frequency)")

points(x=1,y=log10(geneTable[4,2]),pch="-",cex=2,col="darkmagenta")
points(x=2,y=log10(geneTable[16,2]),pch="-",cex=2,col="darkmagenta")
points(x=3,y=log10(geneTable[22,2]),pch="-",cex=2,col="darkmagenta")
points(x=4,y=log10(geneTable[28,2]),pch="-",cex=2,col="darkmagenta")
points(x=5,y=log10(geneTable[18,2]),pch="-",cex=2,col="darkmagenta")
points(x=6,y=log10(geneTable[8,2]),pch="-",cex=2,col="darkmagenta")
points(x=7,y=log10(geneTable[24,2]),pch="-",cex=2,col="darkmagenta")
points(x=8,y=log10(geneTable[1,2]),pch="-",cex=2,col="darkmagenta")
points(x=9,y=log10(geneTable[7,2]),pch="-",cex=2,col="darkmagenta")
points(x=10,y=log10(geneTable[6,2]),pch="-",cex=2,col="darkmagenta")
points(x=11,y=log10(geneTable[15,2]),pch="-",cex=2,col="darkmagenta")
points(x=12,y=log10(geneTable[20,2]),pch="-",cex=2,col="darkmagenta")
points(x=13,y=log10(geneTable[25,2]),pch="-",cex=2,col="darkmagenta")
points(x=14,y=log10(geneTable[2,2]),pch="-",cex=2,col="darkmagenta")
points(x=15,y=log10(geneTable[27,2]),pch="-",cex=2,col="darkmagenta")
points(x=16,y=log10(geneTable[13,2]),pch="-",cex=2,col="darkmagenta")
points(x=17,y=log10(geneTable[19,2]),pch="-",cex=2,col="darkmagenta")
points(x=18,y=log10(geneTable[11,2]),pch="-",cex=2,col="darkmagenta")
points(x=19,y=log10(geneTable[31,2]),pch="-",cex=2,col="darkmagenta")
points(x=20,y=log10(geneTable[32,2]),pch="-",cex=2,col="darkmagenta")
points(x=21,y=log10(geneTable[5,2]),pch="-",cex=2,col="darkmagenta")
points(x=22,y=log10(geneTable[10,2]),pch="-",cex=2,col="darkmagenta")
points(x=23,y=log10(geneTable[29,2]),pch="-",cex=2,col="darkmagenta")
points(x=24,y=log10(geneTable[26,2]),pch="-",cex=2,col="darkmagenta")
points(x=25,y=log10(geneTable[12,2]),pch="-",cex=2,col="darkmagenta")
points(x=26,y=log10(geneTable[14,2]),pch="-",cex=2,col="darkmagenta")
points(x=27,y=log10(geneTable[23,2]),pch="-",cex=2,col="darkmagenta")
    
rect(0.6,-6,1.4,CFTRzeros,col="grey")
rect(1.6,-6,2.4,IDUAzeros,col="grey")
rect(2.6,-6,3.4,POLGzeros,col="grey")
rect(3.6,-6,4.4,SMPD1zeros,col="grey")
rect(4.6,-6,5.4,LAMB3zeros,col="grey")
rect(5.6,-6,6.4,FAHzeros,col="grey")
rect(6.6,-6,7.4,PPT1zeros,col="grey")
rect(7.6,-6,8.4,ASPAzeros,col="grey")
rect(8.6,-6,9.4,ERCC8zeros,col="grey")
rect(9.6,-6,10.4,DHCR7zeros,col="grey")
rect(10.6,-6,11.4,HSD17B4zeros,col="grey")
rect(11.6,-6,12.4,PEX7zeros,col="grey")
rect(12.6,-6,13.4,PRF1zeros,col="grey")
rect(13.6,-6,14.4,ASS1zeros,col="grey")
rect(14.6,-6,15.4,SMARCAL1zeros,col="grey")
rect(15.6,-6,16.4,GBE1zeros,col="grey")
rect(16.6,-6,17.4,NPC1zeros,col="grey")
rect(17.6,-6,18.4,GALCzeros,col="grey")
rect(18.6,-6,19.4,TK2zeros,col="grey")
rect(19.6,-6,20.4,TPP1zeros,col="grey")
rect(20.6,-6,21.4,CLN5zeros,col="grey")
rect(21.6,-6,22.4,GAAzeros,col="grey")
rect(22.6,-6,23.4,STARzeros,col="grey")
rect(23.6,-6,24.4,SLC22A5zeros,col="grey")
rect(24.6,-6,25.4,GANzeros,col="grey")
rect(25.6,-6,26.4,HEXAzeros,col="grey")
rect(26.6,-6,27.4,POMGNT1zeros,col="grey")
axis(side=4,at=c(-6,-5.75,-5.5,-5.25,-5),labels=c("0","25%","50%","75%","100%"),las=1,cex.axis=0.7)
dev.off()

############################################################################################
################################        Figure S1           ################################
################################  Compare Counsyl and ExAC  ################################
############################################################################################

#setEPS()
#postscript("FigS1.eps",width=6,height=6)
tiff(filename = "FigS1.tif", width = 10, height = 10, units = "in",res = 450)
plot(counsyl$freqEur_Counsyl,counsyl$EuroExAC,xlim = c(0,max(counsyl$freqEur_Counsyl)),ylim=c(0,max(counsyl$freqEur_Counsyl)),xlab="Allele frequency in Counsyl",ylab="Allele frequency in ExAC",pch=16,cex.lab=1.6,cex.axis=1.6)
abline(0,1,col=4,lty=2)
dev.off()          

###########################################################################################
#############################               Fig S2           ##############################
#############################   Comparing FIN and Tennessen  ##############################
###########################################################################################

#setEPS()
#postscript("FigS2.eps",width=5.5,height=6.7)
tiff(filename = "FigS2.tif", width = 5.5, height = 6.7, units = "in",res = 450)
par(mfrow=c(3,1))

SmallN<-14376 #This is the smallest N used for SIM rescalling considering the difference in mutation rate
LargeN<-1006933 #This is the largest N used for SIM (corresponds to current European pop size) rescalling considering the difference in mutation rate
plot(0:1,0:1,type="n",xlim=c(0,LargeN),ylim=c(0,FIN(LargeN)),axes=FALSE,ann=FALSE)
plot(FIN,SmallN,LargeN,add=TRUE)
axis(side=1)
axis(side=2)
title(ylab="Deleterious allele frequency",xlab="Population size")
abline(v=size(SIMmean),col="red")
title(outer=F,adj=0,main="A",cex=10,col="black",font=4,line=1)
HconstSize<-hist(constSize$X0,plot=F,breaks=28)
Hsims_AverU<-hist(sims_AverUExons$V1[1:100000],plot=F,breaks = 14)
Hsims_AverU$counts=Hsims_AverU$counts+1
Hsims_AverU$counts=log10(Hsims_AverU$counts)
Hsims_AverU$counts=Hsims_AverU$counts/max(Hsims_AverU$counts)
plot(Hsims_AverU,ylab='Log10(Density)',xlab='Allele frequency',main="",col="grey",yaxt = "n",xlim=c(0,0.007))
title(outer=FALSE,adj=0,main="B",cex=2,col="black",font=2,line=1)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1))
HconstSize<-hist(constSize$X0[1:100000],plot=F,breaks=28)
HconstSize$counts=HconstSize$counts+1
HconstSize$counts=log10(HconstSize$counts)
HconstSize$counts=HconstSize$counts/max(HconstSize$counts)
plot(HconstSize,ylab='Log10(Density)',xlab='Allele frequency',main="",col="grey",yaxt = "n",xlim=c(0,0.007))
title(outer=FALSE,adj=0,main="C",cex=2,col="black",font=2,line=1)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1))

ks.test(sims_AverUExons$V1[1:100000],constSize$X0[1:100000])

dev.off()

##########################################################################################
######################             Fig S3             	##################################
######################      Selective effect in hets    ##################################
##########################################################################################

#setEPS()
#postscript("FigS3.eps",width=5.5,height=6.7)
tiff(filename = "FigS3.tif", width = 5.5, height = 6.7, units = "in",res = 450)
par(mfrow=c(2,1))

sims_AverULOG=replace(sims_AverUExons$V1,sims_AverUExons$V1==0,0.5e-6)
Hsims_AverU<-hist(log10(sims_AverULOG[1:100000]),plot=F)
Hsims_AverU$counts=Hsims_AverU$counts+1
Hsims_AverU$counts=log10(Hsims_AverU$counts)
Hsims_AverU$counts=Hsims_AverU$counts/max(Hsims_AverU$counts)
plot(Hsims_AverU,ylab='Log10(Density)',xlab='Log10(Allele frequency)',main="h = 0",col="grey",yaxt = "n")
abline(v=log10(mean(sims_AverUExons$V1)),col="red",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1))
sims_EffectHetLOG=replace(sims_EffectHet$V1,sims_EffectHet$V1==0,0.5e-6)
Hsims_EffectHet<-hist(log10(sims_EffectHetLOG),plot=F) 
Hsims_EffectHet$counts=Hsims_EffectHet$counts+1
Hsims_EffectHet$counts=log10(Hsims_EffectHet$counts)
Hsims_EffectHet$counts=Hsims_EffectHet$counts/max(Hsims_EffectHet$counts)
plot(Hsims_EffectHet,ylab='Log10(Density)',xlab='Log10(Allele frequency)',main="h = 1%",col="grey",yaxt = "n",xlim=c(-6.301030, -2.258098))
abline(v=log10(mean(sims_EffectHet$V1)),col="red",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1))

dev.off()

##############################################################################################
######################                Fig S4                ##################################
######################      Consider more growth for SIM    ##################################
##############################################################################################

howManyObs <- function(x)  x[which(x != 0)]
Tobs<-howManyObs(sim1MPOI$V1)
Bobs<-howManyObs(sim2MPOI$V1)
Cobs<-howManyObs(sim5MPOI$V1)
Dobs<-howManyObs(sim10MPOI$V1)
Eobs<-howManyObs(simIGPOI$V1)
Fobs<-howManyObs(sims_LargerAverUPOI$V1)
eobs<-howManyObs(ExAC$eur_freq)

howManyZeros <- function(x, xo) (((length(x)-length(xo))/length(x))*1)-7
Tzero<-howManyZeros(sim1MPOI$V1,Tobs)
Bzero<-howManyZeros(sim2MPOI$V1,Bobs)
Czero<-howManyZeros(sim5MPOI$V1,Cobs)
Dzero<-howManyZeros(sim10MPOI$V1,Dobs)
Ezero<-howManyZeros(simIGPOI$V1,Eobs)
Fzero<-howManyZeros(sims_LargerAverUPOI$V1,Fobs)
ezero<-howManyZeros(ExAC$eur_freq,eobs)


#setEPS()
#postscript("FigS4.eps",width=6,height=7)
tiff(filename = "FigS4.tif", width = 6, height = 7, units = "in",res = 450)
plot(0:1,0:1,type="n",xlim=c(0.5,7.5),ylim=c(-7,-2),axes=FALSE,ann=FALSE)
vioplot(log10(Tobs),log10(Bobs),log10(Cobs),log10(Dobs),log10(Eobs),log10(Fobs),log10(eobs),drawRect=FALSE,col=c("grey","grey","grey","grey","grey","grey","grey","white"),add=TRUE)
axis(side=1,at=1:7,labels=c("A","B","C","D","E","F","ExAC"))
axis(side=2,at=c(-6.65,-5:-2),labels=c(0,-5:-2))
axis.break(axis = 2,breakpos = -5.7)
title(ylab="Log10(Frequency)")
rect(0.6,-7,1.4,Tzero,col="grey")
rect(1.6,-7,2.4,Bzero,col="grey")
rect(2.6,-7,3.4,Czero,col="grey")
rect(3.6,-7,4.4,Dzero,col="grey")
rect(4.6,-7,5.4,Ezero,col="grey")
rect(5.6,-7,6.4,Fzero,col="grey")
rect(6.6,-7,7.4,ezero,col="white")
points(x=1,y=log10(mean(sim1MPOI$V1)), pch="-",cex=2,col="red")
points(x=2,y=log10(mean(sim2MPOI$V1)), pch="-", cex=2,col="red")
points(x=3,y=log10(mean(sim5MPOI$V1)), pch="-", cex=2,col="red")
points(x=4,y=log10(mean(sim10MPOI$V1)), pch="-", cex=2,col="red")
points(x=5,y=log10(mean(simIGPOI$V1)), pch="-", cex=2,col="red")
points(x=6,y=log10(mean(sims_LargerAverUPOI$V1)), pch="-", cex=2,col="red")
points(x=7,y=log10(mean(ExAC$eur_freq)), pch="-", cex=2,col="red")

x<-format(round(((2000000-(length(Tobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(1,-6.65,paste(x,"%",sep=""),cex=1)

x<-format(round(((2000000-(length(Bobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(2,-6.65,paste(x,"%",sep=""),cex=1)

x<-format(round(((2000000-(length(Cobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(3,-6.65,paste(x,"%",sep=""),cex=1)

x<-format(round(((2000000-(length(Dobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(4,-6.65,paste(x,"%",sep=""),cex=1)

x<-format(round(((2000000-(length(Eobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(5,-6.65,paste(x,"%",sep=""),cex=1)

x<-format(round(((2000000-(length(Fobs)))/2000000), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(6,-6.65,paste(x,"%",sep=""),cex=1)

x<-format(round(((nrow(ExAC)-(length(eobs)))/nrow(ExAC)), 2), nsmall = 2)
x<-as.numeric(x)
x<-x*100
text(7,-6.73,paste(x,"%",sep=""),cex=1)

dev.off()

############################################################################################
######################                Fig S5                	############################
######################      (equivalent to Fig 2 w/o 2 genes    ############################
############################################################################################

exac_wo2genes<-subset(ExAC_noGCI,gene!="CFTR" & gene!="DHCR7")
XnonCpGtv <- subset(exac_wo2genes,type=="nonCpGtv")
XnonCpGti <- subset(exac_wo2genes,type=="nonCpGti")
XCpGtv <- subset(exac_wo2genes,type=="CpGtv")
XCpGti <- subset(exac_wo2genes,type=="CpGti")

# K 
obsA<-nrow(XCpGti)
obsB<-nrow(XCpGtv)
obsC<-nrow(XnonCpGti)
obsD<-nrow(XnonCpGtv)

#How many repetitions
x=100000

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

#Significance of the deviation
perc.rank <- function(x, xo)  (1-(length(x[x <= xo])/length(x)))*2
perc.rank(x = sig_CpGti, xo = mean(XCpGti$eur_freq) )   # 0.79808
perc.rank(x = sig_CpGtv, xo = mean(XCpGtv$eur_freq) )   # 0.06582
perc.rank(x = sig_nonCpGti, xo = mean(XnonCpGti$eur_freq) ) # 0.00082
perc.rank(x = sig_nonCpGtv, xo = mean(XnonCpGtv$eur_freq) ) # 8e-05

#setEPS()
#postscript("FigS5.eps",width=7,height=7)
tiff(filename = "FigS5.tif", width = 7, height = 7, units = "in",res = 450)
par(mfrow=c(2,2))

SIMsig<-hist(sig_CpGti,plot=F,breaks=24)
br=SIMsig$breaks
SIMsig$counts=SIMsig$counts+1
SIMsig$counts=log10(SIMsig$counts)
SIMsig$counts=SIMsig$counts/max(SIMsig$counts)
plot(SIMsig,ylab='Log10(Density)',xlab='Mean allele frequencies',main="CpGti (K=80)\np-val=0.80",cex.main=0.9,col="grey",yaxt = "n")#,xlim = c(0,1e-4))
abline(v=mean(XCpGti$eur_freq),col="blue",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1))
 
SIMsig<-hist(sig_CpGtv,plot=F,breaks=24)
SIMsig$counts=SIMsig$counts+1
SIMsig$counts=log10(SIMsig$counts)
SIMsig$counts=SIMsig$counts/max(SIMsig$counts)
plot(SIMsig,ylab='Log10(Density)',xlab='Mean allele frequencies',main="CpGtv (K=9)\np-val=0.07",cex.main=0.9,col="grey",yaxt = "n")#,xlim = c(0,1e-4))
abline(v=mean(XCpGtv$eur_freq),col="blue",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1))

SIMsig<-hist(sig_nonCpGti,plot=F,breaks=24)
SIMsig$counts=SIMsig$counts+1
SIMsig$counts=log10(SIMsig$counts)
SIMsig$counts=SIMsig$counts/max(SIMsig$counts)
xuplim=mean(XnonCpGti$eur_freq)*1.1
plot(SIMsig,ylab='Log10(Density)',xlab='Mean allele frequencies',main="nonCpGti (K=83)\np-val<1e-2",col="grey",cex.main=0.9,yaxt = "n")# xlim = c(0,xuplim))
abline(v=mean(XnonCpGti$eur_freq),col="blue",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1))

SIMsig<-hist(sig_nonCpGtv,plot=F,breaks=24)
SIMsig$counts=SIMsig$counts+1
SIMsig$counts=log10(SIMsig$counts)
SIMsig$counts=SIMsig$counts/max(SIMsig$counts)
xuplim=mean(XnonCpGtv$eur_freq)*1.1
plot(SIMsig,ylab='Log10(Density)',xlab='Mean allele frequencies',main="nonCpGtv (K=59)\np-val<1e-3",col="grey",cex.main=0.9,yaxt = "n")# xlim = c(0,xuplim),yaxt = "n")
abline(v=mean(XnonCpGtv$eur_freq),col="blue",lwd = 2)
axis(2,labels=c(0,1,10,100),at=c(0,0.33,0.66,1))

dev.off()

########################################################################################
##################        Simulations of ascertainment bias        #####################
##################        			and Fig S6   				               #####################
########################################################################################

# A function that simulates genotypes with certain inbreeding coefficient for a mutation at frequency "freq" 
geno_sim_w_inbreeding = function(freq, n_ind, inbreed_coef){
  allele_1 = rbinom(n_ind, size=1, prob=freq)
  inbreeding_prob= runif(n_ind, min=0, max=1)
  allele_2 = numeric()
  allele_2[which(inbreeding_prob<=inbreed_coef)] <- allele_1[which(inbreeding_prob<=inbreed_coef)] # a random individual is IBD with probability of inbreed_coef
  allele_2[which(inbreeding_prob>inbreed_coef)] <- rbinom(sum(inbreeding_prob>inbreed_coef), size=1, prob=freq)
  sim_geno<- allele_1+allele_2
  return(c(sum(sim_geno==2),sum(sim_geno==1),sum(sim_geno==0)))
}

# A function that simulates the ascertainment status for a vector of frequencies
ascertain_or_not = function(freqs, n, f){
  n_homo <- sapply(freqs, function(x) geno_sim_w_inbreeding(x, n, f)[1])
  return(n_homo>=1) # If and only if there is at least one homozygote, the mutation is ascertained
}

# q represents the output of forward simulations (i.e., a vector of population frequencies of deleterious alleles)
qa <- sims_fourTypes$CpGti
qb <- sims_fourTypes$CpGtv
qc <- sims_fourTypes$nonCpGti
qd <- sims_fourTypes$nonCpGtv  

qa <-qa[1:100000]
qb <-qb[1:100000]
qc <-qc[1:100000]
qd <-qd[1:100000]

# Set parameters for ascertainment set
n_asc=100000 #32881 # number of INDIVIDUALS in the ascertainment set
f_asc=1/6 # average inbreeding coefficient for the ascertainment set, 1/16 corresponds to offspring of first cousins
ascertain_statusa <- ascertain_or_not(qa, n_asc, f_asc)
ascertain_statusb <- ascertain_or_not(qb, n_asc, f_asc)
ascertain_statusc <- ascertain_or_not(qc, n_asc, f_asc)
ascertain_statusd <- ascertain_or_not(qd, n_asc, f_asc)

# Set parameters for ExAC set
n_ExAC=65762

n_allele_ExACa = sapply(qa, function(x) rbinom(n_ExAC, n=1, prob=x))
n_allele_ExACb = sapply(qb, function(x) rbinom(n_ExAC, n=1, prob=x))
n_allele_ExACc = sapply(qc, function(x) rbinom(n_ExAC, n=1, prob=x))
n_allele_ExACd = sapply(qd, function(x) rbinom(n_ExAC, n=1, prob=x))

freq_ExACa = n_allele_ExACa/n_ExAC
freq_ExACb = n_allele_ExACb/n_ExAC
freq_ExACc = n_allele_ExACc/n_ExAC
freq_ExACd = n_allele_ExACd/n_ExAC

# Compare results with or without ascertainment bias
# Distribution of allele frequencies without ascertainment bias 
summary(freq_ExAC)
sum(freq_ExAC==0)/length(freq_ExAC)
mean(freq_ExAC,na.rm=T)
# Distribution of allele frequencies conditional on being ascertained
summary(freq_ExAC[ascertain_status])
sum(freq_ExAC[ascertain_status]==0)/length(freq_ExAC[ascertain_status])
mean(freq_ExAC[ascertain_status],na.rm=T)
# Ascertainment probability
sum(ascertain_statusa)/length(ascertain_statusa)
sum(ascertain_statusb)/length(ascertain_statusb)
sum(ascertain_statusc)/length(ascertain_statusc)
sum(ascertain_statusd)/length(ascertain_statusd)
# Fold change in allele frequency due to ascertainment bias
mean(freq_ExACa[ascertain_statusa],na.rm=T)/mean(freq_ExACa,na.rm=T)
mean(freq_ExACb[ascertain_statusb],na.rm=T)/mean(freq_ExACb,na.rm=T)
mean(freq_ExACc[ascertain_statusc],na.rm=T)/mean(freq_ExACc,na.rm=T)
mean(freq_ExACd[ascertain_statusd],na.rm=T)/mean(freq_ExACd,na.rm=T)


### Function to calculate the probability of being ascertained given frequency of mutation, sample size and inbreeding coeficient

fbase=1/16 #first cousins 
freqbase=mean(sims_AverUExons$V1)
nbase=10000

p_asc = function (freq,n_ind,inbreed_coef){
  prob_asc = 1-((1-freq)*(1+freq-(freq*inbreed_coef)))^n_ind
  return(prob_asc)
}
p_ascA = function (freq){
  prob_asc = 1-((1-freq)*(1+freq-(freq*1/16)))^32881
  return(prob_asc)
}
p_ascB = function (n_ind){
  prob_asc = 1-((1-SIMmean)*(1+SIMmean-(SIMmean*1/16)))^n_ind
  return(prob_asc)
}
p_ascC = function (inbreed_coef){
  prob_asc = 1-((1-SIMmean)*(1+SIMmean-(SIMmean*inbreed_coef)))^32881
  return(prob_asc)
}

#setEPS()
#postscript("FigS6.eps",width=8,height=3)
tiff(filename = "FigS6.tif", width = 8, height = 3, units = "in",res = 450)
par(mfrow=c(1,3))
plot.function(p_ascA,min(sims_AverUExons$V1),max(sims_AverUExons$V1),xlab="q",ylab=expression('P'['asc']))
plot.function(p_ascB,100,5000000,xlab=expression('n'['a']),ylab=expression('P'['asc']))
plot.function(p_ascC,0,1/2,xlab=expression('F'['a']),ylab=expression('P'['asc']))
dev.off()

#######################################################################
##################              Coverage          #####################
##################         	     Fig S7   		  #####################
#######################################################################

#setEPS()
#postscript("FigS7.eps",width=6,height=7)
tiff(filename = "FigS7.tif", width = 6, height = 7, units = "in",res = 450)
zero<-read.table("zeros_coverage.txt",header=FALSE)
nonzero<-read.table("nonZeros_coverage.txt",header=FALSE)
data<-data.frame(zero)
boxplot(zero,at=2,xlim=c(0,3))
boxplot(nonzero,at=1,add=TRUE)
axis(side = 1,at = c(1,2),labels = c("A","B"))
title(ylab="Coverage")
mean(zero$V1)
mean(nonzero$V1)
dev.off()
