#..........................................#
#........Turnover_dendritic_network........#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt                                                    #
#... Author of the script: Eric Harvey & Isabelle Gounand  - contact: eric.harvey@eawag.ch                                                                                                   #
#..........................................................................................................................................#

#########################################################################
################# DATA STRUCTURE
#########################################################################

##################
#Clear any variables from the R environment 
##################
rm(list=ls())


##################
#Load packages
##################
library(vegan)
library(plyr)
library(sciplot)
library(nlme)
library(Hmisc)

##################
#Directory paths 
##################
datapath = "~/Documents/Research/Eawag/Projects/8.Blue-Green/Turnover_dendritic_network/Experiment/"
graphpath = "~/Documents/Research/Eawag/Projects/8.Blue-Green/4.Results/"

##################
#Load data in R environment
##################
Prot =  read.delim(paste0(datapath,"Turnover_dendritic_network_(20170320).txt")) #Protist

##################
#Fix Protist data structure
##################
Prot$date = with(Prot,revalue(date, c("16-05-02"="20160502", "16-05-23"="20160523","16-05-31"="20160531","5/17/2016"="20160517","5/9/2016"="20160509")))
Prot$date = factor(Prot$date,levels=c("20160502","20160509","20160517","20160523","20160531"))
Prot = Prot[order(Prot$date),]
Prot$day = c(rep(0,22),rep(7,158),rep(15,158),rep(21,158),rep(29,158))
Prot = Prot[which(Prot$Treatment=="Isolated" & Prot$day!=0),]


##################
#Calculate protist TOTAL abundance per microcosm (as opposed to density per volume)
##################
colnames(Prot)[23:29] = c("Rot.ul","Spi.ul","Ble.ul","Pca.ul","Col.ul","Chi.ul","Tet.ul")
Prot$Rot.all = Prot$Rot.ul*Prot$Size*1000
Prot$Spi.all = Prot$Spi.ul*Prot$Size*1000
Prot$Ble.all = Prot$Ble.ul*Prot$Size*1000
Prot$Pca.all = Prot$Pca.ul*Prot$Size*1000
Prot$Col.all = Prot$Col.ul*Prot$Size*1000
Prot$Chi.all = Prot$Chi.ul*Prot$Size*1000
Prot$Tet.all = Prot$Tet.ul*Prot$Size*1000
Prot$bioarea.all = Prot$bioarea_per_ul*Prot$Size*1000
Prot$Prot.tot.ab = rowSums(Prot[,32:38])

#Turn Size back to a factor
Prot$Size = as.factor(Prot$Size)
Prot = droplevels(Prot)

##################
#Calculate protist diversity
##################
Prot$Prot.rich = specnumber(Prot[,32:38])
date.rep.int = interaction(Prot$day,Prot$Replicate)
(gamma = specnumber(Prot[,32:38],groups=date.rep.int)) #gamma diversity does not chagne


#########################################################################
################# Figures & Analyses  
#########################################################################

##################
# Figures
##################

#Protist richness as a function of patch size for each experimental day
lineplot.CI(Size,Prot.rich,day,data=Prot,xlab="Experimental days",ylab="Species richness")
#Protist total abundance as a function of patch size for each experimental day
lineplot.CI(Size,Prot.tot.ab,day,data=Prot,xlab="Experimental days",ylab="Species richness")

##################
# Stats
##################
Mod = lme(Prot.rich ~ Size*day, ~ date|Replicate,data=Prot,method="ML",control=lmeControl(optimMethod="REML",maxIter=100,opt="optim"))
summary(Mod)$tTable
plot(density(Mod$residuals))

#Summary information 
#....Species richness Mean +- SD for day 7
tapply(Prot$Prot.rich[which(Prot$day==7)],Prot$Size[which(Prot$day==7)],mean)
tapply(Prot$Prot.rich[which(Prot$day==7)],Prot$Size[which(Prot$day==7)],sd)
#.....Species richness Mean +- SD for day 15 and 21
tapply(Prot$Prot.rich[which(Prot$day==15)],Prot$Size[which(Prot$day==15)],mean)
tapply(Prot$Prot.rich[which(Prot$day==15)],Prot$Size[which(Prot$day==15)],sd)
tapply(Prot$Prot.rich[which(Prot$day==21)],Prot$Size[which(Prot$day==21)],mean)
tapply(Prot$Prot.rich[which(Prot$day==21)],Prot$Size[which(Prot$day==21)],sd)
#.....Species richness Mean +- SD for day 29
tapply(Prot$Prot.rich[which(Prot$day==29)],Prot$Size[which(Prot$day==29)],mean)
tapply(Prot$Prot.rich[which(Prot$day==29)],Prot$Size[which(Prot$day==29)],sd)


##################
# Explore changes in beta-diversity
##################

Comp.mat1 = Prot[,32:38]
Comp.mat1[Comp.mat1>0] = 1

Prot.diss.null = raupcrick(Comp.mat1,nsimul=999) 
Prot.diss = vegdist(Comp.mat1,"jaccard") 

beta.null = betadisper(Prot.diss.null,Prot$Size)
beta = betadisper(Prot.diss,Prot$Size)


mean.effect.null = 0
mean.effect = 0
sd.effect.null = 0
sd.effect = 0
for(i in 1:10000){
  
  beta.effect = log(beta$distances[which(Prot$Size=="45")]/sample(beta$distances[which(Prot$Size=="7.5")],44))
  beta.effect.null = log(beta.null$distances[which(Prot$Size=="45")]/sample(beta.null$distances[which(Prot$Size=="7.5")],44))
  
  mean.effect[i] = mean(beta.effect)
  sd.effect[i] = sd(beta.effect)
  mean.effect.null[i] = mean(beta.effect.null)
  sd.effect.null[i] = sd(beta.effect.null)
}


error <- qnorm(0.975)*mean(sd.effect)/sqrt(44)
error.null <- qnorm(0.975)*mean(sd.effect.null)/sqrt(44)

#pdf(paste0(fig.path,"Beta_Alpha_Patch.pdf"),width=5,height = 5)
plot(mean.effect,type="n",ylab="Effect size",xaxt="n",ylim=c(-1,1))
abline(h=0,lwd=1,col="gray")
errbar(x=2000,mean(mean.effect),mean(mean.effect)+error,mean(mean.effect)-error,add=T,errbar.col="blue",col="blue",lty=1,pch=16)
errbar(x=8000,mean(mean.effect.null),mean(mean.effect.null)+error.null,mean(mean.effect.null)-error.null,add=T,errbar.col="orange",col="orange",lty=1,pch=16)


