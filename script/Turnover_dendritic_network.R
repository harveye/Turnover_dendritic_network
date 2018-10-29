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
library(tidyverse)
library(vegan)
library(sciplot)
library(nlme)
library(Hmisc)

##################
#Directory paths 
##################
to.data <- "./data/"
to.script <- "./scripts/"
to.output <- "./output/"
to.figs <- "./figs/"
to.R <- "./R/"

##################
#Load data in R environment
##################
Prot.raw =  read.delim(paste0(to.data,"Turnover_dendritic_network_(20170320).txt")) #Protist

##################
#Fix Protist data structure
##################
species = c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other") #place species names in a vector

Prot <- Prot.raw %>% 
  mutate(prot.rich = specnumber(Prot.raw[,species]),
         day = ifelse(date=="16-05-02",0, ifelse(date=="5/9/2016",7, ifelse(date=="5/17/2016",15,
                     ifelse(date=="16-05-23",21, ifelse(date=="16-05-31",29, NA))))),
         date= recode_factor(date,"16-05-02"="20160502",
                            "16-05-23"="20160523",
                            "16-05-31"="20160531",
                            "5/17/2016"="20160517",
                            "5/9/2016"="20160509"),
         Size = as.factor(Size)) %>% 
  filter(Treatment=="Isolated" & day !=0)

Prot <- droplevels(Prot)

##################
#Calculate protist diversity
##################

#Look for changes in gamma diversity
date.rep.int = interaction(Prot$day,Prot$Replicate)
(gamma = specnumber(Prot[,species],groups=date.rep.int)) #gamma diversity does not chagne


#########################################################################
################# Figures & Analyses  
#########################################################################

##################
# Figures
##################

#######
#Figure 2 (Experimental only)
X <- Prot %>% 
  dplyr::filter(Treatment=="Isolated", day==29)
  ggplot(mapping=aes(x=as.factor(Size),y=prot.rich,fill=Size),X) +
  geom_boxplot() 

#Information for figure 2 legend
# Samping unit (N) for each path size
length(Prot$prot.rich[which(Prot$day==29 & Prot$Size==7.5)])
length(Prot$prot.rich[which(Prot$day==29 & Prot$Size==13)])
length(Prot$prot.rich[which(Prot$day==29 & Prot$Size==22.5)])
length(Prot$prot.rich[which(Prot$day==29 & Prot$Size==45)])


#######
#Supplementary Figure 2
pdf(paste0(to.figs,"FigS2.pdf"),width=10,height=5)
#Protist richness as a function of patch size for each experimental day
lineplot.CI(Size,prot.rich,day,data=Prot,xlab="Experimental days",ylab="Species richness")
dev.off()


##################
# Stats
##################

#Model for all experimental days
Mod = lme(prot.rich~ Size*day, ~ date|Replicate,data=Prot,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
(table = summary(Mod)$tTable)
plot(density(Mod$residuals))
write.csv(table,paste0(to.output,"Mod.ALL.csv"))

#Model for last experimental day only (day 29 - used in all main figures)
Mod.day29 = lme(prot.rich~ Size, ~ 1|Replicate,data=Prot[which(Prot$day==29),],method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
(table.day29 = summary(Mod.day29)$tTable)
plot(density(Mod.day29$residuals))
write.csv(table.day29,paste0(to.output,"Mod.day29.csv"))

#Summary information (in text - Results section) 
#....Species richness Mean +- SD for day 7
tapply(Prot$prot.rich[which(Prot$day==7)],Prot$Size[which(Prot$day==7)],mean)
tapply(Prot$prot.rich[which(Prot$day==7)],Prot$Size[which(Prot$day==7)],sd)
#.....Species richness Mean +- SD for day 15 and 21
tapply(Prot$prot.rich[which(Prot$day==15)],Prot$Size[which(Prot$day==15)],mean)
tapply(Prot$prot.rich[which(Prot$day==15)],Prot$Size[which(Prot$day==15)],sd)
tapply(Prot$prot.rich[which(Prot$day==21)],Prot$Size[which(Prot$day==21)],mean)
tapply(Prot$prot.rich[which(Prot$day==21)],Prot$Size[which(Prot$day==21)],sd)
#.....Species richness Mean +- SD for day 29
tapply(Prot$prot.rich[which(Prot$day==29)],Prot$Size[which(Prot$day==29)],mean)
tapply(Prot$prot.rich[which(Prot$day==29)],Prot$Size[which(Prot$day==29)],sd)


##################
# Explore changes in composition (Fig. S11)
##################

#...Build the dendritic landscapes from connectivity matrices and calculate network metrics of interest
library(igraph)
source(paste0(to.R,"Network_metric_bg.R"))

#...Set the color palette for each species
library(RColorBrewer)
#display.brewer.all(n=7,type="qual",colorblindFriendly = T)
pal = brewer.pal(n=7,name="Dark2")
rbPal <- colorRampPalette(pal,space="Lab")

#...Add those colors to the properties of the igraph file
V(RepA.g)$pie.color=list(rbPal(7))
V(RepB.g)$pie.color=list(rbPal(7))
V(RepC.g)$pie.color=list(rbPal(7))
V(RepD.g)$pie.color=list(rbPal(7))

#...Factors to loop on
dates = as.factor(c(20160509,20160531))
day = as.factor(c(7,29))

pdf(paste0(to.figs,"PROT_NET_COMP",".pdf"),width=5,height=10)

layout(matrix(c(1,0,0,0,2:9),ncol=3,nrow=4,byrow=F),width=c(0.5,1,1,1))
#layout.show(n=9)

#..plot legends
par(mar = c(4,3,3,4)) #c(bottom, left, top, right) 
image(matrix(1:7,nrow=1,ncol=7), col = rbPal(7), axes =F)
axis(2, at = seq(0.001,0.999,length=7),font=2, labels =species[-8],tick=F,cex.axis=0.8)

#..plot pie chart networks
par(mar = c(1,1,1,1))
for(i in 1:nlevels(dates)) {
  for(j in 1:nlevels(Prot$Replicate)) {
    
    #.Extract density data for each species
    rel.ab = Prot[which(Prot$date == levels(dates)[i] & Prot$Treatment=="Isolated" & Prot$Replicate==levels(Prot$Replicate)[j]),species[-8]]
    rel.ab = rel.ab+0.00000000000000000000000000000001
    rel.ab.list <- as.list(as.data.frame(t(rel.ab)))
    #abundance = Prot.b$Pcount[which(Prot.b$date == levels(Prot.b$date)[i] & Prot.b$Treatment==levels(Prot.b$Treatment)[k] & Prot.b$Replicate==levels(Prot.b$Replicate)[j])]
    Psize = Prot$Size[which(Prot$date == levels(dates)[i] & Prot$Treatment=="Isolated" & Prot$Replicate==levels(Prot$Replicate)[j])]
    
    #.Plot the figure
    
    #.Generate figure layout corresponding to replicate landscape
    coords1 <- layout_on_grid(eval(parse(text= paste("Rep",levels(Prot$Replicate)[j],".","g",sep=""))),width=6,height=6) #set the rectangular layout
    coords1[,2] <- max(coords1[,2])-coords1[,2]
    
    #.Plot
    plot.igraph(eval(parse(text= paste("Rep",levels(Prot$Replicate)[j],".","g",sep=""))),
                layout=coords1,
                vertex.shape="pie",
                vertex.pie=rel.ab.list,
                vertex.size=20,#(abundance)*magn Use this if you want the size of the pie to be proportional to total abundance in the patch
                vertex.label=Psize,
                vertex.label.dist=2.2,
                vertex.label.degree=-pi/2,
                vertex.label.color="black",
                vertex.label.cex=0.8,
                vertex.label.font=2,
                main=paste0("Day"," ",levels(day)[i]," ; ","Rep"," ",levels(Prot$Replicate)[j]))
  }
}

dev.off()




#THE END#######


