#... Blue-Green Project                                              
#... co-authors: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt
#... Author of the script: Isabelle Gounand
#... date: 06-07-2017
#... Script description: functions for the analysis of the C output: plot the landscapes (of 6x6 patches)

#install.packages("plotrix")
library(plotrix)
#install.packages("mapplots")
library(mapplots)

#... Function to calculate the species richness at one time step 
#...............................................................
# x are the densities at one timestep (a vector with the abundances of each species x in each patch)
# p_num are the number of patches
# d_thresh is the density under which a species is considered extinct
# test values
#x = dyn.list[[1]][[1]][10,-1]
#p_num = p["np"]; d_tresh = p["death_thresh"]

get.patchRichness = function(x,p_num,d_tresh){
  m0 = matrix(as.numeric(x), ncol = p_num, byrow = F)
  m = ifelse(m0 < d_tresh,0,1)
  return(apply(m,2,sum))
}



#... Function to calculate abundance at one time step 
#....................................................
# x are the densities at one timestep (a vector with the abundances of each species x in each patch)
# p_num are the number of patches
# test values
# x = dyn.list[[1]][[1]][10,-1]
# p_num = p["np"];

get.patchAbundance = function(x,p_num){
  m = matrix(as.numeric(x), ncol = p_num, byrow = F)
  return(apply(m,2,sum))
}



#... Plot the vertices (connectivity between patches)
#.....................
# provide the connectivity matrix of size patch number ^2

plot.vertices = function(connectivity.mat)
{
  co = connectivity.mat
  co[upper.tri(co)]=0
  for(i in 1:nrow(co))for(j in 1:ncol(co))if(co[i,j]==1)segments(pTOL[i,"x"],pTOL[i,"y"],pTOL[j,"x"],pTOL[j,"y"])
}



#... Scale an array
#..................
# x0 the array to scale
# minmax the boundaries (array of 2 number) between which the array has to be scaled (uniform transformation)

scale.func=function(x0,minmax){
  ran.x=max(x0,na.rm=T)-min(x0,na.rm=T)
  ran.s=max(minmax,na.rm=T)-min(minmax,na.rm=T)
  x1=x0 * ran.s / ran.x
  return(x1+(min(minmax,na.rm=T)-min(x1,na.rm=T)))
}



#... Plot the patches with a colour gradient
#...........................................
# values: array of patch values to scale to the colour gradient (species richness)
# minmaxV: min and max values to scale to the colour extrema (array of 2 numbers)
# colour.grad: colour gradient
# r: radius of the patches

plot.patch.oneVariable = function(values,minmaxV,colour.grad,r=0.2)
{
  values.scaled = scale.func(x0 = c(minmaxV,values),minmax = c(1:length(colour.grad)))[-c(1,2)]
  for(i in 1:length(values)){
    diff = abs(c(1:length(colour.grad))-values.scaled[i])
    colo = colour.grad[which(diff == min(diff))]
    draw.circle(x=pTOL[i,"x"],y=pTOL[i,"y"],radius=r,nv=100,border=NULL,col=colo,lty=1,lwd=1)
  }
}

#... Plot the patches with a colour palette
#...........................................
# values: array of patch values to scale to the colour gradient (species richness)
# minmaxV: min and max values to scale to the colour extrema (array of 2 numbers)
# colour.grad: colour gradient
# r: radius of the patches

plot.patch.oneVariable.pies = function(values,colour.palet,r=0.2)
{
  for(i in 1:ncol(values)){
    if(sum(values[,i],na.rm=T)>0)add.pie(x=pTOL[i,"x"],y=pTOL[i,"y"],z=values[,i],labels = NA,radius=r,col = colour.palet,lty=0)
  }
}




#... Plot the patches with patch size and colour gradients
#.........................................................
# value1: array of patch values to scale to the colour gradient
# minmaxV1: min and max values to scale to the colour extrema (array of 2 numbers)
# value2: array of patch values to scale to the pie size
# minmaxV2: min and max values to scale to the pie size extrema (array of 2 numbers)
# colour.grad: colour gradient
# minmaxRay: min and max of pie radius (array of 2 numbers)

plot.patch.twoVariables = function(values1,minmaxV1,values2,minmaxV2,colour.grad,minmaxRay)
{
  values.scaled.col = scale.func(x0 = c(minmaxV1,values1),minmax = c(1:length(colour.grad)))[-c(1,2)]
  values.scaled.ray = scale.func(x0 = sqrt(sqrt(c(minmaxV2,values2))),minmax = minmaxRay)[-c(1,2)]
  for(i in 1:length(values1)){
    diffcol = abs(c(1:length(colour.grad))-values.scaled.col[i])
    colo = colour.grad[which(diffcol == min(diffcol))]
    draw.circle(x=pTOL[i,"x"],y=pTOL[i,"y"],radius=values.scaled.ray[i],nv=100,border=NULL,col=colo,lty=1,lwd=1)
  }
}

#... Plot the patches with patch size and colour gradients
#.........................................................
# value1sp: mat of species x patch percents
# minmaxV1: min and max values to scale to the colour extrema (array of 2 numbers)
# value2: array of patch values to scale to the pie size
# minmaxV2: min and max values to scale to the pie size extrema (array of 2 numbers)
# colour.palet: colours representing each species
# minmaxRay: min and max of pie radius (array of 2 numbers)

plot.patch.twoVariables.pies = function(values1sp,minmaxV1,values2,minmaxV2,colour.palet,minmaxRay)
{
  values.scaled.ray = scale.func(x0 = sqrt(sqrt(c(minmaxV2,values2))),minmax = minmaxRay)[-c(1,2)]
  for(i in 1:ncol(values1sp)){
    if(sum(values[,i],na.rm=T)>0)add.pie(x=pTOL[i,"x"],y=pTOL[i,"y"],z=values1sp[,i],labels = NA,radius=values.scaled.ray[i],col = colour.palet,lty=0)
  }
}


#... Plot the landscape with colour gradient (species richness)
#..............................................................
# L.connectivity: matrix of landscape connectivity
# values: array of patch values to scale to the colour gradient
# legendCol: array of the values to put in colour legend (including min and max)
# colour.gradient: colour gradient
# ray: radius of the patches
# tit: title of the plot
# Test values
#L.connectivity = connect.list[[1]]
#values = runif(36,0.2,5.9); values = rich

plot.landscape.oneVariable = function(L.connectivity,values,legendCol,colour.gradient,ray,tit)
{
  par(mar=c(0,0,0,0))
  plot(NA,xlim=c(0,8),ylim=c(0,8),xlab ="",ylab="", xaxt="n",yaxt="n",frame.plot = F)
  text(3.5,7.5,labels = tit)
  plot.vertices(L.connectivity)
  plot.patch.oneVariable(values,range(legendCol),colour.gradient,r=ray)
  gradient.rect(xleft=7.1,ybottom = 1, xright = 7.5, ytop=6,col = colour.gradient, gradient = "y")
  yticks = scale.func(x0 = legendCol,minmax = c(1,6))
  for(i in 1:length(yticks)){
    segments(7.5,yticks[i],7.7,yticks[i])
    text(7.8,yticks[i],labels=legendCol[i],adj=0,cex=0.8)
  }
  text(7.5,6.5,labels = "S")
}


#... Plot the landscape with colour palette (species percent)
#............................................................
# L.connectivity: matrix of landscape connectivity
# values: matrix of species percent (lines) for each patch (col)
# legendCol_names: anames of species
# colour.palet: array of colour for each species
# ray: radius of the patches
# tit: title of the plot
# Test values
#L.connectivity = connect.list[[1]]
#values = runif(36,0.2,5.9); values = rich

plot.landscape.oneVariable.pies = function(L.connectivity,values,legendCol_names,colour.palet,ray,tit)
{
  par(mar=c(0,0,0,0))
  plot(NA,xlim=c(0,8),ylim=c(0,8),xlab ="",ylab="", xaxt="n",yaxt="n",frame.plot = F)
  text(3.5,7.5,labels = tit)
  plot.vertices(L.connectivity)
  plot.patch.oneVariable.pies(values,colour.palet,r=ray)
  
  yrecCol = scale.func(x0 = 1:length(legendCol_names),minmax = c(3,6))
  for(i in 1:length(yrecCol)){
    x1=7.1; x2=7.4; h=0.15 ; y1=yrecCol[i]+h/2; y2=yrecCol[i]-h/2
    polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=colour.palet[i],lty=0)
    text(7.5,yrecCol[i],labels=legendCol_names[i],adj=0,cex=0.8)
  }
  text(7.4,6.5,labels = "S")
}



#... Plot the landscape with colour gradient (species richness)
#..............................................................
# L.connectivity: matrix of landscape connectivity
# values1: array of patch values to scale to the colour gradient
# legendCol: array of the values to put in colour legend (including min and max)
# colour.gradient: colour gradient
# values2 : array of patch values to scale to the pie radius scale
# legendPie: array of 6 values to put in pie size legend (including min and max)
# minmaxPie: min and max values for the size scale of pies
# tit: title of the plot
# Test values
#L.connectivity = connect.list[[1]]
#values1 = runif(36,0.2,5.9); values2 = runif(36,0,89); values1 = rich ; values2 = ab


plot.landscape.twoVariables = function(L.connectivity,values1,legendCol,colour.gradient,values2, legendPie,minmaxPie,tit)
{
  par(mar=c(0,0,0,0))
  plot(NA,xlim=c(0,9),ylim=c(-1,8),xlab ="",ylab="", xaxt="n",yaxt="n",frame.plot = F)
  text(3.5,7.5,labels = tit)
  plot.vertices(L.connectivity)
  plot.patch.twoVariables(values1,range(legendCol),values2,range(legendPie),colour.gradient,minmaxPie)
  
  gradient.rect(xleft=7.1,ybottom = 1, xright = 7.5, ytop=6,col = colour.gradient, gradient = "y")
  yticksCol = scale.func(x0 = legendCol,minmax = c(1,6))
  for(i in 1:length(yticksCol)){
    segments(7.5,yticksCol[i],7.7,yticksCol[i])
    text(7.8,yticksCol[i],labels=legendCol[i],adj=0,cex=0.8)
  }
  text(7.5,6.5,labels = "S")
  
  centerLegendPie = 1:6
  radiusLegendPie = scale.func(x0 = sqrt(sqrt(legendPie)),minmax = minmaxPie)
  for(i in 1:length(legendPie)){
    draw.circle(x=centerLegendPie[i],y=-0.2,radius=radiusLegendPie[i],nv=100,border=NULL,col="grey",lty=1,lwd=1)
    text(x=centerLegendPie[i],y=-0.8,labels=legendPie[i],cex=0.8)
  }
}

#... Plot the landscape with colour palette (species proportions)
#................................................................
# L.connectivity: matrix of landscape connectivity
# values1sp: matrix of speices percent (lines) for each patch (col)
# legendCol_names: names of species
# colour.palet: colour palette (size of number of species)
# values2 : array of patch values to scale to the pie radius scale
# legendPie: array of 6 values to put in pie size legend (including min and max)
# minmaxPie: min and max values for the size scale of pies
# tit: title of the plot
# Test values
#L.connectivity = connect.list[[1]]
#values1 = runif(36,0.2,5.9); values2 = runif(36,0,89); values1 = rich ; values2 = ab


plot.landscape.twoVariables.pies = function(L.connectivity,values1sp,legendCol_names,colour.palet,values2, legendPie,minmaxPie,tit)
{
  par(mar=c(0,0,0,0))
  plot(NA,xlim=c(0,9),ylim=c(-1,8),xlab ="",ylab="", xaxt="n",yaxt="n",frame.plot = F)
  text(3.5,7.5,labels = tit)
  plot.vertices(L.connectivity)
  plot.patch.twoVariables.pies(values1sp,legendCol_names,values2,range(legendPie),colour.palet,minmaxPie)
  
  yrecCol = scale.func(x0 = 1:length(legendCol_names),minmax = c(3,6))
  for(i in 1:length(yrecCol)){
    x1=7.1; x2=7.4; h=0.15 ; y1=yrecCol[i]+h/2; y2=yrecCol[i]-h/2
    polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=colour.palet[i])
    text(7.5,yrecCol[i],labels=legendCol_names[i],adj=0,cex=0.8)
  }
  text(7.4,6.5,labels = "S")
  
  centerLegendPie = 1:6
  radiusLegendPie = scale.func(x0 = sqrt(sqrt(legendPie)),minmax = minmaxPie)
  for(i in 1:length(legendPie)){
    draw.circle(x=centerLegendPie[i],y=-0.2,radius=radiusLegendPie[i],nv=100,border=NULL,col="grey",lty=1,lwd=1)
    text(x=centerLegendPie[i],y=-0.8,labels=legendPie[i],cex=0.8)
  }
}

