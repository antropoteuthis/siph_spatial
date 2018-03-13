#Spatial Predator-Prey
rm(list=ls())

library(tidyverse)
library(reshape2)
library(rmutil)
library(vegan)
library(fields)
library(gstat)
library(statnet)
library(spatstat)
library(sp)

zdata = read.csv("~/Dropbox/siphonophore_phylogeny_2017/character_coding/Siphonophore_depth_pruned.tsv", sep='\t')
zdata[,-1] = lapply(zdata[,-1],as.numeric)
rownames(zdata) = zdata$Species

zSiphs = zdata[1:16,]
zPrey = zdata[17:29,]
zPrey$Species = letters[1:length(zPrey$Species)]

Siphs  = list()
for(i in 1:nrow(zSiphs)){
  zi = rnorm(zSiphs$Count[i], zSiphs$Median[i], zSiphs$StdDev[i])
  Siphs[[i]] = zi[which(zi>zSiphs$Min[i], zi<zSiphs$Max[i])]
}
zSiphsX=rep(0,max(zSiphs$Count))
for(i in 1:length(Siphs)){
  #zni = c(Siphs[[i]],rep(NA, (max(zSiphs$Count) - length(Siphs[i]))))
  zni = Siphs[[i]]
  print(i)
  print(length(zni))
  zSiphsX = rbind(zSiphsX,zni)
}
zSiphs_data = melt(cbind(zSiphs$Species,data.frame(zSiphsX[-1,1:ncol(zSiphsX)])))
zSiphs_data = zSiphs_data[which(!(is.na(zSiphs_data$value))),-2]
names(zSiphs_data) = c("Species","Depth")

Prey  = list()
for(i in 1:nrow(zPrey)){
  zi = rnorm(zPrey$Count[i], zPrey$Median[i], zPrey$StdDev[i])
  Prey[[i]] = zi[which(zi>zPrey$Min[i], zi<zPrey$Max[i])]
}
zPreyX=rep(0,max(zPrey$Count))
for(i in 1:length(Prey)){
  #zni = c(Prey[[i]],rep(NA, (max(zPrey$Count) - length(Prey[i]))))
  zni = Prey[[i]]
  print(i)
  print(length(zni))
  zPreyX = rbind(zPreyX,zni)
}
zPrey_data = melt(cbind(zPrey$Species,data.frame(zPreyX[-1,1:ncol(zPreyX)])))
zPrey_data = zPrey_data[which(!(is.na(zPrey_data$value))),-2]
names(zPrey_data) = c("Species","Depth")

###### SHOWTIME ####
distances = as.data.frame(list("SPS", "SPP",0))
colnames(distances) = c("SiphSP", "PreySP", "Zdist")
for(s in 1:nrow(zSiphs_data)){
  zS = zSiphs_data[s,2]
  SPS = zSiphs_data[s,1]
  for(p in 1:nrow(zPrey_data)){
    zP = zPrey_data[p,2]
    SPP = zPrey_data[p,1]
    Zdist = sqrt((zS-zP)^2)
    if(Zdist<1){
      distance_row = as.data.frame(list(SPS, SPP,Zdist))
      colnames(distance_row) = c("SiphSP", "PreySP", "Zdist")
      distances = rbind(distances,distance_row)
      print(nrow(distances))
    }
  }
}
distances = distances[-1,]
Siph_Prey_Matrix = dcast(distances, SiphSP~PreySP, length)
siphnames = Siph_Prey_Matrix[,1]
Siph_Prey_Matrix = as.matrix(sapply(Siph_Prey_Matrix[,-1], as.numeric))
rownames(Siph_Prey_Matrix) = siphnames
Siph_Siph_Matrix = as.matrix(vegdist(Siph_Prey_Matrix, method = "jaccard"))
heatmap(Siph_Prey_Matrix)
heatmap(Siph_Siph_Matrix)
heatmap(as.matrix(vegdist(zSiphs[,-1])))


#### EQUALIZING PREDATORS
Siphs  = list()
for(i in 1:nrow(zSiphs)){
  zi = rnorm(20, zSiphs$Median[i], zSiphs$StdDev[i])
  Siphs[[i]] = zi[which(zi>zSiphs$Min[i], zi<zSiphs$Max[i])]
}
zSiphsX=rep(0,20)
for(i in 1:length(Siphs)){
  #zni = c(Siphs[[i]],rep(NA, (max(zSiphs$Count) - length(Siphs[i]))))
  zni = Siphs[[i]]
  print(i)
  print(length(zni))
  zSiphsX = rbind(zSiphsX,zni)
}
zSiphs_data_eq = melt(cbind(zSiphs$Species,data.frame(zSiphsX[-1,1:ncol(zSiphsX)])))
zSiphs_data_eq = zSiphs_data_eq[,-2]
names(zSiphs_data_eq) = c("Species","Depth")

distances_eq = as.data.frame(list("SPS", "SPP",0))
colnames(distances_eq) = c("SiphSP", "PreySP", "Zdist")
for(s in 1:nrow(zSiphs_data_eq)){
  zS = zSiphs_data_eq[s,2]
  SPS = zSiphs_data_eq[s,1]
  for(p in 1:nrow(zPrey_data)){
    zP = zPrey_data[p,2]
    SPP = zPrey_data[p,1]
    Zdist = sqrt((zS-zP)^2)
    if(Zdist<1){
      distance_row = as.data.frame(list(SPS, SPP,Zdist))
      colnames(distance_row) = c("SiphSP", "PreySP", "Zdist")
      distances_eq = rbind(distances_eq,distance_row)
      print(nrow(distances_eq))
    }
  }
}
distances_eq = distances_eq[-1,]
Siph_Prey_Matrix_EQ = dcast(distances_eq, SiphSP~PreySP, length)
siphnames = Siph_Prey_Matrix_EQ[,1]
Siph_Prey_Matrix_EQ = as.matrix(sapply(Siph_Prey_Matrix_EQ[,-1], as.numeric))
rownames(Siph_Prey_Matrix_EQ) = siphnames
Siph_Siph_Matrix_EQ = as.matrix(vegdist(Siph_Prey_Matrix_EQ, method = "jaccard"))
heatmap(Siph_Prey_Matrix_EQ)
heatmap(Siph_Siph_Matrix_EQ)
mantel(Siph_Siph_Matrix, Siph_Siph_Matrix_EQ)


### Equalizing prey ###
Prey  = list()
for(i in 1:nrow(zPrey)){
  zi = rnorm(mean(zPrey$Count), zPrey$Median[i], zPrey$StdDev[i])
  Prey[[i]] = zi[which(zi>zPrey$Min[i], zi<zPrey$Max[i])]
}
zPreyX=rep(0,mean(zPrey$Count))
for(i in 1:length(Prey)){
  #zni = c(Prey[[i]],rep(NA, (max(zPrey$Count) - length(Prey[i]))))
  zni = Prey[[i]]
  print(i)
  print(length(zni))
  zPreyX = rbind(zPreyX,zni)
}
zPrey_data_eq = melt(cbind(zPrey$Species,data.frame(zPreyX[-1,1:ncol(zPreyX)])))
zPrey_data_eq = zPrey_data_eq[,-2]
names(zPrey_data) = c("Species","Depth")

distances_eq2 = as.data.frame(list("SPS", "SPP",0))
colnames(distances_eq2) = c("SiphSP", "PreySP", "Zdist")
for(s in 1:nrow(zSiphs_data_eq)){
  zS = zSiphs_data_eq[s,2]
  SPS = zSiphs_data_eq[s,1]
  for(p in 1:nrow(zPrey_data_eq)){
    zP = zPrey_data_eq[p,2]
    SPP = zPrey_data_eq[p,1]
    Zdist = sqrt((zS-zP)^2)
    if(Zdist<1){
      distance_row = as.data.frame(list(SPS, SPP,Zdist))
      colnames(distance_row) = c("SiphSP", "PreySP", "Zdist")
      distances_eq2 = rbind(distances_eq2,distance_row)
      print(nrow(distances_eq2))
    }
  }
}
distances_eq2 = distances_eq2[-1,]
Siph_Prey_Matrix_EQ2 = dcast(distances_eq2, SiphSP~PreySP, length)
siphnames = Siph_Prey_Matrix_EQ2[,1]
Siph_Prey_Matrix_EQ2 = as.matrix(sapply(Siph_Prey_Matrix_EQ2[,-1], as.numeric))
rownames(Siph_Prey_Matrix_EQ2) = siphnames
Siph_Siph_Matrix_EQ2 = as.matrix(vegdist(Siph_Prey_Matrix_EQ2, method = "jaccard"))
heatmap(Siph_Prey_Matrix_EQ2)
heatmap(Siph_Siph_Matrix_EQ2)
mantel(Siph_Siph_Matrix_EQ2, Siph_Siph_Matrix)

##Adding XY dimensions
#xyzSiphData  = cbind(zSiphs_data, )

