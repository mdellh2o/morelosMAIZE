######################
# This script creates the SNP per sample data presented in Supplementary Information 5.
# Date: 11/24/2020
# Author: Francis Denisse McLean-Rodriguez
######################

#load packages, set working directory, import .vcf and .csv files and create genind object and data frame based on first script (Fig2a_SI7_Phylogenetic_tree.R)

#load files created in first script
load(file="geninI.RData")
load(file="dfinfoI.RData")

#Assign populations
geninI<-geninI[order(indNames(geninI)),] #Sort before assigning populations
geninI@pop<-as.factor(dfinfoI$pop) #label each of the 260 seedlings based on their sample number (26 samples 10 seedlings each)

#Separate populations
geninIpop <- seppop(geninI)

#Calculate MAF for each subgroup
#MAF frequencies for each SNP are estimated in two points in time, ex situ and in situ.
for (i in 1:26) { 
  y<- paste("MAFpop",i,sep="")
  assign(y,cbind(as.data.frame(as.numeric(minorAllele(geninIpop[[i]]))),i))
}

#Combine samples into list
MAFpopmat <- lapply(c("MAFpop1","MAFpop2","MAFpop3","MAFpop4","MAFpop5","MAFpop6","MAFpop7","MAFpop8","MAFpop9","MAFpop10","MAFpop11","MAFpop12","MAFpop13","MAFpop14","MAFpop15","MAFpop16","MAFpop17","MAFpop18","MAFpop19","MAFpop20","MAFpop21","MAFpop22","MAFpop23","MAFpop24","MAFpop25","MAFpop26"), get )

#Combine samples into pairs
for (i in 1:13) { 
  x<- paste("MAFpair",i,sep="")
  assign(x,rbind(MAFpopmat[[(i*2)-1]],MAFpopmat[[i*2]]))
}

#Combine pairs of samples into a list
MAFpairmat <- lapply(c("MAFpair1","MAFpair2","MAFpair3","MAFpair4","MAFpair5","MAFpair6","MAFpair7","MAFpair8","MAFpair9","MAFpair10","MAFpair11","MAFpair12","MAFpair13"), get )


#obtain physical positions of all markers
populations.snps2<-read.vcfR("/input/snp_i_N30.vcf")
positionsI<-populations.snps2@fix[,1:3] #extracting positions
positionsI<-cbind(positionsI, c(1:nrow(positionsI)),NA) #add dummy column needed for merging
colnames(positionsI)[4]<-"index"

pos_DAPC<-as.data.frame(positionsI) #adjusting variables' format
pos_DAPC[,1]<-as.numeric(as.character(pos_DAPC[,1]))
pos_DAPC[,2]<-as.numeric(as.character(pos_DAPC[,2]))
pos_DAPC[,4]<-as.numeric(as.character(pos_DAPC[,4]))

#Obtain SNP positions
pos<-as.data.frame(positionsI[,c(1,2,4,5)])
for (i in 1:4){
  pos[,i]<-as.numeric(as.character(pos[,i]))
}

#Bind SNP positions with MAF in samples (MAFpopmat) and pairs (MAFpairmat)
for (i in 1:26) { 
  MAFpopmat[[i]]<-cbind(pos,MAFpopmat[[i]])
}

#Duplicate to bind with sample pairs
pos<-rbind(pos,pos) 

#Bind SNP positions with MAF
for (i in 1:13) { 
  MAFpairmat[[i]]<-cbind(pos,MAFpairmat[[i]])
}

#Separate SNPS from monomorphic sites (snps with a MAF of 1 or 0)
pairs<-c("M32","M33","M34","M35","M39","M44","M45","M46","M47","M49","M50","M75","M87")
MAFpairSNP<-list()
MAFpairNSNP<-list()
for (j in 1:13){
  colnames(MAFpairmat[[j]])[5]<-"MAF" #rename column
  #relabel in situ and ex situ columns
  MAFpairmat[[j]]$i[MAFpairmat[[j]]$i== (j*2)-1]<-"Ex-situ" 
  MAFpairmat[[j]]$i[MAFpairmat[[j]]$i== j*2 ]<-"In-situ"
  MAFpairmat[[j]]$i<-as.factor(MAFpairmat[[j]]$i)
  #separate in two lists 
  nam<-paste(pairs[j],"SNP",sep="")
  assign(nam, MAFpairmat[[j]][-c(which(MAFpairmat[[j]]$MAF %in% c(0,1,NA))),])
  MAFpairSNP[[j]]<-get(nam)
  nam2<-paste(pairs[j],"NSNP",sep="")
  assign(nam2, MAFpairmat[[j]][c(which(MAFpairmat[[j]]$MAF %in% c(0,1,NA))),])
  MAFpairNSNP[[j]]<-get(nam2)
}

#Chi-squared
p.val<-c()

for (j in 1:13){ 
  MAFpairSNP[[j]]$MAF_cat[MAFpairSNP[[j]]$MAF < 0.01]<- "1. Rare < 1%"
  MAFpairSNP[[j]]$MAF_cat[MAFpairSNP[[j]]$MAF>=0.01 & MAFpairSNP[[j]]$MAF <= 0.05]<- "2. Low frequency 1-5%"
  MAFpairSNP[[j]]$MAF_cat[MAFpairSNP[[j]]$MAF>0.05]<- "3. Common > 5%"
  MAFpairSNP[[j]]$MAF_cat<-as.factor(MAFpairSNP[[j]]$MAF_cat)
  
  temp<-chisq.test(MAFpairSNP[[j]]$i,MAFpairSNP[[j]]$MAF_cat)
  p.val[j]<-temp$p.value
}

SNP_freq<-list()
for (j in 1:13){ 
  nam<-paste("SNP_",pairs[j],sep="")
  assign(nam,rbind(t(table(MAFpairNSNP[[j]]$i)),t(table(MAFpairSNP[[j]]$i,MAFpairSNP[[j]]$MAF_cat))))
  print(colSums(get(nam)))
  SNP_freq[[j]]<-get(nam)
}

#SNPs by samples in 13 pairs
SNP_13pairs<-list()
chisq_13pairs<-list()
for (j in 1:13){
  SNP_13pairs[[j]]<-rbind(SNP_freq[[j]][1,],colSums(SNP_freq[[j]][2:3,]))
  chisq_13pairs[[j]]<-chisq.test(SNP_13pairs[[j]])
}

SNP_13pairs_sum<-list()
for (j in 1:13){
  SNP_13pairs_sum[[j]]<-as.data.frame(t(c(pairs[j],SNP_13pairs[[j]][2,1],SNP_13pairs[[j]][2,2],chisq_13pairs[[j]]$p.value)))
}
SNP_13pairs_sum<-do.call(rbind,SNP_13pairs_sum)
SNP_13pairs_sum

#For clarity headers were edited in Excel