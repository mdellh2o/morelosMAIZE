######################
# This script creates the SNPs by MAF frequency data presented in Table 2.
# Date: 11/24/2020
# Author: Francis Denisse McLean-Rodriguez
######################

#load packages, set working directory, import .vcf and .csv files and create genind object and data frame based on first script (Fig2a_SI7_Phylogenetic_tree.R)

#load files created in first script
load(file="geninI.RData")
load(file="dfinfoI.RData")

#Assign populations
geninI<-geninI[order(indNames(geninI)),] #Sort before assigning populations
geninI@pop<-as.factor(dfinfoI$in_or_ex) #label each of the 260 seedlings as "in situ" or "ex situ"

#Separate populations
grp_inex<-seppop(geninI)

#Calculate MAF for each subgroup
#MAF frequencies for each SNP are estimated in two points in time, ex situ and in situ.
for (i in 1:2) { 
  y<- paste("MAF_inex",i,sep="")
  assign(y,cbind(as.data.frame(as.numeric(minorAllele(grp_inex[[i]]))),i))
}

#Bind the two groups to compare
MAF_inex<-rbind(MAF_inex1,MAF_inex2)

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

#Duplicate to bind with population pairs
pos<-rbind(pos,pos) 

#Bind SNP positions with MAF
MAFgrpmat<-cbind(pos,MAF_inex)

###############SNPs Chi-square
  colnames(MAFgrpmat)[5]<-"MAF"
  colnames(MAFgrpmat)[6]<-"i"
  MAFgrpmat$i[MAFgrpmat$i== 1]<-"Ex situ"
  MAFgrpmat$i[MAFgrpmat$i== 2]<-"In situ"
  MAFgrpmat$i<-as.factor(MAFgrpmat$i)
  MAFgrpmat$SNP <-"SNP"
  MAFgrpmat$SNP[MAFgrpmat$MAF == 0 | MAFgrpmat$MAF == 1 | MAFgrpmat$MAF == "NA" ]<- "NSNP"
  MAFgrpmat$SNP<-as.factor(MAFgrpmat$SNP)

  addmargins(table(MAFgrpmat$i,MAFgrpmat$SNP))
  prop.table(table(MAFgrpmat$i,MAFgrpmat$SNP),1)
  print(chisq.test(MAFgrpmat$i,MAFgrpmat$SNP))
  
################MAF Chi-square

#Separate SNPS from monomorphic sites (snps with a MAF of 1 or 0) (labels above)
MAFgrpSNP<-MAFgrpmat[-c(which(MAFgrpmat$MAF %in% c(0,1,NA))),]
MAFgrpNSNP<-MAFgrpmat[c(which(MAFgrpmat$MAF %in% c(0,1,NA))),]

#Chi-squared
  breaks <- seq(0, 0.5, by=0.05)
  breaks <- c(0,0.01,0.05,0.5)
  MAFgrpSNP$MAF_cat4 <- cut(MAFgrpSNP$MAF, breaks=breaks, include.lowest=TRUE, right=FALSE)
  MAFgrpSNP$MAF_cat4<-as.factor(MAFgrpSNP$MAF_cat4)
  chisq.test(MAFgrpSNP$MAF_cat4,MAFgrpSNP$i) 
  addmargins(table(MAFgrpSNP$MAF_cat4,MAFgrpSNP$i))
  addmargins(table(MAFgrpSNP$MAF_cat4,MAFgrpSNP$i))/74739
  
#For clarity values were rearranged manually in Excel