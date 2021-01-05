######################
# This script creates the SNP density graph by chromosome presented in Supplementary Information 4.
# Date: 11/24/2020
# Author: Francis Denisse McLean-Rodriguez
######################

#load packages and set working directory based on first script (Fig2a_SI7_Phylogenetic_tree.R)

#Import .vcf file
populations.snps2<-read.vcfR("/input/snp_i_N30.vcf")
#Extract positions
pos<-as.data.frame(populations.snps2@fix[,1:3])
pos<-pos[,c(3,1,2)]
pos[,1]<-as.character(pos[,1])
pos[[2]]<-as.numeric(as.character(pos[[2]]))
pos[[3]]<-as.numeric(as.character(pos[[3]]))
#Estimate SNP density accross all genome
CMplot(pos,plot.type="d",bin.size=1e6,chr.den.col=c("yellow","orange","red"),file="jpg",memo="",dpi=300)

