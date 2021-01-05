######################
# This script creates the phylogenetic tree presented in Figure 2a and Supplementary Information 7.
# Date: 11/24/2020
# Author: Francis Denisse McLean-Rodriguez
######################

#################################################################################################################

######################This section will be needed to work on all other scripts
#Load needed packages for the complete workflow followed in the paper
library(adegenet)
library(ade4)
library(ape)
library(apex)
library(car)
library(CMplot)
library(cowplot)
library(ctv)
library(fields)
library(genetics)
library(ggfortify)
library(gghighlight)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggtree)
library(gridExtra)
library(gtools)
library(hierfstat)
library(lattice)
library(LDheatmap)
library(mmod)
library(pegas)
library(phytools)
library(plyr)
library(poppr)
library(qqman)
library(reshape)
library(reshape2)
library(scales)
library(stringr)
library(treeio)
library(vcfR)
library(zoo)

options(stringsAsFactors = F)

#Set the appropriate working directory
wd<-"./output"
setwd(wd)

#Import .vcf file and create genind object
populations.snps <- read.vcf("/input/snp_i_N30.vcf",which.loci = 1:90000)
#Convert the dataframe to a matrix for fast character replacement
populations.snps.mg<-as.matrix(populations.snps)
populations.snps.mg[populations.snps.mg=="./."]<-"NA/NA"
#Produce a genind object
colnames(populations.snps.mg)<-1:ncol(populations.snps.mg)
geninI<-df2genind(populations.snps.mg, sep = "/", NA.char = "NA")
save(geninI,file="geninI.RData")

#Import .csv file with seedling data, bind columns with genind object
insamplesI<-indNames(geninI) #get seedling names 
insamplesI<-cbind(insamplesI, NA) #add dummy column needed for merging
snp_i_N30<-read.csv("./snp_i_N30.csv") #create a dataframe in which we associate each seedling to its metadata
dfinfoI<-merge(insamplesI, snp_i_N30, by.x="insamplesI", by.y="insamplesI") #merge tables
dfinfoI<-dfinfoI[,-2] #eliminate NA column
save(dfinfoI,file="dfinfoI.RData") #save data frame

#################################################################################################################

#Assign populations
geninI<-geninI[order(indNames(geninI)),] #Sort before assigning populations
geninI@pop<-as.factor(dfinfoI$in_or_ex) #label 260 seedlings as "in situ" or "ex situ"

######################Phylogenetic tree
######################For all 260 seedlings
#Neighbor-Joining Tree Estimation
Dind <- dist(geninI@tab,method = "euclidean",diag = FALSE, upper = FALSE, p = 2)
treind <- nj(Dind)

#Define tree tip colors and shapes to distinguish seedlings' in situ or ex situ origin and race, respectively
colsadhoc<-c("#f9d232","#48a0a7")
cols_tree<-rep("#000000",260)
bg_tree<-rep(colsadhoc,each=130)
legend_tree<-matrix(1:26,ncol=2,byrow=T)
lw_tree<-c(rep(1,260))
lgtext_tree<-gsub(c("ORE"),"",levels(as.factor(dfinfoI$pop)))
  lgtext_tree<-gsub("ex","",lgtext_tree)
  lgtext_tree<-gsub("in","",lgtext_tree)
pch_tree<-c(rep(23,20),rep(22,10),rep(21,20),rep(24,10),rep(25,10),rep(24,20),rep(25,20),rep(24,20),
            rep(23,20),rep(22,10),rep(21,20),rep(24,10),rep(25,10),rep(24,20),rep(25,20),rep(24,20)) #selecting shape based on in situ/ex situ origin and race
cex_tree<-c(rep(1.5,30),rep(1.2,10),rep(2.1,10),rep(1.5,80),
            rep(1.5,30),rep(1.2,10),rep(2.1,10),rep(1.5,80))

#Export to .pdf file
pdf("Phylogenetic_tree.pdf",width = 12.5)
plot(treind, type="unrooted",show.tip.lab = F,edge.color = "#696969", edge.width = 0.6, cex=0.5,use.edge.length = T,font=0.5,no.margin=T)
tiplabels(col = cols_tree,  pch=pch_tree, bg=bg_tree, cex=cex_tree,lw=lw_tree)
legend ("left",c("Ex situ", "In situ"),text.font=3, bty="n",xjust = 1, 
        col = rep("#000000",2),
        pch=c(22,22),pt.bg=rep(colsadhoc), pt.lwd = c(1,1), pt.cex=c(3,3),ncol=1,cex = 1.5)
legend ("topleft",c(expression(italic("Ancho")*" - Same seed lot"),
                    expression(italic("Ancho")*" - Different seed lot"),
                    expression(italic("Chalqueño")),
                    expression(italic("Elotes Cónicos")),
                    expression(italic("Cónico")),
                    expression(italic("Pepitilla"))),
        title=expression(bold("Race")),title.adj = 0.05, bty="n",col = "#000000", pch=c(24,25,23,22,21,21),pt.bg="black",
        pt.cex = c(2,2,2,2,1.7,2.6),ncol=1, cex=1.5)
dev.off()

######################For a subsample of 5 out of 10 seedlings per sample for clearer representaton (130 seedlings)
#Create reduced genind object
geninI_xs<-geninI #copy object
drop<-rep(c(1,1,0,1,0,0,0,1,1,0),26) #create list of seedlings to drop
geninI_xs@pop<-as.factor(drop) #drop seedlings
geninI_xs<-seppop(geninI_xs)
geninI_xs<-geninI_xs[[1]]

#Neighbor-Joining Tree Estimation
Dind_xs <- dist(geninI_xs@tab,method = "euclidean",diag = FALSE, upper = FALSE, p = 2) 
treind_xs <- nj(Dind_xs)

#Create reduced parameter lists to graph tree from parameter lists created above
cols_tree_xs<-cols_tree[which(drop==1)]
pch_tree_xs<-pch_tree[which(drop==1)]
bg_tree_xs<-bg_tree[which(drop==1)]
cex_tree_xs<-cex_tree[which(drop==1)]*1.2 #increase tip size
lw_tree_xs<-lw_tree[which(drop==1)]

#Export to .pdf file
pdf("Fig2a_Phylogenetic_tree.pdf",width = 12.5)
plot(treind_xs, type="unrooted",show.tip.lab = F,edge.color = "#696969", edge.width = 0.6, cex=0.5,use.edge.length = T,font=0.5,no.margin=T)
tiplabels(col = cols_tree_xs,  pch=pch_tree_xs, bg=bg_tree_xs, cex=cex_tree_xs,lw=lw_tree_xs)
legend ("left",c("Ex situ", "In situ"),text.font=3, bty="n",xjust = 1, 
        col = rep("#000000",2),
        pch=c(22,22),pt.bg=rep(colsadhoc), pt.lwd = c(1,1), pt.cex=c(3,3),ncol=1,cex = 1.5)
legend ("topleft",c(expression(italic("Ancho")*" - Same seed lot"),
                    expression(italic("Ancho")*" - Different seed lot"),
                    expression(italic("Chalqueño")),
                    expression(italic("Elotes Cónicos")),
                    expression(italic("Cónico")),
                    expression(italic("Pepitilla"))),
        title=expression(bold("Race")),title.adj = 0.05, bty="n",col = "#000000", pch=c(24,25,23,22,21,21),pt.bg="black",
        pt.cex = c(2,2,2,2,1.7,2.6),ncol=1, cex=1.5)
dev.off()
#For clarity, a "Conservation method" label to the ex situ and in situ legend was added manually.

######################Supplementary Information 7
#Define tree tip colors and shapes to distinguish seedlings' in situ or ex situ origin and race, respectively
colsadhoc_SI7<-c("darkred","red","orange","yellow","green","darkgreen","aquamarine4","aquamarine",
                 "darkblue","blue","purple4","darkmagenta","magenta")
cols_tree_SI7<-rep(c(colsadhoc_SI7,rep("#000000",13)),each=5)
bg_tree_SI7<-rep(c(colsadhoc_SI7,colsadhoc_SI7),each=5)
legend_tree_SI7<-matrix(1:26,ncol=2,byrow=T)
lw_tree_SI7<-c(rep(c(2,1),each=65))
lgtext_tree_SI7<-gsub(c("ORE"),"",levels(as.factor(dfinfoI$pop)))
lgtext_tree_SI7<-gsub("ex","",lgtext_tree_SI7)
lgtext_tree_SI7<-gsub("in","",lgtext_tree_SI7)
pch_tree_SI7<-c(rep(5,10),rep(0,5),rep(1,10),rep(2,5),rep(6,5),rep(2,10),rep(6,10),rep(2,10),
                rep(23,10),rep(22,5),rep(21,10),rep(24,5),rep(25,5),rep(24,10),rep(25,10),rep(24,10)) #selecting shape based on in situ/ex situ origin and race
cex_tree_SI7<-c(rep(c(rep(1.5,15),rep(1.2,5),rep(2.1,5),rep(1.5,40)),2))

#Export to .pdf file
pdf("SupplementaryInformation7_Phylotree.pdf",width = 6)
plot(treind_xs, type="unrooted",show.tip.lab = F,edge.color = "#696969", edge.width = 0.6, cex=0.5,use.edge.length = T,font=0.5,no.margin=T)
tiplabels(col = cols_tree_SI7,  pch=pch_tree_SI7, bg=bg_tree_SI7, cex=cex_tree_SI7,lw=lw_tree_SI7)
dev.off()