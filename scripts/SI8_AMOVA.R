######################
# This script creates the AMOVA tests presented in Supplementary Information 8.
# Date: 11/24/2020
# Author: Francis Denisse McLean-Rodriguez
######################

#load packages, set working directory, import .csv file and create data frame based on first script (Fig2a_SI7_Phylogenetic_tree.R)
#import .vcf file and create genind object based on script Fig2c_SI9_SI10_Fst.R

#load files created in previous script
load(file="geninI_MAF10.RData")
load(file="dfinfoI.RData")

#AMOVA using seedlings' samples as grouping variables (26 groups) 
geninI_MAF10<-geninI_MAF10[order(indNames(geninI_MAF10)),]
geninI_MAF10@pop<-as.factor(dfinfoI$pop) #assign SAMPLES populations

strata(geninI_MAF10)<-(data.frame(pop(geninI_MAF10)))
nameStrata(geninI_MAF10)<-~accessions
head(strata(geninI_MAF10))

AMOVA_ACC<-poppr.amova(geninI_MAF10,hier=~accessions, missing="loci",cutoff=0.3, method = "ade4")

#AMOVA using seedlings' sample pairs as grouping variables (13 groups) 
geninI_MAF10<-geninI_MAF10[order(indNames(geninI_MAF10)),]
geninI_MAF10@pop<-as.factor(dfinfoI$popMORE) #assign PAIRS populations

strata(geninI_MAF10)<-(data.frame(pop(geninI_MAF10)))
nameStrata(geninI_MAF10)<-~pairs
head(strata(geninI_MAF10))

AMOVA_PAIRS<-poppr.amova(geninI_MAF10,hier=~pairs, missing="loci",cutoff=0.3, method = "ade4")

#AMOVA using seedlings' in situ or ex situ origins as grouping variables (2 groups) 
geninI_MAF10<-geninI_MAF10[order(indNames(geninI_MAF10)),]
geninI_MAF10@pop<-as.factor(dfinfoI$in_or_ex) #assign IN SITU/EX SITU populations

strata(geninI_MAF10)<-(data.frame(pop(geninI_MAF10)))
nameStrata(geninI_MAF10)<-~inex
head(strata(geninI_MAF10))

AMOVA_IN_EX<-poppr.amova(geninI_MAF10,hier=~inex, missing="loci",cutoff=0.3, method = "ade4")

#AMOVA using seedlings' races as grouping variables (5 groups) 
geninI_MAF10<-geninI_MAF10[order(indNames(geninI_MAF10)),]
geninI_MAF10@pop<-as.factor(dfinfoI$race_1st) #assign PRIMARY RACE populations

strata(geninI_MAF10)<-(data.frame(pop(geninI_MAF10)))
nameStrata(geninI_MAF10)<-~p_race
head(strata(geninI_MAF10))

AMOVA_1RACE<-poppr.amova(geninI_MAF10,hier=~p_race, missing="loci",cutoff=0.3, method = "ade4")

#AMOVA using seedlings' municipalities as grouping variables (4 groups) 
geninI_MAF10<-geninI_MAF10[order(indNames(geninI_MAF10)),]
geninI_MAF10@pop<-as.factor(dfinfoI$municipality) #assign MUNICIPALITY populations

strata(geninI_MAF10)<-(data.frame(pop(geninI_MAF10)))
nameStrata(geninI_MAF10)<-~municipality
head(strata(geninI_MAF10))

AMOVA_MUN<-poppr.amova(geninI_MAF10,hier=~municipality, missing="loci",cutoff=0.3, method = "ade4")

AMOVAtable<-cbind(AMOVA_ACC$componentsofcovariance,AMOVA_PAIRS$componentsofcovariance,AMOVA_IN_EX$componentsofcovariance,
             AMOVA_1RACE$componentsofcovariance,AMOVA_MUN$componentsofcovariance)

write.csv(AMOVAtable,"SupplementaryInformation8.csv") 
#For clarity, headers were edited in Excel