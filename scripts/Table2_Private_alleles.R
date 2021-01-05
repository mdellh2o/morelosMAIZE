######################
# This script creates the private allele data presented in Table 2.
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

#importing minor and major allele frequencies
#MAF frequencies take for all ex situ and all in situ accessions separately are used, as we need the global MAF, 
#but it is the same allele sampled in two points in time. This is the same as in GENPaper_MAF through different methods
GenoSummaryAlleles_ex<-read.csv("/input/GenoSummaryAlleles_ex.csv")
GenoSummaryAlleles_in<-read.csv("/input/GenoSummaryAlleles_in.csv")

#subsetting a reduced allele frequencies dataset
allfreq_ex <- GenoSummaryAlleles_ex[,c("Site.Name","Major.Allele","Major.Allele.Frequency","Minor.Allele","Minor.Allele.Frequency")]
allfreq_in <- GenoSummaryAlleles_in[,c("Site.Name","Major.Allele","Major.Allele.Frequency","Minor.Allele","Minor.Allele.Frequency")]


#finding private alleles

  PA_FALSE<-private_alleles(geninI, form = alleles ~ ., report = "table", level = "population", count.alleles = FALSE, drop = TRUE)
  #obtain index of private alleles
  temp<-gsub(".*\\.","", colnames(PA_FALSE)) #returns the alleles in colname vector
  colnames(PA_FALSE)<-as.numeric(gsub("\\..*", "", colnames(PA_FALSE))) #returns the index on colname vector
  PA_FALSE<-as.data.frame(t(PA_FALSE))
  PA_FALSE<-PA_FALSE[c(which((rowSums(PA_FALSE))>0)),]
  PA_FALSE<-cbind(PA_FALSE, rownames(PA_FALSE),temp) #binding allele and index column
  colnames(PA_FALSE)[3:4]<-c("index","allele") 
  PA_FALSE$pa_freq<-0 #empty column for private allele frequency
  PA_FALSE$pa_freq_cat<-0 #empty column for private allele frequency category

  x<-c("ex","in")
  
  for (j in 1:2){
  
  temp <- PA_FALSE[PA_FALSE[,j]==1,]
  temp2<- get(paste("allfreq_",x[j],sep = ""))
  
  for (i in 1:nrow(temp)){ #filling private allele frequency column based on which allele was private at each loci
    
    if (temp$allele[i]==temp2$Major.Allele[as.numeric(temp$index[i])])
    {temp$pa_freq[i] <- temp2$Major.Allele.Frequency[as.numeric(temp$index[i])]}
    
    if (temp$allele[i]==temp2$Minor.Allele[as.numeric(temp$index[i])])
    {temp$pa_freq[i] <- temp2$Minor.Allele.Frequency[as.numeric(temp$index[i])]}
  }
  temp$pa_freq_cat[temp$pa_freq < 0.01] <- "1. Rare < 1%"
  temp$pa_freq_cat[temp$pa_freq >= 0.01 & temp$pa_freq <= 0.05]<- "2. Low frequency 1-5%"
  temp$pa_freq_cat[temp$pa_freq > 0.05]<- "3. Common > 5%"
  temp$pa_freq_cat <- as.factor(temp$pa_freq_cat)
  assign(paste("PA_",x[j],sep=""),temp)
  }

addmargins(rbind(table(PA_ex$pa_freq_cat),table(PA_in$pa_freq_cat)))
addmargins(rbind(table(PA_ex$pa_freq_cat),table(PA_in$pa_freq_cat))) / 74739
print(chisq.test(rbind(table(PA_ex$pa_freq_cat),table(PA_in$pa_freq_cat)))) #chi square of private allele proportion by frequency category

#chi square test of ex situ vs in situ PA proportion compared to common alleles
geninex<-popsub(geninI,sublist = c("exsitu"))
geninin<-popsub(geninI,sublist = c("insitu"))#to find number of alleles per group
dim(geninex$tab)[2] #Alleles ex situ
dim(geninin$tab)[2] #Alleles in situ
temp3<-as.data.frame(rbind(cbind(nrow(PA_ex),nrow(PA_in)),c(dim(geninex$tab)[2],dim(geninin$tab)[2]))) #binding ex situ and in situ number of PA along with total number of alleles
temp3[3,]<-c(temp3[2,]-temp3[1,]) #estimating number of common alleles
temp3<-temp3[-2,] #dropping row with total number of alleles
print(chisq.test(temp3)) #test

#For clarity values were rearranged manually in Excel