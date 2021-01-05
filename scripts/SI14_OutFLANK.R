##################
# This script implements the outlier identification method with OutFLANK presented in Supplementary Information 14 and as part of Table 4 (only SNPs overlapping with BayeScan analysis).
# OutFLANK is distributed through GitHub (https://github.com/whitlock/OutFLANK)
# Date: /11/2020
# Author: Matteo Dell'Acqua
##################

#add missing function
getFSTs_diploids = function(popNameList, SNPDataColumn){  
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]
  
  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0
  
  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }
  
  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)} 
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }
  
  out = WC_FST_Diploids_2Alleles(HetCounts)	
  return(out)
}

library(devtools)
library(qvalue)
library(vcfR)
#devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)

#get vcf
vcf <- read.vcfR("https://raw.githubusercontent.com/mdellh2o/morelosMAIZE/main/input/snp_i_N30_MAF1.vcf", verbose=FALSE)

#get name of loci
locinames <- paste(vcf@fix[,"CHROM"], vcf@fix[,"POS"], sep="_")
head(locinames)

#convert vcf to file handled by OutFLANK
gen_table <- as.matrix(vcf@gt[,-1])
gen_table[gen_table =="0/0"] <- 0
gen_table[gen_table =="0/1"] <- 1
gen_table[gen_table =="1/0"] <- 1
gen_table[gen_table =="1/1"] <- 2
gen_table[is.na(gen_table)]<-9
gen_table<-data.frame(gen_table)
gen_table<-apply(gen_table,2,as.numeric)

colnames(gen_table)<-colnames(vcf@gt[,-1])
rownames(gen_table)<-locinames
gen_table[1:4,1:4]

#clean the genotyping from bad markers
todrop<-apply(gen_table,1,sum)


#get SNP data for outFLANK
SNPdata <- data.frame(t(gen_table))

#get info file
ind<-read.csv("https://raw.githubusercontent.com/mdellh2o/morelosMAIZE/main/input/snp_i_N30.csv")

#sort datasets
ind<-ind[order(ind[,"insamplesI"]),]
SNPdata<-SNPdata[order(rownames(SNPdata)),]
stopifnot(all(rownames(SNPdata) == ind[,"insamplesI"]))

#set the criterion to define populations (setting column name of ind)
idx<-"in_or_ex"

#get data for outFLANK
FstDataFrame <- MakeDiploidFSTMat(SNPdata,locinames,ind[,idx])
head(FstDataFrame)

#check data for markers with unusual Fst values
plot(FstDataFrame$FST, FstDataFrame$FSTNoCorr, 
     xlim=c(-0.01,0.3), ylim=c(-0.01,0.3),
     pch=20)
abline(0,1)

# Note the large FST values for loci with low heterozygosity (He < 0.1)
hist(FstDataFrame$FSTNoCor, breaks=300)

# Run outflank
outlier <- OutFLANK(FstDataFrame, NumberOfSamples=length(unique(ind[,idx])))

# Visualize results
OutFLANKResultsPlotter(outlier, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

sum(outlier$results$qvalues<0.01, na.rm=TRUE)


### Note how OutFLANK identifies potential outliers at He < 0.1, even though
### these loci were excluded in the trimming algorithm

# Create logical vector for top candidates
top_candidates <- outlier$results$qvalues<0.05 & outlier$results$He>0.1

plot(outlier$results$He, outlier$results$FST, pch=20, col="grey")
points(outlier$results$He[top_candidates], outlier$results$FST[top_candidates], pch=21, col="blue")
topcan <- outlier$results[top_candidates,]
topcan[order(topcan$LocusName),]

write.table(topcan, file="SupplementaryInformation14.csv", row.names=F, quote=F)


# clean up by removing low freq variants
keep <- outlier$results$He>0.1 & !is.na(outlier$results$He)
plot(vcf@fix[keep,"POS"], outlier$results$FST[keep], 
     col=grey((as.numeric(vcf@fix[keep,"CHROM"])%%2+1)/3),
     pch=20
) 
points(vcf@fix[top_candidates,"POS"], outlier$results$FST[top_candidates], 
       pch=21, cex=2, col=2)


