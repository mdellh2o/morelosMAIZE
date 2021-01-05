######################
# This script creates the heterozygosity per sample data presented in Supplementary Information 5.
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

#Create summaries for populations
summ26 <- do.call("c", lapply(geninIpop, function(x) summary(x))) #Long process to run on a personal computer

#BARLETT TESTS
accessions<-c(1:26)
bt26<-data.frame(accessions)
#Observed vs expected heterozygosity in each sample
for (i in 1:26) { 
  assign("temp",bartlett.test(list(summ26[[(i*7)-1]],summ26[[i*7]])))
  bt26$p_val_bt[i]<-get("temp")$p.value
  bt26$var[bt26$p_val_bt>=0.05] <- "homogeneus"
  bt26$var[bt26$p_val_bt<0.05] <- "heterogeneus"
}

#Observed heterozygosity between each ex situ and in situ pair
pairs<-c(1:13)
bt13Hopairs<-data.frame(pairs)
summ26Ho<-summ26[c(grep("Ho",names(summ26)))]
for (i in 1:13) { 
  assign("temp",bartlett.test(list(summ26Ho[[(i*2)-1]],summ26Ho[[i*2]])))
  bt13Hopairs$p_val_bt[i]<-get("temp")$p.value
  bt13Hopairs$var[bt13Hopairs$p_val_bt>=0.05] <- "homogeneus"
  bt13Hopairs$var[bt13Hopairs$p_val_bt<0.05] <- "heterogeneus"
}


#Expected heterozygosity between each ex situ and in situ pair
pairs<-c(1:13)
bt13Hepairs<-data.frame(pairs)
summ26He<-summ26[c(grep("He",names(summ26)))]
for (i in 1:13) { 
  assign("temp",bartlett.test(list(summ26He[[(i*2)-1]],summ26He[[i*2]])))
  bt13Hepairs$p_val_bt[i]<-get("temp")$p.value
  bt13Hepairs$var[bt13Hepairs$p_val_bt>=0.05] <- "homogeneus"
  bt13Hepairs$var[bt13Hepairs$p_val_bt<0.05] <- "heterogeneus"
}


#T TESTS
#Observed vs expected heterozygosity in each sample
accessions<-c(1:26)
tt26<-data.frame(accessions)
for (i in 1:26) { 
  if (bt26[i,3]=="homogeneus") {
    assign("temp",t.test(summ26[[(i*7)-1]],summ26[[i*7]],pair=T,var.equal=TRUE))
  } else {
    assign("temp",t.test(summ26[[(i*7)-1]],summ26[[i*7]],pair=T,var.equal=FALSE))
  }
  tt26$p_val_tt[i]<-get("temp")$p.value
  tt26$sig[tt26$p_val_tt>=0.05] <- "non-significant"
  tt26$sig[tt26$p_val_tt<0.05] <- "significant"
  tt26$mean_diff[i]<-get("temp")$estimate
}

#Observed heterozygosity between each ex situ and in situ pair
pairs<-c(1:13)
tt13Hopairs<-data.frame(pairs)
for (i in 1:13) { 
  if (bt13Hopairs[i,3]=="homogeneus") {
    assign("temp",t.test(summ26Ho[[(i*2)-1]],summ26Ho[[i*2]],pair=F,var.equal=TRUE)) 
  } else {
    assign("temp",t.test(summ26Ho[[(i*2)-1]],summ26Ho[[i*2]],pair=F,var.equal=FALSE)) 
  }
  tt13Hopairs$p_val_tt[i]<-get("temp")$p.value
  tt13Hopairs$sig[tt13Hopairs$p_val_tt>=0.05] <- "non-significant"
  tt13Hopairs$sig[tt13Hopairs$p_val_tt<0.05] <- "significant"
  tt13Hopairs$mean_m.x[i]<-get("temp")$estimate[1]
  tt13Hopairs$mean_m.y[i]<-get("temp")$estimate[2]
}

#Expected heterozygosity between each ex situ and in situ pair
pairs<-c(1:13)
tt13Hepairs<-data.frame(pairs)
for (i in 1:13) { 
  if (bt13Hepairs[i,3]=="homogeneus") {
    assign("temp",t.test(summ26He[[(i*2)-1]],summ26He[[i*2]],pair=F,var.equal=TRUE)) 
  } else {
    assign("temp",t.test(summ26He[[(i*2)-1]],summ26He[[i*2]],pair=F,var.equal=FALSE)) 
  }
  tt13Hepairs$p_val_tt[i]<-get("temp")$p.value
  tt13Hepairs$sig[tt13Hepairs$p_val_tt>=0.05] <- "non-significant"
  tt13Hepairs$sig[tt13Hepairs$p_val_tt<0.05] <- "significant"
  tt13Hepairs$mean_m.x[i]<-get("temp")$estimate[1]
  tt13Hepairs$mean_m.y[i]<-get("temp")$estimate[2]
}

bt26
bt13Hepairs
bt13Hopairs

tt26
tt13Hepairs
tt13Hopairs

#For clarity values were rearranged manually in Excel