######################
# This script creates the heterozygosity data presented in Table 2.
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

#Create summaries for populations
summ_inex <- do.call("c", lapply(grp_inex, function(x) summary(x))) #Long process to run on a personal computer

#######################start from here ALL SAMPLES
#BARLETT TESTS
inex<-c(1:2)
bt_inex<-data.frame(inex)
#Observed vs expected heterozygosity
for (i in 1:2) { 
  assign("temp",bartlett.test(list(summ_inex[[(i*7)-1]],summ_inex[[i*7]])))
  bt_inex$p_val_bt[i]<-get("temp")$p.value
  bt_inex$var[bt_inex$p_val_bt>=0.05] <- "homogeneus"
  bt_inex$var[bt_inex$p_val_bt<0.05] <- "heterogeneus"
}

#Observed heterozygosity in situ vs ex situ
inex<-1
bt_inexHo<-data.frame(inex)

summ_inexHo<-summ_inex[c(grep("Ho",names(summ_inex)))]
  
temp<-bartlett.test(list(summ_inexHo[[1]],summ_inexHo[[2]]))
bt_inexHo$p_val_bt[1]<-get("temp")$p.value
bt_inexHo$var[bt_inexHo$p_val_bt>=0.05] <- "homogeneus"
bt_inexHo$var[bt_inexHo$p_val_bt<0.05] <- "heterogeneus"

#Expected heterozygosity in situ vs ex situ
inex<-1
bt_inexHe<-data.frame(inex)

summ_inexHe<-summ_inex[c(grep("He",names(summ_inex)))]

temp<-bartlett.test(list(summ_inexHe[[1]],summ_inexHe[[2]]))
bt_inexHe$p_val_bt[1]<-get("temp")$p.value
bt_inexHe$var[bt_inexHe$p_val_bt>=0.05] <- "homogeneus"
bt_inexHe$var[bt_inexHe$p_val_bt<0.05] <- "heterogeneus"


#T TESTS
#Observed vs expected heterozygosity
inex<-c(1:2)
tt_inex<-data.frame(inex)
for (i in 1) { 
  if (bt_inex[i,3]=="homogeneus") {
    assign("temp",t.test(summ_inex[[(i*7)-1]],summ_inex[[i*7]],pair=T,var.equal=TRUE))
  } else {
    assign("temp",t.test(summ_inex[[(i*7)-1]],summ_inex[[i*7]],pair=T,var.equal=FALSE))
  }
  tt_inex$p_val_tt[i]<-get("temp")$p.value
  tt_inex$sig[tt_inex$p_val_tt>=0.05] <- "non-significant"
  tt_inex$sig[tt_inex$p_val_tt<0.05] <- "significant"
  tt_inex$mean_diff[i]<-get("temp")$estimate}

#Observed heterozygosity in situ vs ex situ
inexp<-1
tt_inexHo<-data.frame(inexp)

  if (bt_inexHo[1,3]=="homogeneus") {
    temp<-t.test(summ_inexHo[[1]],summ_inexHo[[2]],pair=F,var.equal=TRUE)
  } else {
    temp<-t.test(summ_inexHo[[1]],summ_inexHo[[2]],pair=F,var.equal=FALSE)
  }
  tt_inexHo$p_val_tt[1]<-get("temp")$p.value
  tt_inexHo$sig[tt_inexHo$p_val_tt>=0.05] <- "non-significant"
  tt_inexHo$sig[tt_inexHo$p_val_tt<0.05] <- "significant"
  tt_inexHo$mean_m.x[1]<-get("temp")$estimate[1]
  tt_inexHo$mean_m.y[1]<-get("temp")$estimate[2]

#Expected heterozygosity in situ vs ex situ
  inexp<-1
  tt_inexHe<-data.frame(inexp)
  
  if (bt_inexHe[1,3]=="homogeneus") {
    temp<-t.test(summ_inexHe[[1]],summ_inexHe[[2]],pair=F,var.equal=TRUE)
  } else {
    temp<-t.test(summ_inexHe[[1]],summ_inexHe[[2]],pair=F,var.equal=FALSE)
  }
  tt_inexHe$p_val_tt[1]<-get("temp")$p.value
  tt_inexHe$sig[tt_inexHe$p_val_tt>=0.05] <- "non-significant"
  tt_inexHe$sig[tt_inexHe$p_val_tt<0.05] <- "significant"
  tt_inexHe$mean_m.x[1]<-get("temp")$estimate[1]
  tt_inexHe$mean_m.y[1]<-get("temp")$estimate[2]

bt_inex
bt_inexHe
bt_inexHo

tt_inex
tt_inexHe
tt_inexHo

#For clarity values were rearranged manually in Excel