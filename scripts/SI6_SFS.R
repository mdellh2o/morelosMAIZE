##################
# This script creates the Site Frequency Spectrum graph presented in Supplementary Information 6.
# Date: /11/2020
# Author: Matteo Dell'Acqua
##################

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/MAIZE_morelos/publication/heredity/rev1"
setwd(wd)

#load vcf2sfs function suite
source("https://raw.githubusercontent.com/shenglin-liu/vcf2sfs/master/vcf2sfs.r")
library(vcfR)
library(ggplot2)
library(scales)

#get info file
ind<-read.csv("https://raw.githubusercontent.com/mdellh2o/morelosMAIZE/main/input/snp_i_N30.csv")
#derive a file to be used with vcf2sfs
popmap<-ind[,c("insamplesI", "in_or_ex")]
write.table(popmap, file="popmap.txt", row.names = F, col.names = F, quote=F, sep="\t")

#set pointer to vcf file and popmap file
vcf <- "../../../data/snp_i_N30.vcf"
popfile <- "popmap.txt"

#produce a gt file to handle populations definition
mygt<-vcf2gt(vcf, popfile)

#calculate SFS
exs<-data.frame(gt2sfs.raw(mygt, c("exsitu")))
ins<-data.frame(gt2sfs.raw(mygt, c("insitu")))

mysfs<-gt2sfs.raw(mygt, c("exsitu", "insitu")) 

#write it and plot it
#write.sfs.dadi(mysfs, f.output="trial.txt")
plot.sfs(mysfs)
lines(x = c(0,100), y = c(0,100))

# #make a ggplot about the SFS comparison across populations
# bisfs<-data.frame(mysfs)
# bisfs[which(bisfs$Freq==0),3]<-NA
# ggplot(data=bisfs) + geom_tile(aes(x=exsitu, y=insitu, fill=Freq)) + 
#         scale_fill_gradientn(limits=c(1, 10000), colors=c("red", "yellow", "green", "blue"),  values = scales::rescale(c(1,10,50,100,10871)),na.value = "transparent") +
#         xlab("Ex situ") + ylab("In situ") + 
#         theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_line(colour = "black"))


#make a ggplot comparing the SFS across the two populations
colnames(exs)[1]<-"Class"
exs$Conservation<-"Ex Situ"

colnames(ins)[1]<-"Class"
ins$Conservation<-"In Situ"

toplot<-rbind(exs, ins)
toplot[,1]<-as.factor(toplot[,1])
toplot[,3]<-as.factor(toplot[,3])
head(toplot)

toplot[,1]<-as.numeric(toplot[,1])
toplot<-toplot[!(toplot$Freq==0),] #removing monomorphic markers

#make the plot
sfsplot<-ggplot(toplot, aes(x=Class, y=Freq, fill=Conservation))+
  xlab("Frequency class") + ylab("Frequency") +
  geom_bar(position = "dodge",stat = "identity") + 
  scale_y_log10() +
  scale_fill_manual(values = c("#f9d232","#48a0a7"))+
  theme_bw(base_size = 20) +
  theme(legend.text = element_text(face = "italic"),legend.position="top")+
  labs(fill = "Conservation method")
sfsplot

ggsave(sfsplot, file="SupplementaryInformation6.pdf")

#Kolmogorov-Smirnov test
exsitu<-rowSums(mygt$genotype[,c(1:120,251:260)],na.rm = TRUE)
insitu<-rowSums(mygt$genotype[,121:250],na.rm = TRUE)
ks.test(exsitu,insitu)

#make the plot in a boxplot
ggplot(toplot, aes(y=Freq, fill=Conservation))+ geom_boxplot()+ scale_y_log10()+theme_bw()
totest<-toplot
totest$Freq<-log10(totest$Freq)
t.test(Freq~Conservation, data=totest)
