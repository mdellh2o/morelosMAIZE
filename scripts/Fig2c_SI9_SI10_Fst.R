######################
# This script creates the Fst graph presented in Figure 2c, the table in Supplementary Information 9 and the tree in Supplementary Information 10.
# Date: 11/24/2020
# Author: Francis Denisse McLean-Rodriguez
######################

#load packages, set working directory, import .csv file and create data frame based on first script (Fig2a_SI7_Phylogenetic_tree.R)

#load file created in first script
load(file="dfinfoI.RData")

#A reduced SNP subset excluding SNPs of minor allele frequency <10% is used
#Import .vcf file and create genind object
populations.snps <- read.vcf("/input/snp_i_N30_MAF10.vcf",which.loci = 1:170000)
#convert the dataframe to a matrix for fast character replacement
populations.snps.mg<-as.matrix(populations.snps)
populations.snps.mg[populations.snps.mg=="./."]<-"NA/NA"
#Produce a genind object
colnames(populations.snps.mg)<-1:ncol(populations.snps.mg)
geninI_MAF10<-df2genind(populations.snps.mg, sep = "/", NA.char = "NA")
save(geninI_MAF10,file="geninI_MAF10.RData")

#Assign populations
geninI_MAF10<-geninI_MAF10[indNames(geninI_MAF10),] #Sort before assigning populations
geninI_MAF10@pop<-as.factor(dfinfoI$pop2) c
######################Fst
#Convert genind to hierfstat object
hierfstatI <- genind2hierfstat(geninI_MAF10,pop=geninI_MAF10@pop) 
#Estimating Fst based on Weir and Cockerham (1984)
WC84Fst<-genet.dist(hierfstatI,diploid=TRUE,method="WC84") #long process to run on a personal computer
save(WC84Fst,file="WC84Fst.RData")

#Discard lower triangle and rename populations for graph (for clarity, horizontal "-Ex situ-" and "-In situ-" labels were added maually using an image editor)
y<-as.matrix(WC84Fst)
y[lower.tri(y)] <- NA
colnames(y)<-gsub("IN","In ",gsub("EX","Ex ",gsub("ORE","",colnames(y))))
rownames(y)<-gsub("IN","In ",gsub("EX","Ex ",gsub("ORE","",rownames(y))))
y<-melt(y)

#Export to .pdf file
pdf("Fig2c_Fst.pdf",height = 8,width = 10)
print(ggplot(y, aes(x=Var1, y=Var2, fill=value))+
        geom_tile() +
        geom_vline(xintercept = 13.5) + geom_hline(yintercept = 13.5) +
        ggtitle("Pairwise Fst") + labs(fill = "Fst") +
        scale_fill_gradient(low = "yellow", high = "red",na.value = "white") +
        theme_classic(base_size = 20)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(), axis.text.x = element_text(angle = 90)))
dev.off()

######################Supplementary Information 9
#Get lower triangle from Fst matrix
Fst<-round(as.matrix(WC84Fst),3)
Fst[upper.tri(Fst)] <- NA

#Estimating 95% confidence interval using boostrap
pop2<-rep(1:26, each=10) #create a numeric population vector
bootFst<-boot.ppfst(data.frame(as.numeric(pop2),hierfstatI[,-c(1)],nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=6)) #boostrapping. Long process to run on a personal computer
save(bootFst, file="bootFst.RData")

#Get lower limit
ll<-round(bootFst$ll,3)
ll[lower.tri(ll)] <- t(ll)[lower.tri(ll)]
ll[upper.tri(ll)] <- NA

#Get upper limit
ul<-round(bootFst$ul,3)
ul[lower.tri(ul)] <- t(ul)[lower.tri(ul)]
ul[upper.tri(ul)] <- NA

#Create parentheses and hyphen matrix
m1<-matrix("(",26,26)
m2<-matrix("-",26,26)
m3<-matrix(")",26,26)

#Bind matrixes
SI9 <- matrix(rbind(Fst,m1,ll,m2,ul,m3),nrow =26, byrow = F)
colnames(SI9)<-rep(colnames(Fst),each=6)
rownames(SI9)<-rownames(Fst)
write.csv(SI9,"SupplementaryInformation9.csv") 
#For clarity, column with and labels were edited and empty rows and columns were eliminated in Excel

######################Supplementary Information 10
#From script Fig2a_SI7_Phylogenetic_tree.R
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

#Create Fst-based tree
treind_fst <- nj(WC84Fst)
treind_fst$tip.label<-gsub("ORE","",gsub("IN","\n   ",gsub("EX","\n   ",treind_fst$tip.label)))
treind_fst$tip.label<-str_pad(str_pad(treind_fst$tip.label,8, side = c("right"), pad = " "),9, side = c("right"), pad = " ")
cols_tree_fst<-cols_tree[seq(1, length(cols_tree), 10)]
pch_tree_fst<-pch_tree[seq(1, length(pch_tree), 10)]
bg_tree_fst<-bg_tree[seq(1, length(bg_tree), 10)]
cex_tree_fst<-cex_tree[seq(1, length(cex_tree), 10)]
lw_tree_fst<-lw_tree[seq(1, length(lw_tree), 10)]
rotate

#Export to .pdf file
pdf("SupplementaryInformation10.pdf",width = 12.5)
plot(treind_fst, type="unrooted",show.tip.lab = T, rotate.tree=0, edge.color = "#696969", edge.width = 0.6, cex=0.8,use.edge.length = T,font=0.5,no.margin=T)
tiplabels(col = cols_tree_fst,  pch=pch_tree_fst, bg=bg_tree_fst, cex=cex_tree_fst,lw=lw_tree_fst)
legend (x=0.1, y=0.01,c("Ex situ", "In situ"),title=expression(bold("Conservation method")),text.font=3, bty="n",xjust = 0.5, 
        col = rep("#000000",2),title.adj = 7,
        pch=c(22,22),pt.bg=rep(colsadhoc), pt.lwd = c(1,1), pt.cex=c(3,3),ncol=1,cex = 1.5)
legend ("right",c(expression(italic("Ancho")*" - Same seed lot"),
                  expression(italic("Ancho")*" - Different seed lot"),
                  expression(italic("Chalqueño")),
                  expression(italic("Elotes Cónicos")),
                  expression(italic("Cónico")),
                  expression(italic("Pepitilla"))),
        title=expression(bold("Race")),title.adj = 0.05, bty="n",col = "#000000", pch=c(24,25,23,22,21,21),pt.bg="black",
        pt.cex = c(2,2,2,2,1.7,2.6),ncol=1, cex=1.5)
dev.off()