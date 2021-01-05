######################
# This script creates the PCA-based plot presented in Figure 2b and Supplementary Information 7.
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

######################PCA
#Replace missing values using the accesor tab
x.geninI<-tab(geninI,freq=TRUE,NA.method="mean")
x.geninI<-x.geninI[order(rownames(x.geninI)),]
#Compute the PCA
pcaI<-prcomp(x.geninI)
#get variances
eigI<-(pcaI$sdev)^2
varsI<-eigI*100/sum(eigI)
varsI<-round(varsI,2)

#Define point colors and shapes to distinguish seedling' in situ or ex situ origin and race, respectively (same as in Fig2a_Phylogenetic_tree.pdf)
colsadhoc<-c("#f9d232","#48a0a7")
cols_tree<-rep(colsadhoc,each=130) #points have no outline in this graph
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
pdf("Fig2c_PCA.pdf",width = 8,height = 5)
par(oma=c(1,1,1,10),mar=c(4, 4, 1, 4))
plot(pcaI$x[,1:2], col=transp(cols_tree,0.7), pch=pch_tree, bg=transp(bg_tree,0.7), cex=cex_tree, 
     xlab=paste("PC1 ", varsI[1], "%", sep=""), ylab=paste("PC2 ", varsI[2], "%", sep=""), cex.axis=1.2, cex.lab=1.2)
abline(h=0,col = "black", lty = "solid")
abline(v=0,col = "black", lty = "solid")
par(fig=c(0,1,0,1),oma=c(1,1,1,0),mar=c(0,26,0,1), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", col = "n")
legend ("bottomleft",c("Ex situ", "In situ"),text.font=3, bty="n",xjust = 1, 
        title=expression(bold("Conservation method")),title.adj = 0.05,
        col = rep("#000000",2),pch=c(22,22),pt.bg=rep(colsadhoc), pt.lwd = c(0,0), pt.cex=c(3.5,3.5),ncol=1,cex =1.2)
legend ("topleft",c(expression(italic("Ancho")*" - Same seed lot"),
                      expression(italic("Ancho")*" - Different seed lot"),
                      expression(italic("Chalqueño")),
                      expression(italic("Elotes Cónicos")),
                      expression(italic("Cónico")),
                      expression(italic("Pepitilla"))),
        title=expression(bold("Race")),title.adj = 0.05, bty="n",col = "#000000", pch=c(24,25,23,22,21,21),pt.bg="black",
        pt.cex = c(1.7,1.7,1.7,1.7,1.4,2.3),ncol=1, cex=1.2)
dev.off()

######################Supplementary Information 7
#Define point colors and shapes to distinguish seedling' in situ or ex situ origin and race, respectively (same as in SupplementaryInformation7_Phylotree.pdf)
colsadhoc_SI7_PC<-c("darkred","red","orange","yellow","green","darkgreen","aquamarine4","aquamarine",
                    "darkblue","blue","purple4","darkmagenta","magenta")
cols_tree_SI7_PC<-rep(c(colsadhoc_SI7_PC,rep("#000000",13)),each=10)
bg_tree_SI7_PC<-rep(c(colsadhoc_SI7_PC,colsadhoc_SI7_PC),each=10)
legend_tree_SI7_PC<-matrix(1:26,ncol=2,byrow=T)
lw_tree_SI7_PC<-c(rep(c(2,1),each=130))
lgtext_tree_SI7_PC<-gsub(c("ORE"),"",levels(as.factor(dfinfoI$pop)))
lgtext_tree_SI7_PC<-gsub("ex","",lgtext_tree_SI7_PC)
lgtext_tree_SI7_PC<-gsub("in","",lgtext_tree_SI7_PC)
pch_tree_SI7_PC<-c(rep(5,20),rep(0,10),rep(1,20),rep(2,10),rep(6,10),rep(2,20),rep(6,20),rep(2,20),
                   rep(23,20),rep(22,10),rep(21,20),rep(24,10),rep(25,10),rep(24,20),rep(25,20),rep(24,20)) #selecting shape based on in situ/ex situ origin and race
cex_tree_SI7_PC<-c(rep(1,30),rep(0.9,10),rep(1.4,10),rep(1,80),
                   rep(1,30),rep(0.9,10),rep(1.4,10),rep(1,80))

#Export to .pdf file
pdf("SupplementaryInformation7_PCA.pdf",width = 8,height = 5)
par(oma=c(1,1,1,10),mar=c(4, 4, 1, 4))
plot(pcaI$x[,1:2], col=transp(cols_tree_SI7_PC,0.7), pch=pch_tree_SI7_PC, bg=transp(bg_tree_SI7_PC,0.7), cex=cex_tree_SI7_PC, 
     xlab=paste("PC1 ", varsI[1], "%", sep=""), ylab=paste("PC2 ", varsI[2], "%", sep=""), cex.axis=1, cex.lab=1)
abline(h=0,col = "black", lty = "solid")
abline(v=0,col = "black", lty = "solid")
par(fig=c(0,1,0,1),oma=c(1,1,1,0),mar=c(0,26,0,1), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",col="n")
legend ("topleft",lgtext_tree_SI7_PC[legend_tree_SI7_PC], title= expression(italic("Ex situ   In situ")), bty="n",
        col = c(colsadhoc_SI7_PC,rep("#000000",13)),
        pch=c(5,5,0,1,1,2,6,2,2,6,6,2,2,23,23,22,21,21,24,25,24,24,25,25,24,24),
        pt.bg=c(colsadhoc_SI7_PC,colsadhoc_SI7_PC), pt.lwd = c(rep(2,13),rep(1,13)), pt.cex =c(1,1,1,0.9,1.4,1,1,1,1,1,1,1,1,1,1,1,0.9,1.4,1,1,1,1,1,1,1,1),ncol=2)
legend ("bottomleft",c(expression(italic("Ancho")*" - Same seed lot"),
                       expression(italic("Ancho")*" - Different seed lot"),
                       expression(italic("Chalqueño")),
                       expression(italic("Elotes Cónicos")),
                       expression(italic("Cónico")),
                       expression(italic("Pepitilla"))),
        title=expression(bold("Race")),title.adj = 0.05, bty="n",col = "#000000", pch=c(24,25,23,22,21,21),pt.bg="black",
        pt.cex = c(1,1,1,1,0.9,1.4),ncol=1, cex=1)
dev.off()