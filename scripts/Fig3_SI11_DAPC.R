######################
# This script creates the DAPC-based plots presented in Figure 3, as well as the BIC graph in Supplementary Information 11.
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

######################DAPC
set.seed(20160308)
#Find clusters with k-mean algorithm. 
grp<-find.clusters(geninI,max.n.clust = 40) 
#Choose the number PCs to retain (>=1): 260 (all)
#Choose the number of cluster to retain (>=2): 2 (minimum). 
#Value of BIC versus number of cluster graph is included in manuscript as Supplementary Information 11
grp$size #check how seedlings were divided into the selected k groups

#Identify optimum number of principal components to retain in DAPC
dapctestI<-dapc(geninI,grp$grp,n.pca=260,n.da = 100)
temp<-optim.a.score(dapctestI) 
#Optimal number of PCs retained based on a-score optimisation - spline interpolation = 15
dapcI<-dapc(geninI,grp$grp,n.da=100,n.pca=15) 

#Extract coordinates
x<-as.data.frame(dapcI$ind.coord)
x<-cbind(x,dapcI$grp)
colnames(x)[2]<-"Assigned groups"

#Export to .pdf file
pdf("Figure3a_DAPC.pdf", height=3.3, width=9)
print(ggdensity(x, x = "LD1", rug = TRUE,color = "Assigned groups", fill = "Assigned groups",palette = c("#4c5559","#9ec759"),
                xlab = "Discriminant funcion 1",ylab = "Density",alpha=1) +
        theme(legend.position="none"
              ,axis.title.y = element_text(margin = margin(t = 0, r = 10, b =0, l = 0))
              ,axis.title.x = element_text(margin = margin(t = 10, r = 0, b =0, l = 0))) +
        coord_cartesian(xlim =c(min(x$LD1), max(x$LD1)), ylim = c(0, 1)))
dev.off()

######################Posterior group assignment graph based on script from https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
dapc.results <- as.data.frame(dapcI$posterior) #extract results
dapc.results$pop <- as.factor(dfinfoI$pop) #add colum with population labels
dapc.results$indNames <- rownames(dapc.results) #add columns with seedling labels (individual names)

dapc.results <- melt(dapc.results) #reorganize data frame
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned","Posterior") #rename columns

dapc.results$pop<-gsub("in","",(gsub("ex","",gsub("MORE","M",dapc.results$Original_Pop)))) #simplify population names

dapc.results$pop_f = factor(dapc.results$pop, levels=c('M87','M75','M50','M49','M47','M46','M45','M44','M39','M35','M34','M33','M32')) #change into a factor variable

dapc.resultsex<-dapc.results[c(grep("ex", dapc.results$Original_Pop)),] #subset ex situ seedlings
dapc.resultsin<-dapc.results[c(grep("in", dapc.results$Original_Pop)),] #subset in situ seedlings

#Export to .pdf file
pdf("Figure3b_DAPC_exsitu.pdf", height=2.2, width=9)
print(ggplot(dapc.resultsex, aes(x=Sample, y=Posterior, fill=Assigned))
      + geom_bar(stat='identity') 
      + scale_fill_manual(values = c("#4c5559","#9ec759"))
      + facet_grid(~pop_f, scales = "free")#, switch="x")
      + theme(plot.title = element_text(hjust = 0, face = "bold.italic"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
      +labs(y="", title="Ex situ \n \n", hjust=0))
dev.off()

pdf("Figure3b_DAPC_insitu.pdf", height=1.8, width=9)
print(ggplot(dapc.resultsin, aes(x=Sample, y=Posterior, fill=Assigned))
      + geom_bar(stat='identity') 
      + scale_fill_manual(values = c("#4c5559","#9ec759"))
      + facet_grid(~pop_f, scales = "free", switch="x")
      + theme(plot.title = element_text(hjust = 0, face = "bold.italic"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
      +labs(y="",title="In situ", hjust=0))
dev.off()
#For clarity, location and race labels were added maually using an image editor

######################Variable loadings
#obtain physical positions of all markers
populations.snps2<-read.vcfR("/input/snp_i_N30.vcf")
positionsI<-populations.snps2@fix[,1:3] #extracting positions
positionsI<-cbind(positionsI, c(1:nrow(positionsI)),NA) #add dummy column needed for merging
colnames(positionsI)[4]<-"index"

pos_DAPC<-as.data.frame(positionsI) #adjusting variables' format
pos_DAPC[,1]<-as.numeric(as.character(pos_DAPC[,1]))
pos_DAPC[,2]<-as.numeric(as.character(pos_DAPC[,2]))
pos_DAPC[,4]<-as.numeric(as.character(pos_DAPC[,4]))

y<-as.data.frame(dapcI$var.contr) #extract DAPC variable contribution
y<-cbind(y,rownames(y)) #bind index
y[,2]<-as.numeric(gsub("\\..*", "", y[,2])) #clean index
colnames(y)[2]<-"index"
y<-y[which(duplicated(as.numeric(gsub("\\..*", "", rownames(y)))) %in% y[,2]),] #eliminate duplicate rows

#binding positions w/DAPC loads
pos_DAPC2<-merge(pos_DAPC,y,by.x = "index",by.y = "index")

#add the following two lines to create a vector with the possition of the highest loading SNPs or SNP regions
z<-pos_DAPC2[which(pos_DAPC2[,6] >= 0.0004),c(1:4,6)] #get 5 highest loading SNPs or SNP regions. Different thresholds were tested until 5 SNP regions were obtained
pos_DAPC<-z$ID

#Export to .pdf file
pdf("Figure3c_DAPC.pdf",width = 11,height = 3)
par(mar=c(5,9,2,2))
print(ggplot(pos_DAPC2,aes(x=POS/1000000, y= LD1, color=as.factor(CHROM))) +
        geom_point() +
        #add the line below if the highest loading SNPs or SNP regions will also be highlighted
        geom_point(data=pos_DAPC2[match(pos_DAPC,pos_DAPC2$ID),], aes(x=POS/1000000, y=LD1), colour="yellow", size=0.4)+
        facet_grid(.~CHROM, scales="free_x") +
        scale_color_manual(values = c(rep(c("black","darkgray"),5)))+
        xlab("Genomewide positions (Mb)") + ylab("Variable loadings") +
        theme_classic(base_size = 16) +
        theme(legend.position = "none",axis.text.x.bottom = element_text(angle=60), panel.spacing=unit(0, "lines")) +
        scale_x_continuous(breaks=seq(0,250,100))
)
dev.off()
#For clarity, circles around potential SNPs under selection were added maually using an image editor