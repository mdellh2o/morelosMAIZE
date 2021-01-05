######################
# This script creates the LD plots presented in Figure 4, as well as the table in Supplementary Information 12.
# Date: 11/24/2020
# Author: Francis Denisse McLean-Rodriguez
######################

#load packages, set working directory, import .vcf and .csv files and create genind object and data frame based on first script (Fig2a_SI7_Phylogenetic_tree.R)

#A reduced SNP subset excluding SNPs of minor allele frequency <10% is used
#Import hapmap file 
data<-read.delim("/input/snp_i_N30_MAF10.hmp.txt")
#trim hapmap
data<-data[,-5:-11]
rownames(data)<-data[,1]
data<-data[,-1:-2]

#revert to mbp
data[,2]<-data[,2]/1000000

#converd to matrix for faster substitutions
dat1<-as.matrix(data)

#substitue characters
dat1<-sub("A", "A/A", dat1)
dat1<-sub("C", "C/C", dat1)
dat1<-sub("G", "G/G", dat1)
dat1<-sub("T", "T/T", dat1)
dat1<-sub("R", "A/G", dat1)
dat1<-sub("Y", "C/T", dat1)
dat1<-sub("S", "G/C", dat1)
dat1<-sub("W", "A/T", dat1)
dat1<-sub("K", "G/T", dat1)
dat1<-sub("M", "A/C", dat1)
dat1<-sub("N", "<NA>", dat1)

data<-data.frame(dat1) #revert to data.frame

data<-t(data) #transposing

gen<-makeGenotypes(data) #Make genotypes. 
save(gen, file="LDgenotypes.Rdata")

#remove snp position
snpspos<-data.frame(t(data[1:2,]))
snpspos[,1]<-as.numeric(as.character(snpspos[,1]))
snpspos[,2]<-as.numeric(as.character(snpspos[,2]))
save(snpspos, file="snpspos.Rdata")

#create heatmaps of pairwise LD for all SNPs, iterating per chromosome. Long process to run on a personal computer
snpchr<-split(snpspos, snpspos[,1])
for (c in length(snpchr):1){
  print(c)
  tmp<-gen[,colnames(gen) %in% rownames(snpchr[[c]])] #extract from the data frame the SNP data from each chromosome
  #remove all non-genotype columns for now
  zz<-lapply(tmp, class)
  notgen<-grep("character", zz)
  if(length(notgen)>0){
	tmp<-tmp[,-notgen]
	}
  #remove also from snpchr
  snpchr[[c]]<-snpchr[[c]][rownames(snpchr[[c]]) %in% colnames(tmp),]
  print(paste("SNP removed", length(notgen), "out of", ncol(tmp)))
  pdf(paste("CHR", c, "LDheatmap.pdf", sep = "."), width = 9, height = 7)
	tmp<-LDheatmap(tmp, genetic.distances=snpchr[[c]][,2], distance="physical",color=heat.colors(20))
  dev.off()
  save(tmp, file=paste("CHR", c, "LDheatmap.Rdata", sep = "."))
}

######################
#restart from here
options(stringsAsFactors = F)

#print out easier to handle files
load("snpspos.RData")
for (i in 1:10){
  print(i)
  load(paste("CHR", i, "LDheatmap.Rdata", sep = "."))
  tmpld<-as.matrix(tmp$LDmatrix)

  #keep positions and rearrange the dataframe / create a position data frame
  tmpos<-snpspos[rownames(snpspos) %in% colnames(tmpld),]
  tmpos$marker<-rownames(tmpos)
  colnames(tmpos)<-c("chr","position","marker")
  
  #collapse the LD matrix
  m<-cbind(which(!is.na(tmpld),arr.ind = TRUE),na.omit(as.vector(tmpld)))
  colnames(m)<-c("m1","m2","r2")
  tmpos<-tmpos[order(tmpos[,2]),]
  tmpos$idx<-1:nrow(tmpos)
  
  #add position information. For each pairwise comparison this will yield marker names, positions and r2
  newpos<-merge(m, tmpos, by.x="m1", by.y="idx", all.x=T)
  newpos1<-merge(newpos, tmpos, by.x="m2", by.y="idx", all.x=T)
  #drop chromosome data
  newpos1<-subset(newpos1, select=-c(chr.x,chr.y))
  
  
  dist<-abs(newpos1$position.x-newpos1$position.y) #calculate distance as the difference in Mbp between two markers
  newpos1<-cbind(newpos1,dist) #add distance vector to matrix
  #rearrange columns in matrix
  m<-newpos1[,c("m1","position.x","m2","position.y","r2","dist")]
  #order based on dist
  m<-m[order(m[,"dist"]),]
  save(m, file=paste("collapsed_matrix_LD_chr",i,"Rdata", sep="."))

}

#set interpolation function based on Hill and Weir (1988) equation based on Marroni (2011) script available at https://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/
Er<-function(C_,d){
  length(d)
  res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2))/(n*(2+C_*d)*(11+C_*d)))
  return(res)
}

#set additional parameters
n=260
exp <- expression(italic(paste(displaystyle(r^2))))

flist<-list()
datashd<-list()

#obtain fpoints by chromosomes
for (j in 1:10){
  print(j)
  load(paste("collapsed_matrix_LD_chr",j,"Rdata", sep="."))
  m1<-m
  m1<-m1[order(m1$dist),]
  #drop variables that will not be used
  m1<-m1[,c("r2","dist")]
  ld<-m1[,"r2"]
  d<-m1[,"dist"]
  nlm<-nls(ld~Er(C_,d[order(d)]),start=c(C_=0.1),control=nls.control(maxiter=100))
  C_<-summary(nlm)$parameters[1]
  fpoints<-Er(C_,d[order(d)])
  tmp<-data.frame(d,fpoints,j)
  flist[[j]]<-data.frame(d,fpoints,j)
  ldd<-0.1
  xpos<-tmp[which((tmp[,2]-ldd)<=0)[1],1]
  datashd[[j]]<-c(max(tmp[,2]), xpos)
}

fpointsall<-do.call(rbind, flist)
fpointsall[,3]<-as.factor(fpointsall[,3])
datashd<-do.call(rbind, datashd)

colnames(fpointsall)<-c("distance","r2","Chromosome")

fpointsall$distanceKb<-fpointsall$distance*1000

#Export maximum LD and LD decay distance per chromosome to .csv file
colnames(datashd)<-c("max", "Mb")
datashd<- as.data.frame(datashd)
datashd$Kb = datashd$Mb * 1000
write.csv(datashd,"SupplementaryInformation12.csv")

#now upgrade the basic plots using sliding windows

#set parameters:
#pairwise r2 is averaged for all surrounding markers within a +/-'win' of each SNP. In the paper, the 'win' parameter is set to 10 Mbp.
win<-10 
#resulting LD are plotted against physical position, averaging values over a sliding window of 'swin'% of each chromosome's markers. In the paper, 'swin' parameter is set to 5%.
swin<-floor(length(unique(m[,1]))*.05) 

#import centromeres' positions
plotter<-list()
centromeres<-read.csv("centromeres.csv")
centromeres$chromosome<-as.factor(centromeres$chromosome)
colnames(centromeres)[3]<-"Chromosome"

#estimate average LD per SNP
for (i in 1:10){
  load(paste("collapsed_matrix_LD_chr",i,"Rdata", sep="."))
  #remove long distances (the interval is on each side of the current SNP)
  m1<-m[which(m$dist<win),] 
  #sort by pos
  m1<-m1[order(m1$position.x),]
  #binding x and y positions
  m2<-m1[,-c(1:2)]
  colnames(m2)[1:2]<-c("m","position")
  m1<-m1[,-c(3:4)]
  colnames(m1)[1:2]<-c("m","position")
  m1<-rbind(m1,m2)
  #get mean for each unique SNP
  m1<-aggregate(m1['r2'], by=m1['position'], mean) #mean of r2 by position for SNPs
  #go through a sliding window
  #noting SNP positions
  md<-rollapply(m1$r2,width=swin, function(x) mean(x), partial=T) #mean of r2 by position for SNPs that are close by
  #create df for plotting
  md<-cbind(m1$position, md, i)
  plot(md[,1],md[,2], type="l", xlab="Mb", ylim=c(0,0.5))
  points(x=centromeres[i,1], y=0, col="red")
  plotter[[i]]<-md
}

#organize data frame
plotter<-as.data.frame(do.call(rbind, plotter))
plotter[,3]<-as.factor(plotter[,3])
colnames(plotter)<-c("position","mean_r2","Chromosome")

#Define graph parameters to represent each chromosome
cols<-brewer.pal(10,"Paired")
datashd$cols<-cols

#Export to .pdf file
pdf("Fig4_LD_main.pdf", width = 4, height = 3.5)
print(ggplot(fpointsall) + 
        geom_line(aes(y=r2, x=distanceKb, col=Chromosome),size=1) + 
        scale_color_manual(values=cols) +
        xlab("Kb") + ylab(expression(italic(r)*{}^2)) +
        coord_cartesian(ylim=c(0,0.1),xlim=c(0,100)) + #set x and y axes limit
        theme_classic(base_size = 14)+
        theme(legend.position="none")
)
dev.off()

#Export to .pdf file
pdf("Fig4_LD_insert.pdf", width = 8, height = 7)
print(ggplot(plotter, aes(x=position)) + 
        geom_line(aes(y=mean_r2, col=Chromosome),size=0.7) +
        geom_point(data=centromeres, mapping=aes(x=x.position, y=y.position),shape=17, color="black", size=1.5)+
        xlab("Mb") + ylab(expression(italic(r)*{}^2)) +
        scale_color_manual(values=cols) +
        theme_classic(base_size = 12) +
        theme(legend.justification=c(1,0), legend.position=c(1,0.4),strip.text.y = element_blank()) +
        scale_x_continuous(breaks=seq(0,320,20), sec.axis = dup_axis()) +
        facet_grid(Chromosome ~ .))
dev.off()